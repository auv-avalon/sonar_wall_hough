#include "FilterHelper.hpp"

namespace sonar_wall_hough
{
  MinWindow::MinWindow(int kernelSize)
    :kernelSize(kernelSize)
    ,minDeque(kernelSize+1)
    ,idxDeque(kernelSize+1)
    ,lastIdx(0)
  {
  }

  MinWindow::~MinWindow()
  {
  }

  uint8_t MinWindow::pushValue()
  {
    //push maximum possible value
    return pushValue(std::numeric_limits<uint8_t>::max());
  }
  
  uint8_t MinWindow::pushValue(uint8_t value)
  {
    //delete all elements being greater than value
    while(!minDeque.empty() && minDeque.back() > value)
    {
      minDeque.pop_back();
      idxDeque.pop_back();
    }
    //insert new value at end
    minDeque.push_back(value);
    idxDeque.push_back(lastIdx);
    
    //pop front if idx is out of kernel
    if(lastIdx - idxDeque.front() >= kernelSize)
    {
      minDeque.pop_front();
      idxDeque.pop_front();
    }
    
    //increase lastIdx for next run
    lastIdx++;
    
    //return minimum which is front of minDeque
    return minDeque.front();
  }
  
  void MinWindow::clear()
  {
    minDeque.clear();
    idxDeque.clear();
    lastIdx = 0;
  }
  
  BeamFilterDst::BeamFilterDst(FilterType type, int kernelSize)
    :type(type)
    ,kernelSize(kernelSize)
    ,minWindow(0)
  {
    if(type == minimum)
      minWindow = new MinWindow(kernelSize);
  }

  BeamFilterDst::~BeamFilterDst()
  {
    if(minWindow != NULL)
      delete minWindow;
  }
    
  base::samples::SonarBeam BeamFilterDst::filter(base::samples::SonarBeam sonarBeam)
  {
    //std::cout << "BeamFilterDst with type = " << type << std::endl;
    base::samples::SonarBeam filteredBeam(sonarBeam);
    if(type == minimum)
      minWindow->clear();
    
    //walk through beam and perform filtering
    int temp;
    std::vector<uint8_t>::iterator in = sonarBeam.beam.begin();
    std::vector<uint8_t>::iterator out = filteredBeam.beam.begin();
    for(; in < sonarBeam.beam.end()-kernelSize/2; in++, out++)
    {
      if(type == minimum)
      {
	//first fill the minWindow
	if(in - sonarBeam.beam.begin() < kernelSize/2)
	{
	  minWindow->pushValue(*in);
	  //do not proceed with out
	  out--;
	}
	else
	{
	  *out = minWindow->pushValue(*in);
	}
      }
      else if(type == sobelGaussDst)
      {
	//write 0 if distance to 0 or beam.size is smaller than half kernelSize
	if((in - sonarBeam.beam.begin()) < kernelSize/2 || (sonarBeam.beam.end() - in) < kernelSize/2)
	  *out = 0;
	else
	{
	  temp = abs(0.25*in[2] + 0.5*in[1] - 0.5*in[-1] - 0.25*in[-2]); //TODO automatic kernel generation
	  *out = temp > 255 ? 255 : temp;
	}
      }
      else if(type == sobelGaussPhi)
      {
	//write 0 if distance to 0 or beam.size is smaller than half kernelSize
	if((in - sonarBeam.beam.begin()) < kernelSize/2 || (sonarBeam.beam.end() - in) < kernelSize/2)
	  *out = 0;
	else
	{
	  temp = abs(0.25*in[2] + in[1] + 1.5*in[0] + in[-1] + 0.25*in[-2]); //TODO automatic kernel generation
	  *out = temp > 255 ? 255 : temp;
	}
      }
    }
    
    //fill rest of out if in minimum type
    /*
    if(type == minimum)
    {
      int i = 0;
      while(out < filteredBeam.beam.end())
      {
	*out = minWindow->pushValue();
	out++;
      }
    }
    */
    return filteredBeam;
  }

  BeamFilterPhi::BeamFilterPhi(FilterType type, int kernelSize)
    :type(type)
    ,kernelSize(kernelSize)
    ,lastbeams(kernelSize)
    ,minWindow(0)
  {
    if(type == minimum)
      minWindow = new MinWindow(kernelSize);
  }

  BeamFilterPhi::~BeamFilterPhi()
  {
    if(minWindow != NULL)
      delete minWindow;
  }

  base::samples::SonarBeam BeamFilterPhi::filter(base::samples::SonarBeam sonarBeam)
  {
    //std::cout << "BeamFilterPhi with type = " << type << std::endl;
    //organize last beams
    //clear deque if beam size has changed
    if(lastbeams.size() > 0 && lastbeams[0].beam.size() != sonarBeam.beam.size())
      lastbeams.clear();
    
    lastbeams.push_back(sonarBeam);
    if(lastbeams.size() > kernelSize)
      lastbeams.pop_front();
    
    //take middle beam as filtering center if enough beams collected
    int idx;
    if(lastbeams.size() == kernelSize)
      idx = kernelSize/2;
    else
      idx = 0;
    base::samples::SonarBeam filteredBeam(lastbeams.at(idx));
    
    //walk through beams and perform filtering
    int temp;
    std::vector<uint8_t>::iterator out = filteredBeam.beam.begin();
    for(int i = 0; out < filteredBeam.beam.end(); out++, i++)
    {
      if(type == minimum)
      {
	minWindow->clear();
	for(int j = 0; j < lastbeams.size() - 1; j++)
	  minWindow->pushValue(lastbeams[j].beam[i]);
	*out = minWindow->pushValue(lastbeams.back().beam[i]);
      }
      else if(type == sobelGaussDst)
      {
	//write 0 if not enough beams collected
	if(lastbeams.size() < kernelSize)
	  *out = 0;
	else
	{
	  temp = abs(0.25*lastbeams[0].beam[i] + lastbeams[1].beam[i] + 1.5*lastbeams[2].beam[i] + lastbeams[3].beam[i] + 0.25*lastbeams[4].beam[i]); //TODO automatic kernel generation
	  *out = temp > 255 ? 255 : temp;
	}
      }
      else if(type == sobelGaussPhi)
      {
	//write 0 if not enough beams collected
	if(lastbeams.size() < kernelSize)
	  *out = 0;
	else
	{	  
	  temp = abs(0.25*lastbeams[0].beam[i] + 0.5*lastbeams[1].beam[i] - 0.5*lastbeams[3].beam[i] - 0.25*lastbeams[4].beam[i]); //TODO automatic kernel generation
	  *out = temp > 255 ? 255 : temp;
	}
      }
    }
    return filteredBeam;
  }

  Filter::Filter(int kernelSize, uint8_t threshold, bool withMinimum)
    :threshold(threshold)
    ,withMinimum(withMinimum)
    ,filterDstMin(minimum, kernelSize)
    ,filterDstSGDst(sobelGaussDst, kernelSize)
    ,filterDstSGPhi(sobelGaussPhi, kernelSize)
    ,filterPhiMin(minimum, kernelSize)
    ,filterPhiSGDst(sobelGaussDst, kernelSize)
    ,filterPhiSGPhi(sobelGaussPhi, kernelSize)
  {
  }

  Filter::~Filter()
  {
  }

  std::vector< SonarPeak > Filter::filter(base::samples::SonarBeam sonarBeam, int minDistance)
  {
    std::vector<SonarPeak> peaks;
    if(sonarBeam.beam.size() < minDistance)
      return peaks;
    
    //first perform a minimum filtering
    base::samples::SonarBeam minFilteredBeam;
    if(withMinimum)
      minFilteredBeam = filterPhiMin.filter(filterDstMin.filter(sonarBeam));
    else
      minFilteredBeam = sonarBeam;
    
    //then perform a sobel + gaussian filtering
    //Phi Directed sobel
    base::samples::SonarBeam fullFilteredDst = filterPhiSGDst.filter(filterDstSGDst.filter(minFilteredBeam));
    //Dst Directed sobel
    base::samples::SonarBeam fullFilteredPhi = filterPhiSGPhi.filter(filterDstSGPhi.filter(minFilteredBeam));
    
    //apply a threshold and push into SonarPeak vector if vaule is local maximum
    uint8_t strength, strengthPre = 0, strengthPrePre = 0;
    for(int i = minDistance; i < fullFilteredDst.beam.size(); i++)
    {
      strength = sqrt(fullFilteredDst.beam[i] * fullFilteredDst.beam[i] + fullFilteredPhi.beam[i] * fullFilteredPhi.beam[i]);
      //std::cout << (int)strength << std::endl;
      if(strengthPre > threshold && strengthPre > strength && strengthPre >= strengthPrePre)
	peaks.push_back(SonarPeak(fullFilteredDst.bearing, i-1, strengthPre));
      
      //push strengths
      strengthPrePre = strengthPre;
      strengthPre = strength;
    }
    
    return peaks;
  }

}//ending namespace































