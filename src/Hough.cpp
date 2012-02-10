#include "Hough.hpp"

namespace sonar_wall_hough
{

Hough::Hough(const Config& config)
  : config(config)
  , houghspace(config)
  , startAngle(-1.0)
  , scanDirection(0)
  , lastAngle(-1.0)
  , filter(5, config.filterThreshold)
  , lastSpatialResolution(0.0)
{
  accMax = 10;
  double accAngle10 = M_PI / 4; //45°
  int accD10 = 16;
  accAsq = accAngle10*accAngle10/(-2* log(0.1));
  accDsq = accD10*accD10/(-2*log(0.1));
  accAngleIdx = houghspace.angle2Idx(-2*accAsq*log(1/accMax));
  accDIdx = houghspace.dst2Idx(-2*accDsq*log(1/accMax));
  
  localRangeDstIdx = houghspace.dst2Idx(10);
  localRangeAngleIdx = houghspace.angle2Idx(M_PI/8);
  
  angleTolerance = M_PI / 32;
}

Hough::~Hough()
{
}

void Hough::accumulate(std::vector<SonarPeak> peaks)
{
  int dstCenterIdx;
  int angleCenterIdx = houghspace.angle2Idx(peaks[0].alpha.rad); //all peaks have the same angle
  double angle = peaks[0].alpha.rad;
  boost::uint8_t* ptr;
  boost::uint16_t newValue;
  
  for(int angleIdx = angleCenterIdx-accAngleIdx; angleIdx < angleCenterIdx+accAngleIdx; angleIdx++)
  {
    for(std::vector<SonarPeak>::iterator it = peaks.begin(); it < peaks.end(); it++)
    {
      //calculate dst for this alpha and this peak
      dstCenterIdx = houghspace.dst2Idx((it->distance * cos(houghspace.idx2Angle(angleIdx) - angle)));
      //std::cout << "dstCenterIdx = " << dstCenterIdx << std::endl;
      
      for(int dstIdx = dstCenterIdx-accDIdx; dstIdx < dstCenterIdx+accDIdx; dstIdx++)
      {
	//get pointer to this houghspace bin
	if((ptr = houghspace.at(angleIdx, dstIdx)) != NULL)
	{
	  int accVal = accumulateValue(angle - houghspace.idx2Angle(angleIdx), houghspace.idx2Dst(dstCenterIdx) - houghspace.idx2Dst(dstIdx), *it);
	  //std::cout << "at angleIdx " << angleIdx << ", dstIdx " << dstIdx << ": " << accVal<< std::endl;
	  //check for overflow first
	  newValue = accVal + (*ptr);
	  //std::cout << "making " << (int)*ptr << "to " << (int)newValue << std::endl;
	  if(newValue > std::numeric_limits<boost::uint8_t>::max())
	    *ptr = std::numeric_limits<boost::uint8_t>::max();
	  else
	    *ptr = (boost::uint8_t)newValue;
	}
      }
    }
  }
}

void Hough::registerBeam(base::samples::SonarBeam beam)
{  
  //std::cout << "Last analysis angle = " << lastAnalysisAngle << std::endl;
  //std::cout << "registering beam with angle = " << beam.bearing << ".\n";
  std::vector<SonarPeak> peaks = filter.filter(beam, (int)(1.5/beam.getSpatialResolution())); //TODO 1.5 auslagern
  lastSpatialResolution = beam.getSpatialResolution();
  //append peaks to allPeaks
  allPeaks.insert(allPeaks.end(), peaks.begin(), peaks.end());
  
  //accumulate peaks to houghspace (all peaks have the same bearing)
  if(!peaks.empty())
    accumulate(peaks);
  
  //do we have a full 360° scan or a change of direction (ping-pong)?
  double actAngle = beam.bearing.getRad();
  if(actAngle < 0.0)
    actAngle += 2*M_PI;
  
  if(lastAngle <= -1.0)
  {
    startAngle = actAngle;
  }
  else
  {
    //what is the actual Direction?
    int actDir = (int)copysign(1.0, actAngle-lastAngle);
    if(fabs(actAngle-lastAngle) > M_PI)
      actDir *= -1;
    
    if(scanDirection == 0)
      scanDirection = actDir;
    else
    {
      //has direction cchanged?
      if(actDir != scanDirection)
      {
	scanDirection = actDir;
	startAngle = actAngle;
	analyzeHoughspace();
      }
      else if(startAngle == actAngle)
      {
	//we have a full 360° scan
	analyzeHoughspace();
      }
    }
  }  
  lastAngle = actAngle;
}

void Hough::analyzeHoughspace()
{
  std::cout << "analysis started" << std::endl;
  actualLines.clear();
  
  for(int dstIdx = -houghspace.getnumberOfPosDistances(); dstIdx <= houghspace.getnumberOfPosDistances(); dstIdx++)
  {
    for(int angleIdx = 0; angleIdx < houghspace.getNumberOfAngles(); angleIdx++)
    {
      if(*(houghspace.uncheckedAt(angleIdx, dstIdx)) > config.minLineVotesRatio*allPeaks.size() && isLocalMaximum(angleIdx, dstIdx))
      {
	//found valid line
	actualLines.push_back(Line(houghspace.idx2Angle(angleIdx), houghspace.idx2Dst(dstIdx), *(houghspace.uncheckedAt(angleIdx, dstIdx))));
	std::cout << "Found Line: alpha = " << actualLines.back().alpha << "[" << actualLines.back().alpha*180.0/M_PI << "°], d = " << actualLines.back().d << ", votes: " << actualLines.back().votes << std::endl;
	std::cout << "indices are angleIdx = " << angleIdx << ", dstIdx = " << dstIdx << std::endl;
      }
    }
  }
  std::cout << "found " << actualLines.size() << " Lines" << std::endl;
  postprocessLines();
  this->clear();
}

int Hough::accumulateValue(double deltaAngle, int deltaDst, SonarPeak& peak)
{
  return 0.01 * peak.strength * accMax * exp(-0.5 * deltaAngle*deltaAngle/accAsq) * exp(-0.5 * deltaDst*deltaDst/accDsq);
}

Houghspace* Hough::getHoughspace()
{
  return &houghspace;
}

std::vector< SonarPeak >* Hough::getAllPeaks()
{
  return &allPeaks;
}

std::vector< Line >* Hough::getActualLines()
{
  return &actualLines;
}

void Hough::clear()
{
  //clear Houghspace
  houghspace.clear();
  //clear allPeaks
  allPeaks.clear();
}

bool Hough::isLocalMaximum(int angleIdx, int dstIdx)
{
  boost::uint8_t thisValue = *(houghspace.at(angleIdx, dstIdx));
  for(int angleIdxCounter = angleIdx - localRangeAngleIdx; angleIdxCounter <= angleIdx + localRangeAngleIdx; angleIdxCounter++)
  {
    for(int dstIdxCounter = dstIdx - localRangeDstIdx; dstIdxCounter <= dstIdx + localRangeDstIdx; dstIdxCounter++)
    {
      //can be outside houghspace
      if(houghspace.at(angleIdxCounter, dstIdxCounter) == NULL)
	continue;
      
      if(*(houghspace.at(angleIdxCounter, dstIdxCounter)) > thisValue)
      {
	return false;
      }
      else if(*(houghspace.at(angleIdxCounter, dstIdxCounter)) == thisValue)
      {
	if(dstIdxCounter < dstIdx || (dstIdxCounter == dstIdx && angleIdxCounter < angleIdx))
	  return false;
      }
    }
  }
  return true;
}

void Hough::postprocessLines()
{
  std::vector<int> validDistances;
  validDistances.push_back(520);
  
  actualLines = Line::selectLines(actualLines, validDistances);
}


}//end namespace