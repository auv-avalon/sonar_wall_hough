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

void Hough::accumulate(SonarPeak peak)
{
  int dstCenterIdx;
  int angleCenterIdx = houghspace.angle2Idx(peak.alpha.rad);
  
  for(int angleIdx = angleCenterIdx-accAngleIdx; angleIdx < angleCenterIdx+accAngleIdx; angleIdx++)
  {
    //calculate dst for this alpha
    dstCenterIdx = houghspace.dst2Idx((peak.distance * cos(houghspace.idx2Angle(angleIdx) - peak.alpha.rad)));
    //std::cout << "dstCenterIdx = " << dstCenterIdx << std::endl;
    for(int dstIdx = dstCenterIdx-accDIdx; dstIdx < dstCenterIdx+accDIdx; dstIdx++)
    {
      //get pointer to this houghspace bin
      boost::uint8_t* ptr;
      if((ptr = houghspace.at(angleIdx, dstIdx)) != NULL)
      {
	int accVal = accumulateValue(peak.alpha.rad - houghspace.idx2Angle(angleIdx), houghspace.idx2Dst(dstCenterIdx) - houghspace.idx2Dst(dstIdx));
	//std::cout << "at angleIdx " << angleIdx << ", dstIdx " << dstIdx << ": " << accVal<< std::endl;
	//check for overflow first
	//std::cout << "adding " << 0.01 * accVal * peak.strength << std::endl;
	if((int)*ptr + 0.01 * accVal * peak.strength > std::numeric_limits<boost::uint8_t>::max())
	  *ptr = std::numeric_limits<boost::uint8_t>::max();
	else
	  *ptr += (int)(0.01 * accVal * peak.strength);
      }
    }
  }
}

void Hough::registerBeam(base::samples::SonarBeam beam)
{  
  //std::cout << "Last analysis angle = " << lastAnalysisAngle << std::endl;
  //std::cout << "registering beam with angle = " << beam.bearing << ".\n";
  std::vector<SonarPeak> peaks = filter.filter(beam, (int)(1.5/beam.getSpatialResolution()));
  //append peaks to allPeaks
  allPeaks.insert(allPeaks.end(), peaks.begin(), peaks.end());
  for(int i = 0; i < (int)peaks.size(); i++)
  {
    accumulate(peaks.at(i));
  }
  
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
      if(*(houghspace.uncheckedAt(angleIdx, dstIdx)) > config.minLineVotes && isLocalMaximum(angleIdx, dstIdx))
      {
	//found valid line
	actualLines.push_back(Line(houghspace.idx2Angle(angleIdx), houghspace.idx2Dst(dstIdx), *(houghspace.uncheckedAt(angleIdx, dstIdx))));
	std::cout << "Found Line: alpha = " << actualLines.back().alpha << "[" << actualLines.back().alpha*180.0/M_PI << "°], d = " << actualLines.back().d << ", votes: " << actualLines.back().votes << std::endl;
	std::cout << "indices are angleIdx = " << angleIdx << ", dstIdx = " << dstIdx << std::endl;
      }
    }
  }
  std::cout << "found " << actualLines.size() << " Lines" << std::endl;
    
  this->clear();
}

int Hough::accumulateValue(double deltaAngle, int deltaDst)
{
  return accMax * exp(-0.5 * deltaAngle*deltaAngle/accAsq) * exp(-0.5 * deltaDst*deltaDst/accDsq);
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

void Hough::postprocessLinesPool(std::vector< Line >& lines, double poolArea)
{/*
  //find all parallel and opposite lying line pairs
  std::vector<std::pair<int, int> > pairs;
  for(int i = 0; i < (int)lines.size(); i++)
  {
    for(int j = i+1; j < (int)lines.size(); j++)
    {
      //are lines i and j parallel and opposite?
      if((fabs(lines[i].alpha-lines[j].alpha) < angleTolerance && copysign(1.0,lines[i].d) != copysign(1.0,lines[j].d)) ||
	(fabs(fabs(lines[i].alpha-lines[j].alpha)-M_PI) < angleTolerance && copysign(1.0,lines[i].d) == copysign(1.0,lines[j].d))     )
      {
	pairs.push_back(std::pair<int, int>(i, j));
      }
    }
  }
  std::cout << "found " << pairs.size() << " pairs." << std::endl;
  
  //find all corresponding pairs
  std::vector<std::pair<std::pair<int, int>, std::pair<int, int> > > quads;
  for(int i = 0; i < (int)pairs.size(); i++)
  {
    for(int j = i+1; j < (int)pairs.size(); j++)
    {
      //TODO: Correct mean angles (i.e. 175° and 5°)
      double angleI = (lines[pairs[i].first].alpha + lines[pairs[i].second].alpha) /2;
      double angleJ = (lines[pairs[j].first].alpha + lines[pairs[j].second].alpha) /2;
      if(fabs(fmod(fabs(angleI-angleJ), M_PI) - M_PI/2) < angleTolerance)
	quads.push_back(std::pair<std::pair<int, int>, std::pair<int, int> >(pairs[i], pairs[j]));
    }
  }
  std::cout << "found " << quads.size() << " quads." << std::endl;
  
  double minDeltatoPoolSize = INFINITY;
  int bestI = -1;
  for(int i = 0; i < (int)quads.size(); i++)
  {
    //calculate area of rectangle
    int firstLength = lines[quads[i].first.first].d - lines[quads[i].first.second].d;
    int secondLength = lines[quads[i].second.first].d - lines[quads[i].second.second].d;
    double area = firstLength*spatialResolution * secondLength*spatialResolution;
    if(fabs(area - poolArea) < minDeltatoPoolSize)
    {
      minDeltatoPoolSize = fabs(area - poolArea);
      bestI = i;
    }
  }
  std::vector<Line> newLines;
  if(bestI > -1)
  {
    newLines.push_back(lines[quads[bestI].first.first]);
    newLines.push_back(lines[quads[bestI].first.second]);
    newLines.push_back(lines[quads[bestI].second.first]);
    newLines.push_back(lines[quads[bestI].second.second]);
  }
  lines = newLines;*/
}


}//end namespace

/*beckenfläche
winkelschritt berechnen vs parameter
tastendruck als ping-pong trigger ?*/