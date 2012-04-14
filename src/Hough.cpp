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
  , lastOrientation(base::Angle::fromRad(0.0))
  , firstOrientation(base::Angle::fromRad(0.0))
  , actualPosition(std::pair<double,double>(0.0,0.0))
{
  accMax = 10;
  double accAngle10 = M_PI * 3/8; //67.5°
  int accD10 = 12;
  accAsq = accAngle10*accAngle10/(-2* log(0.1));
  accDsq = accD10*accD10/(-2*log(0.1));
  accAngleIdx = houghspace.angle2Idx(-2*accAsq*log(1/accMax));
  accDIdx = houghspace.dst2Idx(-2*accDsq*log(1/accMax));
  //std::cout << "accDIdx = " << accDIdx << std::endl;
  
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
  double downvoteNear;
  
  for(int angleIdx = angleCenterIdx-accAngleIdx; angleIdx < angleCenterIdx+accAngleIdx; angleIdx++)
  {
    //std::cout << "angleIdx = " << angleIdx << std::endl;
    for(std::vector<SonarPeak>::iterator it = peaks.begin(); it < peaks.end(); it++)
    {
      //calculate dst for this alpha and this peak
      dstCenterIdx = houghspace.dst2Idx((it->distance * cos(houghspace.idx2Angle(angleIdx) - angle)));
      
      //std::cout << "angleIdx = " << angleIdx << ", dstCenterIdx = " << dstCenterIdx << std::endl;
      
      //calculate downvote for small peak distance
      //dstCenterIdx is between 0 and houghspace.getnumberOfPosDistances()
      double dstRatio = (double)dstCenterIdx / houghspace.getnumberOfPosDistances();
      //downvoteNear = dstRatio;
      downvoteNear = 1-exp(-0.25*dstRatio*dstRatio/0.0016);
      
      //std::cout << "dstCenterIdx = " << dstCenterIdx << std::endl;
      //std::cout << "downvoteNear = " << downvoteNear << std::endl;
      
      for(int dstIdx = dstCenterIdx-accDIdx; dstIdx < dstCenterIdx+accDIdx; dstIdx++)
      {
	//std::cout << "hough.cpp angleIdx " << angleIdx << ", dstIdx " << dstIdx << std::endl;
	//get pointer to this houghspace bin
	if((ptr = houghspace.at(angleIdx, dstIdx)) != NULL)
	{
	  int accVal = downvoteNear * accumulateValue(angle - houghspace.idx2Angle(angleIdx), houghspace.idx2Dst(dstCenterIdx) - houghspace.idx2Dst(dstIdx), *it);
	  
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
  
  //correct angle for this beam
  double actAngle = beam.bearing.getRad();
  beam.bearing = beam.bearing + (lastOrientation - firstOrientation);
  //base::Angle bla = base::Angle::fromRad(angleCorrect);
  //std::cout << "angle is " << beam.bearing << ", should be " << bla << std::endl;
  
  std::vector<SonarPeak> peaks = filter.filter(beam, (int)(config.minDistance/beam.getSpatialResolution()));
  lastSpatialResolution = beam.getSpatialResolution();
  //append peaks to allPeaks
  allPeaks.insert(allPeaks.end(), peaks.begin(), peaks.end());
  
  //accumulate peaks to houghspace (all peaks have the same bearing)
  if(!peaks.empty())
    accumulate(peaks);
  
  //do we have a full 360° scan or a change of direction (ping-pong)?
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
  
  firstOrientation = lastOrientation;
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
  std::pair<int,int> basinSize(config.basinHeight, config.basinWidth);
  actualLines = Line::selectLines(actualLines, basinSize, lastSpatialResolution, (config.sensorAngularResolution*4/180*M_PI), firstOrientation.getRad(), true, true);
  if(actualLines.size() < 4)
    return;
  
  //calculate position
  double xScale = fabs(actualLines[2].d-actualLines[3].d) / config.basinWidth;
  double yScale = fabs(actualLines[0].d-actualLines[1].d) / config.basinHeight;
  double xPos = -(actualLines[2].d + actualLines[3].d) / 2 / xScale;
  double yPos = -(actualLines[0].d + actualLines[1].d) / 2 / yScale;
  std::cout << "AUV at x = " << xPos << ", y = " << yPos << std::endl;
  std::cout << "Scale deltas are x: "<<(xScale*lastSpatialResolution)<<", y: " << (yScale*lastSpatialResolution) << std::endl;
  
  actualPosition.first = xPos;
  actualPosition.second = yPos;
  
  //x- and y-axis for debugging
  actualLines.push_back(Line(actualLines[2].alpha,-xPos*xScale,0));
  actualLines.push_back(Line(actualLines[0].alpha,-yPos*yScale,0));
}

void Hough::setOrientation(double orientation)
{
  
  this->lastOrientation = base::Angle::fromRad(orientation + (config.angleDelta*M_PI/180.0));
}

base::Angle Hough::getOrientation()
{
  return firstOrientation;
}

std::pair< double, double > Hough::getActualPosition()
{
  return actualPosition;
}



}//end namespace