#include "Hough.hpp"

namespace sonar_wall_hough
{

Hough::Hough(const Config& config)
  : config(config)
  , houghspace(config)
  , startAngle(-1.0)
  , halfDone(false)
  , scanDirection(0)
  , lastAngle(-1.0)
  , filter(5, config.filterThreshold, config.withMinimumFilter)
  , lastSpatialResolution(0.0)
  , lastOrientation(base::Angle::fromRad(0.0))
  , firstOrientation(base::Angle::fromRad(0.0))
  , actualPosition(std::pair<double,double>(0.0,0.0))
  , accMax(config.gain)
  , orientationDrift(0.0)
  , orientationDriftDetected(false)
{
  //accMax = 5;
  double accAngle10 = M_PI * 3/8; //67.5°
  int accD10 = 12;
  accAsq = accAngle10*accAngle10/(-2* log(0.1));
  accDsq = accD10*accD10/(-2*log(0.1));
  accAngleIdx = houghspace.angle2Idx(-2*accAsq*log(1/accMax));
  accDIdx = houghspace.dst2Idx(-2*accDsq*log(1/accMax));
  //std::cout << "accDIdx = " << accDIdx << std::endl;
  
  localRangeDstIdx = houghspace.dst2Idx(10);
  localRangeAngleIdx = houghspace.angle2Idx(M_PI/8);
  
  lastPosition=std::make_pair(0.0,0.0);  
  
  newDataCount=0;
  minError=-1;
  maxError=-1;
  sumError=0;
  count=0;
  
  houghspace.clear();
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

void Hough::accumulate(SonarPeak peak)
{
  int dstCenterIdx;
  int angleCenterIdx = houghspace.angle2Idx(peak.alpha.rad); //all peaks have the same angle
  double angle = peak.alpha.rad;
  boost::uint8_t* ptr;
  boost::uint16_t newValue;
  double downvoteNear;
  
  for(int angleIdx = angleCenterIdx-accAngleIdx; angleIdx < angleCenterIdx+accAngleIdx; angleIdx++)
  {

      //calculate dst for this alpha and this peak
      dstCenterIdx = houghspace.dst2Idx((peak.distance * cos(houghspace.idx2Angle(angleIdx) - angle)));
      
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
	  int accVal = downvoteNear * accumulateValue(angle - houghspace.idx2Angle(angleIdx), houghspace.idx2Dst(dstCenterIdx) - houghspace.idx2Dst(dstIdx), peak);
	  
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

void Hough::registerBeam(base::samples::SonarBeam beam)
{
  //std::cout << "Last analysis angle = " << lastAnalysisAngle << std::endl;
  //std::cout << "registering beam with angle = " << beam.bearing << ".\n";
  //std::cout << "registering beam with " << beam.beam.size() << " buckets" << std::endl;
  
  //correct angle for this beam
  double actAngle = beam.bearing.getRad();
  beam.bearing = beam.bearing + (lastOrientation - firstOrientation);
  //base::Angle bla = base::Angle::fromRad(angleCorrect);
  //std::cout << "angle is " << beam.bearing << ", should be " << bla << std::endl;
  
  std::vector<SonarPeak> peaks = filter.filter(beam, (int)(config.minDistance/beam.getSpatialResolution()));
  //std::cout << "registering beam with " << peaks.size() << " peaks" << std::endl;
  lastSpatialResolution = beam.getSpatialResolution();

  
  //accumulate peaks to houghspace (all peaks have the same bearing)
  if(!peaks.empty()){
    //std::cout << "PEAKS: " << peaks.size() << std::endl;
    
    //for(std::vector<SonarPeak>::iterator it = peaks.begin(); it < peaks.end(); it++)
    //	std::cout << "Distance " <<it->distance << " Str "<< int(it->strength)<< std::endl; 
       
    if(config.poseCorrection){
          
	if(config.correctToFirstPosition){
	  correctPeaks(&peaks); 	  
	  
	}else{
	  posPeaks.push_back(std::make_pair(actualPosition,peaks));
	  newDataCount++;
	  allPeaks.clear();
	}  
 
    }else{  
      //append peaks to allPeaks
      allPeaks.insert(allPeaks.end(), peaks.begin(), peaks.end());
      
      accumulate(peaks);
    }
    
  }  
  
  //do we have a full 360° scan or a change of direction (ping-pong)?
  if(actAngle < 0.0)
    actAngle += 2*M_PI;
  
  //std::cout << "last = " << lastAngle << ", act = " << actAngle << ", start = " << startAngle << std::endl;
  if(lastAngle <= -1.0)
  {
    startAngle = actAngle;
  }
  else
  {
    int pass = 0;
    //what is the actual Direction?
    int actDir = (int)copysign(1.0, actAngle-lastAngle);
    if(fabs(actAngle-lastAngle) > M_PI)
    {
      actDir *= -1;
      pass = 1;
    }
    
    if(scanDirection == 0)
      scanDirection = actDir;
    else
    {
      //has direction changed?
      if(actDir != scanDirection)
      {
	scanDirection = actDir;
	startAngle = actAngle;
	analyzeHoughspace();
      }
      else
      {
	//adjust actAngle
	double actAdjusted = actAngle;
	if(actDir == +1 && actAdjusted < startAngle)
	  actAdjusted += 2*M_PI;
	else if(actDir == -1 && actAdjusted > startAngle)
	  actAdjusted -= 2*M_PI;
	
	if(halfDone)
	{
	  if((actAdjusted - startAngle) * actDir < M_PI)
	  {
	    //we have a full 360° scan
	    startAngle = actAngle;
	    halfDone = false;
	    
	    if(config.poseCorrection && !config.correctToFirstPosition){
	      
	      correctPeaks(posPeaks);	      
	    }
	    
	    
	    analyzeHoughspace();
	  }
	}
	else
	{
	  if((actAdjusted - startAngle) * actDir > M_PI)
	    halfDone = true;
	}
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
  if(!config.poseCorrection || config.correctToFirstPosition) allPeaks.clear();
  
  std::cout << "NewDataCount:" << newDataCount << std::endl;
  posPeaks.erase(posPeaks.begin(),posPeaks.end()-newDataCount);
  //posPeaks.clear();
  //posPeaks.erase(posPeaks.begin(),posPeaks.end());
  
  firstOrientation = lastOrientation;
    
  firstPosition = std::make_pair(lastPosition.first,lastPosition.second);
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
  
  //actualLines = Line::selectLines(actualLines, basinSize, lastSpatialResolution, (config.sensorAngularResolution*4/180*M_PI), firstOrientation.getRad(), true, true);
  double* actualOrientationDrift = NULL;
  actualLines = Line::selectLines2(actualLines, basinSize, lastSpatialResolution, (config.sensorAngularResolution*8/180*M_PI), -firstOrientation.getRad(), houghspace, actualOrientationDrift); //orientation must be negative (heading of auv = - heading of basin in respect of auv)
  if(actualLines.size() < 4)
    return;
  
  if(actualOrientationDrift != NULL)
  {
    orientationDriftDetected = true;
    orientationDrift = *actualOrientationDrift;
    delete actualOrientationDrift;
  }
  
  std::cout << "the lines are:" <<std::endl << "horz:" <<std::endl;
  std::cout << "Line: alpha = " << actualLines[0].alpha << "[" << actualLines[0].alpha*180.0/M_PI << "°], d = " << actualLines[0].d << ", votes: " << actualLines[0].votes << std::endl;
  std::cout << "Line: alpha = " << actualLines[1].alpha << "[" << actualLines[1].alpha*180.0/M_PI << "°], d = " << actualLines[1].d << ", votes: " << actualLines[1].votes << std::endl;
  std::cout << "vert:\nLine: alpha = " << actualLines[2].alpha << "[" << actualLines[2].alpha*180.0/M_PI << "°], d = " << actualLines[2].d << ", votes: " << actualLines[2].votes << std::endl;
  std::cout << "Line: alpha = " << actualLines[3].alpha << "[" << actualLines[3].alpha*180.0/M_PI << "°], d = " << actualLines[3].d << ", votes: " << actualLines[3].votes << std::endl;
  
  calculateError();  
  
  //calculate position
  double xScale = fabs(actualLines[2].d-actualLines[3].d) / config.basinWidth;
  double yScale = fabs(actualLines[0].d-actualLines[1].d) / config.basinHeight;
  double xPos = -(actualLines[2].d + actualLines[3].d) / 2 / xScale;
  double yPos = -(actualLines[0].d + actualLines[1].d) / 2 / yScale;
  std::cout << "AUV at x = " << xPos << ", y = " << yPos << std::endl;
  
  actualPosition.first = xPos;
  actualPosition.second = yPos;
  
  /*
  //x- and y-axis for debugging
  actualLines.push_back(Line(actualLines[2].alpha,-xPos*xScale,0));
  actualLines.push_back(Line(actualLines[0].alpha,-yPos*yScale,0));
  */  
}

void Hough::calculateError()
{
  //mean square error of peaks & line supporters
  std::vector<int> lineSupporters(actualLines.size(), 0);
  meanSqErr = 0.0;
  int nrOfPeaks = allPeaks.size();
  for(int i = 0; i < nrOfPeaks; i++)
  {
    std::pair<int,double> lineDst = allPeaks[i].smallestDstFromLine(actualLines);
    //std::cout << ", " << lineDst.second;
    meanSqErr += (lineDst.second*lineDst.second)/nrOfPeaks;
    if(lineDst.second <= 20) //TODO: property if needed
    {
      lineSupporters[lineDst.first]++;
    }
  }
  //std::cout << std::endl;
  
  //distance of parallel lines
  basinWidthDiff = fabs(actualLines[2].d-actualLines[3].d) * lastSpatialResolution - config.basinWidth;
  basinHeightDiff = fabs(actualLines[0].d-actualLines[1].d) * lastSpatialResolution - config.basinHeight;
  
  supportRatio = (double)(lineSupporters[0]+lineSupporters[1]+lineSupporters[2]+lineSupporters[3])/nrOfPeaks;
  //cout the values
  std::cout << "QUALITY OF MEASUREMENT:" << std::endl;
  std::cout << "difference of parallel walls is " << basinWidthDiff << "m at x and " << basinHeightDiff << "m at y." << std::endl;
  std::cout << "Mean square error of peaks is " << meanSqErr << std::endl;
  std::cout << "the wall's supporters: " << lineSupporters[0] << ", " << lineSupporters[1] << ", " << lineSupporters[2] << ", " <<lineSupporters[3] << std::endl;
 
}

void Hough::setOrientation(double orientation)
{
  
  this->lastOrientation = base::Angle::fromRad(orientation + (config.angleDelta*M_PI/180.0));
  //std::cout << "orientation set to " << lastOrientation << std::endl;
}

base::Angle Hough::getOrientation()
{
  return firstOrientation;
}

std::pair< double, double > Hough::getActualPosition()
{
  return actualPosition;
}

double Hough::getOrientationDrift()
{
  if(orientationDriftDetected)
  {
    return orientationDrift;
  }
  else
    return 0.0;
}

double Hough::getBasinHeightDiff()
{
  return basinHeightDiff;
}
double Hough::getBasinWidthDiff()
{
  return basinWidthDiff;
}
double Hough::getMeanSqErr()
{
  return meanSqErr;
}
double Hough::getSupportRatio()
{
  return supportRatio;
}

void Hough::setPosition(std::pair<double,double> pose){
  lastPosition=pose;
}


void Hough::correctPeaks(std::vector<SonarPeak>* peaks){
      
      //Relative movement of avalon between begin of sonar-scan and now
      double xDiv= cos(firstOrientation.getRad()) * (lastPosition.first - firstPosition.first) 
		    -sin(firstOrientation.getRad()) * (lastPosition.second - firstPosition.second);
      double yDiv= cos(firstOrientation.getRad()) * (lastPosition.second - firstPosition.second)
		    +sin(firstOrientation.getRad()) * (lastPosition.first - firstPosition.first);
      
      double sonarXdiv = cos(lastOrientation.getRad() - firstOrientation.getRad()) * config.avalonSonarPose; 
      double sonarYdiv = sin(lastOrientation.getRad() - firstOrientation.getRad()) * config.avalonSonarPose;
		    
      for(std::vector<SonarPeak>::iterator it = peaks->begin(); it < peaks->end(); it++){
	 
	//Relative Position of the Peak to the actual Position of Avalon
	 double x= cos(it->alpha.getRad()) * it->distance + sonarXdiv; 
	 double y= sin(it->alpha.getRad()) * it->distance + sonarYdiv;
	 	 
	 //calculate new angle and distance of peak
	 double angle = atan2(yDiv+y, xDiv+x);
	 double distance = sqrt(pow((yDiv+y),2.0) + pow((xDiv+x),2.0));
	 
	 it->distance=distance;
	 it->alpha.rad=angle;
	 
	 accumulate(*it);
      }
      
      //append peaks to allPeaks
      allPeaks.insert(allPeaks.end(), peaks->begin(), peaks->end());
}

void Hough::correctPeaks(std::vector <std::pair<std::pair<double,double>,std::vector<SonarPeak> > > peaks){
    std::pair<double,double> lastPos=peaks.back().first;
   
    newDataCount=0;
        
    for(int i=0; i<peaks.size() ; i++){
      std::pair<double,double> accPos=peaks[i].first;
      std::vector<SonarPeak> accPeaks = peaks[i].second;
      
      //Relative Movement of avalon between the time of the sonar_peak and the finish of the sonar-scan
      double xDiv= cos(firstOrientation.getRad()) * (accPos.first - lastPos.first) 
		    -sin(firstOrientation.getRad()) * (accPos.second - lastPos.second);
      double yDiv= cos(firstOrientation.getRad()) * (accPos.second - lastPos.second)
		    +sin(firstOrientation.getRad()) * (accPos.first - lastPos.first);
      
      for(std::vector<SonarPeak>::iterator it = accPeaks.begin(); it< accPeaks.end() ; it++){
	
	//Calculate new angle and distance
	double x=cos(it->alpha.getRad()) * it->distance;
	double y=sin(it->alpha.getRad()) * it->distance;	
	
	it->distance=sqrt(pow(xDiv+x,2.0) + pow(yDiv+y,2.0));
	it->alpha.rad=atan2(yDiv+y,xDiv+x);	
		
	accumulate(*it);	
      }
      
      //apend corrected peaks to allpeaks
      allPeaks.insert(allPeaks.end(),accPeaks.begin(),accPeaks.end());
    }
    
}

void Hough::calcPositionError(){
    double error= sqrt(pow(lastPosition.first-actualPosition.first,2.0)+pow(lastPosition.second-actualPosition.second,2.0));
    
    if(error < minError || minError==-1){
	minError=error;
    }
    if(error > maxError || maxError==-1){
	maxError=error;
    }  
    sumError+=error;
    count++;
} 

double Hough::getMinError(){
  return minError;
}

double Hough::getMaxError(){
  return maxError;
}

double Hough::getAvgError(){
  return sumError/count;
}  

}//end namespace