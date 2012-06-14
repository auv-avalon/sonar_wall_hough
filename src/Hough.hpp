#ifndef _SONAR_WALL_HOUGH_HOUGH_HPP_
#define _SONAR_WALL_HOUGH_HOUGH_HPP_

#include <base/samples/sonar_beam.h>
#include <base/samples/frame.h>
#include <math.h>
#include <limits.h>
#include <algorithm>
#include "Config.hpp"
#include "Line.hpp"
#include "Houghspace.hpp"
#include "SonarPeak.hpp"
#include "FilterHelper.hpp"

namespace sonar_wall_hough 
{

class Hough
{

public:
  Hough(const Config& config = Config());
  ~Hough();
  void accumulate(std::vector<SonarPeak> peaks);
  void registerBeam(base::samples::SonarBeam beam);
  void analyzeHoughspace();
  std::vector<SonarPeak>* getAllPeaks();
  std::vector<Line>* getActualLines();
  std::pair<double,double> getActualPosition();
  base::Angle getOrientation();
  Houghspace* getHoughspace();
  void clear();
  void setOrientation(double orientation);
  double calculateError();
  double getOrientationDrift();
  
  //quality outputs
  double getBasinWidthDiff();
  double getBasinHeightDiff();
  double getMeanSqErr();
  double getSupportRatio();
  
private:
  inline int accumulateValue(double deltaAngle, int deltaDst, SonarPeak& peak);
  bool isLocalMaximum(int angleIdx, int dstIdx);
  void postprocessLines();
  Config config;
  Houghspace houghspace;
  Filter filter;
  double accAsq;
  double accMax;
  double accDsq;
  int dMax;
  std::vector<SonarPeak> allPeaks;
  std::vector<Line> actualLines;
  std::pair<double,double> actualPosition;
  double startAngle;
  int scanDirection;
  bool halfDone;
  double lastAngle;
  int localRangeDstIdx;
  int localRangeAngleIdx;
  double angleTolerance;
  double lastSpatialResolution;
  base::Angle lastOrientation;
  base::Angle firstOrientation;
  double orientationDrift;
  bool orientationDriftDetected;
  
  //the farest delta angle being accumulated (1)
  int accAngleIdx;
  int accDIdx;
  
  //quality outputs
  double basinWidthDiff;
  double basinHeightDiff;
  double meanSqErr;
  double supportRatio;
};

}

#endif // _SONAR_WALL_HOUGH_HOUGH_HPP_