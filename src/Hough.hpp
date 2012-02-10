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
  Houghspace* getHoughspace();
  void clear();
  
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
  double startAngle;
  int scanDirection;
  double lastAngle;
  int localRangeDstIdx;
  int localRangeAngleIdx;
  double angleTolerance;
  double lastSpatialResolution;
  
  //the farest delta angle being accumulated (1)
  int accAngleIdx;
  int accDIdx;
};

}

#endif // _SONAR_WALL_HOUGH_HOUGH_HPP_