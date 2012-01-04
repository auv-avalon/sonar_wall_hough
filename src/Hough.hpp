#ifndef _SONAR_WALL_HOUGH_HOUGH_HPP_
#define _SONAR_WALL_HOUGH_HOUGH_HPP_

#include <base/samples/sonar_beam.h>
#include <base/samples/frame.h>
#include <math.h>
#include <limits.h>
#include "Config.hpp"
#include "Line.hpp"
#include "Houghspace.hpp"
#include "SonarPeak.hpp"

namespace sonar_wall_hough 
{

class Hough
{

public:
  Hough(const Config& config = Config());
  ~Hough();
  void accumulate(SonarPeak peak);
  void registerBeam(base::samples::SonarBeam beam);
  void analyzeHoughspace();
  std::vector<SonarPeak>* getAllPeaks();
  std::vector<Line>* getActualLines();
  Houghspace* getHoughspace();
  void clear();
  
private:
  inline int accumulateValue(double deltaAngle);
  bool isLocalMaximum(int angleIdx, int dstIdx);
  Config config;
  Houghspace houghspace;
  double accAsq;
  double accMax;
  std::vector<SonarPeak> allPeaks;
  std::vector<Line> actualLines;
  double lastAnalysisAngle;
  int localRangeDstIdx;
  int localRangeAngleIdx;
  
  //the farest delta angle being accumulated (1)
  int accAngleIdx;
};

}

#endif // _SONAR_WALL_HOUGH_HOUGH_HPP_