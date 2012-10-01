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

/**
 * This class holds methods to perform the whole hough-localization
 */
class Hough
{

public:
  /**
   * ctor
   * @param config the config
   */
  Hough(const Config& config = Config());
  
  /**
   * dtor
   */
  ~Hough();
  
  /**
   * this method takes a peak and accumulates the hough space appropriately
   * @param peaks a vector of peaks that all have the same angle (for faster accumulation)
   */
  void accumulate(std::vector<SonarPeak> peaks);
  
  /**
   * this method receives a beam and invokes filtering and accumulation
   * it also detects when the hough space should be analyzed (i.e. after a full 360 degree scan)
   * @param beam the incoming sonar beam
   */
  void registerBeam(base::samples::SonarBeam beam);
  
  /**
   * this method searches for local maxima in the hough space with values above the threshold
   * and invokes the post-processing for the so found lines
   */
  void analyzeHoughspace();
  
  /**
   * @return all peaks that have been detected so far
   */
  std::vector<SonarPeak>* getAllPeaks();
  
  /**
   * @return the last found lines of hte basin
   */
  std::vector<Line>* getActualLines();
  
  /**
   * @return the last calculated position of the vehicle
   */
  std::pair<double,double> getActualPosition();
  
  /**
   * @return the actual orientation of the vehicle
   */
  base::Angle getOrientation();
  
  /**
   * @return a pointer to the hough space used here
   */
  Houghspace* getHoughspace();
  
  /**
   * invokes hough space clearing, actual lines clearing and sets the starting orientation to the actual orientation
   */
  void clear();
  
  /**
   * sets the actual orientation
   * @param orientation the actual orientation
   */
  void setOrientation(double orientation);
  
  /**
   * calculates the different values to measure the quality of the measurement
   */
  void calculateError();
  
  /**
   * @return the estimated drift of the input orientation
   */
  double getOrientationDrift();
  
  /**
   * @return the difference of given and measured basin width
   */
  double getBasinWidthDiff();
  
  /**
   * @return the difference of given and measured basin height
   */
  double getBasinHeightDiff();
  
  /**
   * @return the mean square error of the peaks to the detected lines of the basin
   */
  double getMeanSqErr();
  
  /**
   * @return the ratio of peaks that have a distance of less than 20 to the closest wall
   */
  double getSupportRatio();
  
private:
  /**
   * calculates the value the hough space has to be accumulated by
   * @param deltaAngle the difference of the peaks' angle and the angle of the hough space cell
   * @param deltaDst the difference of the peaks' distance and the distance of the hough space cell
   * @param peak the sonar peak
   * @return the accumulation value
   */
  inline int accumulateValue(double deltaAngle, int deltaDst, SonarPeak& peak);
  
  /**
   * determines if the given indices point to a local maximum in the hough space
   * @param angleIdx the angle index
   * @param dstIdx the distance index
   * @return true if we have a local maximum, false if not
   */
  bool isLocalMaximum(int angleIdx, int dstIdx);
  
  /**
   * performes the post processing of the detected lines
   * invokes the line selection, position calculation and quality value determination
   */
  void postprocessLines();
  
  //the config, set by ctor
  Config config;
  
  //the used hough space
  Houghspace houghspace;
  
  //the used filter for the incoming beams
  Filter filter;
  
  //constant for exponential function for accumulation (angle)
  double accAsq;
  
  //maximum accumulation value
  double accMax;
  
  //constant for exponential function for accumulation (distance)
  double accDsq;
  
  //maximum distance for accumulation
  int dMax;
  
  //stores all peaks since last clear
  std::vector<SonarPeak> allPeaks;
  
  //stores the actual lines
  std::vector<Line> actualLines;
  
  //holds the actual position
  std::pair<double,double> actualPosition;
  
  //holds the start angle of the beam to detect full 360 degree scan
  double startAngle;
  
  //holds the scan direction to detect ping-pong mode
  int scanDirection;
  
  //is true if we already did more than 180 degree scan (makes full scan detection easier)
  bool halfDone;
  
  //stores the last sonar beams' angle to detect scan direction
  double lastAngle;
  
  //the distance range within the hough space cells are being accumulated for a peak
  int localRangeDstIdx;
  
  //the angle range within the hough space cells are being accumulated for a peak
  int localRangeAngleIdx;
  
  //stores the spatial resolution of last registered beam
  double lastSpatialResolution;
  
  //stores last given orientation of vehicle
  base::Angle lastOrientation;
  
  //stores first orientation of vehicle after clear
  base::Angle firstOrientation;
  
  //stores the drift of the given orientations vs the measured orientations
  double orientationDrift;
  
  //true, if the drit is above a threshold
  bool orientationDriftDetected;
  
  //the farest delta angle being accumulated (1)
  int accAngleIdx;
  
  //the farest delta distance being accumulated
  int accDIdx;
  
  //difference of given and measured basin width
  double basinWidthDiff;
  
  //difference of given and measured basin height
  double basinHeightDiff;
  
  //the mean square error of the peaks to the detected lines of the basin
  double meanSqErr;
  
  //the ratio of peaks that have a distance of less than 20 to the closest wall
  double supportRatio;
};

}

#endif // _SONAR_WALL_HOUGH_HOUGH_HPP_