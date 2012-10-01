#ifndef _SONAR_WALL_HOUGH_HOUGHSPACE_HPP_
#define _SONAR_WALL_HOUGH_HOUGHSPACE_HPP_

#include "Config.hpp"
//#include "SonarPeak.hpp"
#include <iostream>
#include <base/samples/sonar_beam.h>
//#include "Line.hpp"

namespace sonar_wall_hough
{
  /**
   * This class represents the houghspace and provides different methods to access it
   * It also holds methods to calculate between real parameters and indices used by the houghspace
   */
  class Houghspace
  {   
  public:
    /**
     * ctor, creates the houghspace using different params from the given config
     * @param config the config
     */
    Houghspace(const Config &config);
    
    /**
     * dtor
     */
    ~Houghspace();
    //boost::uint8_t* at(Line line);
    
    /**
     * returns the hough space cell for given parameters
     * @param angle the desired angle as real angle
     * @param dst the desired distance as real distance
     * @return pointer to the desired hough space cell or NULL if it does not exist
     */
    boost::uint8_t* at(double angle, int dst);
    
    /**
     * returns the hough space cell for given indices
     * @param angleIdx the desired angle as angle index
     * @param dstIdx the desired distance as distance index
     * @return pointer to the desired hough space cell or NULL if it does not exist
     */
    boost::uint8_t* at(int angleIdx, int dstIdx);
    
    /**
     * returns the hough space cell for given indices without checking
     * if that cell exists (faster)
     * @param angleIdx the desired angle as angle index
     * @param dstIdx the desired distance as distance index
     * @return pointer to the desired hough space cell
     */
    boost::uint8_t* uncheckedAt(int angleIdx, int dstIdx);
    
    /**
     * clears the whole hough space
     */
    void clear();
    
    /**
     * @return the width of the hough space
     */
    int getWidth();
    
    /**
     * @return the height of the hough space
     */
    int getHeight();
    
    /**
     * @return the number of angles of the hough space
     */
    int getNumberOfAngles();
    
    /**
     * @return the number of positive distances of the hough space
     */
    int getnumberOfPosDistances();
    
    /**
     * calculates a given real angle to an index
     * returns an index between -numberOfAngles and +numberOfAngles (numberOfAngles corresponds to 180 degree)
     * @param angle the real angle in radians
     * @return the corresponding index for this angle
     */
    inline int angle2Idx(double angle)
    {
    //normalize angle to (-M_PI, M_PI]
    /*angle = fmod(angle, 2*M_PI);
    if(angle > M_PI)
      angle -= 2*M_PI;
    */
    while(angle >= M_PI)
      angle -= 2*M_PI;
    while(angle < -M_PI)
      angle += 2*M_PI;
    return (int)(angle * numberOfAngles / M_PI + 0.5);
    }
    
    /**
     * calculates a given real distance to an index
     * returns an index between -numberOfPosDistances and +numberOfPosDistances if abs(dst) < config.maxDistance
     * @param dst the real distance in bins
     * @return the corresponding index for this distance
     */
    inline int dst2Idx(int dst)
    {
      return (int)(dst / config.distancesPerBin + 0.5);
    }
    
    /**
     * calculates a real angle from a given index
     * returns an angle between -M_PI and M_PI if angleIdx in [0, numberOfAngles]
     * @param angleIdx the angle index
     * @return the corresponding angle in radians for this index
     */
    inline double idx2Angle(int angleIdx)
    {
      return M_PI * angleIdx / numberOfAngles;
    }
    
    /**
     * calculates a real distance from a given index
     * returns an signed distance corresponding to the index if abs(dstIdx) < numberOfPosDistances
     * @param dstIdx the distance index
     * @return the corresponding distance in bins for this index
     */
    inline int idx2Dst(int dstIdx)
    {
      return dstIdx * config.distancesPerBin;
    }
    
  private:
    //the config, set in ctor
    Config config;
    
    //the storage for the hough space entries
    boost::uint8_t* space;
    
    //the number of angles
    int numberOfAngles;
    
    //the number of positive distances
    int numberOfPosDistances;
  };  

}
#endif // _SONAR_WALL_HOUGH_HOUGHSPACE_HPP_