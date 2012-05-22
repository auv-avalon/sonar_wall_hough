#ifndef _SONAR_WALL_HOUGH_HOUGHSPACE_HPP_
#define _SONAR_WALL_HOUGH_HOUGHSPACE_HPP_

#include "Config.hpp"
//#include "SonarPeak.hpp"
#include <iostream>
#include <base/samples/sonar_beam.h>
//#include "Line.hpp"

namespace sonar_wall_hough
{
  class Houghspace
  {   
  public:
    Houghspace(const Config &config);
    ~Houghspace();
    //boost::uint8_t* at(Line line);
    boost::uint8_t* at(double angle, int dst);
    boost::uint8_t* at(int angleIdx, int dstIdx);
    boost::uint8_t* uncheckedAt(int angleIdx, int dstIdx);
    void clear();
    int getWidth();
    int getHeight();
    int getNumberOfAngles();
    int getnumberOfPosDistances();
    
    //returns an index between -numberOfAngles and +numberOfAngles (numberOfAngles corresponds to 180 degree)
    inline int angle2Idx(double angle)
    {
    //normalize angle to (-M_PI, M_PI]
    angle = fmod(angle, 2*M_PI);
    if(angle > M_PI)
      angle -= 2*M_PI;

    return (int)(angle * numberOfAngles / M_PI + 0.5);
    }
    
    //returns an index between -numberOfPosDistances and +numberOfPosDistances if abs(dst) < config.maxDistance
    inline int dst2Idx(int dst)
    {
      return (int)(dst / config.distancesPerBin + 0.5);
    }
    //returns an angle between -M_PI and M_PI if angleIdx in [0, numberOfAngles]
    inline double idx2Angle(int angleIdx)
    {
      return M_PI * angleIdx / numberOfAngles;
    }
    //returns an signed distance corresponding to the index if abs(dstIdx) < numberOfPosDistances
    inline int idx2Dst(int dstIdx)
    {
      return dstIdx * config.distancesPerBin;
    }
    
  private:
    Config config;
    boost::uint8_t* space;
    int numberOfAngles;
    int numberOfPosDistances;
  };  

}
#endif // _SONAR_WALL_HOUGH_HOUGHSPACE_HPP_