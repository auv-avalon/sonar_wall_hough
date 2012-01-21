#ifndef _SONAR_WALL_HOUGH_CONFIG_HPP_
#define _SONAR_WALL_HOUGH_CONFIG_HPP_

#include <iostream>
#include <stdint.h>

namespace sonar_wall_hough
{
  struct Config
  {
    double sensorAngularResolution;
    uint8_t filterThreshold;
    int anglesPerBin;
    int maxDistance;
    int distancesPerBin;
    int minLineVotes;

    Config()
	: sensorAngularResolution(0)
	, filterThreshold(0)
	, anglesPerBin(0)
	, maxDistance(0)
	, distancesPerBin(0)
	, minLineVotes(0)
    {
    }   

  };

}
#endif // _SONAR_WALL_HOUGH_CONFIG_HPP_