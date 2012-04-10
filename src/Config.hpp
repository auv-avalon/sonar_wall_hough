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
    double minDistance;
    int distancesPerBin;
    double minLineVotesRatio;
    double angleDelta;
    double basinHeight;
    double basinWidth;

    Config()
	: sensorAngularResolution(0.0)
	, filterThreshold(0)
	, anglesPerBin(0)
	, maxDistance(0)
	, minDistance(0.0)
	, distancesPerBin(0)
	, minLineVotesRatio(0.0)
	, angleDelta(0)
	, basinHeight(0.0)
	, basinWidth(0.0)
    {
    }   

  };

}
#endif // _SONAR_WALL_HOUGH_CONFIG_HPP_