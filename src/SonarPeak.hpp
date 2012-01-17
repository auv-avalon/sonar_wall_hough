#ifndef _SONAR_WALL_HOUGH_SONAR_PEAK_HPP_
#define _SONAR_WALL_HOUGH_SONAR_PEAK_HPP_

#include <base/angle.h>
#include <base/samples/sonar_beam.h>

namespace sonar_wall_hough
{
    class SonarPeak
    {
    public:
      base::Angle alpha;
      boost::uint16_t distance;
      boost::uint8_t strength;
      
      SonarPeak():alpha(base::Angle::fromRad(0)),distance(0),strength(0){};
      SonarPeak(base::Angle alpha, boost::uint16_t distance, boost::uint8_t strength):alpha(alpha),distance(distance),strength(strength){};
    };

}
#endif // _SONAR_WALL_HOUGH_SONAR_PEAK_HPP_
