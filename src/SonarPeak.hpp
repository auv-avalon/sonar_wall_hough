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
      base::Angle sobelDirection;
      static base::samples::SonarBeam preprevious;
      static base::samples::SonarBeam previous;
      static base::samples::SonarBeam actual;
      static base::samples::SonarBeam next;
      
      static base::samples::SonarBeam preprevious2;
      static base::samples::SonarBeam previous2;
      static base::samples::SonarBeam actual2;
      static base::samples::SonarBeam next2;
      
      SonarPeak();
      SonarPeak(base::Angle alpha, boost::uint16_t distance, boost::uint8_t strength, base::Angle sobelDirection);
      
      static std::vector<SonarPeak> preprocessSonarBeam(base::samples::SonarBeam afternext, int closestBin);
      static std::vector<SonarPeak> preprocessSonarBeam2(base::samples::SonarBeam afternext2, int closestBin);
    };

}
#endif // _SONAR_WALL_HOUGH_SONAR_PEAK_HPP_
