#ifndef _SONAR_WALL_HOUGH_SONAR_PEAK_HPP_
#define _SONAR_WALL_HOUGH_SONAR_PEAK_HPP_

#include <base/angle.h>
#include <base/samples/sonar_beam.h>
#include <limits>
#include "Line.hpp"

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
      
      double dstFromLine(const Line& line) const
      {
	return fabs(line.d - distance * sin(M_PI/2 + line.alpha - alpha.getRad()));
      };
      
      std::pair<int,double> smallestDstFromLine(const std::vector<Line>& lines) const
      {
	double min = std::numeric_limits<double>::infinity();
	int whichLine = -1;
	for(int i = 0; i < lines.size(); i++)
	{
	  double dst = dstFromLine(lines[i]);
	  if(dst < min)
	  {
	    min = dst;
	    whichLine = i;
	  }
	}
	return std::pair<int,double>(whichLine, min);
      };
      
      base::Vector3d toCartesian()
      {
	double x = distance * cos(alpha.getRad());
	double y = distance * sin(alpha.getRad());
	return base::Vector3d(x, y, 0.0);
      };
    };

}
#endif // _SONAR_WALL_HOUGH_SONAR_PEAK_HPP_
