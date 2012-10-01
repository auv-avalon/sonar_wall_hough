#ifndef _SONAR_WALL_HOUGH_SONAR_PEAK_HPP_
#define _SONAR_WALL_HOUGH_SONAR_PEAK_HPP_

#include <base/angle.h>
#include <base/samples/sonar_beam.h>
#include <limits>
#include "Line.hpp"

namespace sonar_wall_hough
{
    /**
     * This class holds sonar peaks and some methods for them. A SonarPeak is given in polar-coordinates and has a distinct strength
     */
    class SonarPeak
    {
    public:
      //the angle of the peak
      base::Angle alpha;
      
      //the distance of the peak
      boost::uint16_t distance;
      
      //the strength of the peak
      boost::uint8_t strength;
      
      /**
       * ctor for empty peak, everything is set to zero
       */
      SonarPeak():alpha(base::Angle::fromRad(0)),distance(0),strength(0){};
      
      /**
       * ctor for peak with params
       * @param alpha the angle of the peak
       * @param distance the distance of the peak
       * @param strength the strength of the peak
       */
      SonarPeak(base::Angle alpha, boost::uint16_t distance, boost::uint8_t strength):alpha(alpha),distance(distance),strength(strength){};
      
      /**
       * calculates the distance of the peak from a given line
       * @param line the given line
       * @return the distance to that line
       */
      double dstFromLine(const Line& line) const
      {
	return fabs(line.d - distance * sin(M_PI/2 + line.alpha - alpha.getRad()));
      };
      
      /**
       * calculates the smallest distance from multiple lines
       * @param lines a vector of lines
       * @return a pair. first: the index of the closest line, second: the distance to that line
       */
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
      
      /**
       * calculates the cartesian coordinates for the peak
       * @return the cartesian coordinates
       */
      base::Vector3d toCartesian()
      {
	double x = distance * cos(alpha.getRad());
	double y = distance * sin(alpha.getRad());
	return base::Vector3d(x, y, 0.0);
      };
    };

}
#endif // _SONAR_WALL_HOUGH_SONAR_PEAK_HPP_
