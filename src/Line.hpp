#ifndef _SONAR_WALL_HOUGH_LINE_HPP_
#define _SONAR_WALL_HOUGH_LINE_HPP_

namespace sonar_wall_hough
{
  struct Line
  {   
      double alpha;
      int d;
      int votes;
      
      Line(double alpha, int d, int votes)
	:alpha(alpha)
	,d(d)
	,votes(votes)
      {
      }
      
  };  

}
#endif // _SONAR_WALL_HOUGH_LINE_HPP_
