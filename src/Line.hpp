#ifndef _SONAR_WALL_HOUGH_LINE_HPP_
#define _SONAR_WALL_HOUGH_LINE_HPP_

#include <vector>
#include <algorithm>
#include <math.h>
#include <iostream>
#include "Config.hpp"

namespace sonar_wall_hough
{
  class Line
  {
  public:
    double alpha;
    int d;
    int votes;
    
    Line(double alpha, int d, int votes);
    
    static std::vector<Line> selectLines(std::vector<Line> lines, std::vector<int> validDistances);
  };
  

}
#endif // _SONAR_WALL_HOUGH_LINE_HPP_
