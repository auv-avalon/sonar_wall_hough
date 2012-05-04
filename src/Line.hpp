#ifndef _SONAR_WALL_HOUGH_LINE_HPP_
#define _SONAR_WALL_HOUGH_LINE_HPP_

#include <vector>
#include <map>
#include <algorithm>
#include <math.h>
#include <iostream>
#include "Config.hpp"
#include <base/eigen.h>

namespace sonar_wall_hough
{
  class Line
  {
  public:
    double alpha;
    int d;
    int votes;
    
    Line(double alpha, int d, int votes);
    bool operator<(const Line& other) const;
    std::pair<base::Vector3d, base::Vector3d> toCartesian(const Line& limitA, const Line& limitB);
    
    static std::vector<Line> selectLines(std::vector<Line> lines, std::pair<int,int> basinSize, double spatialResolution, double angleTolerance, double orientation, bool alignLines, bool guessMissing);
  private:
    static std::vector<Line> selectByDistance(std::vector<Line> lines, std::pair<int,int> basinSize, double spatialResolution);
  };
  

}
#endif // _SONAR_WALL_HOUGH_LINE_HPP_