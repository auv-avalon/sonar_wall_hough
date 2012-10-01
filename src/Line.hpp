#ifndef _SONAR_WALL_HOUGH_LINE_HPP_
#define _SONAR_WALL_HOUGH_LINE_HPP_

#include <vector>
#include <map>
#include <algorithm>
#include <math.h>
#include <iostream>
#include "Config.hpp"
#include <base/eigen.h>
#include "Houghspace.hpp"

namespace sonar_wall_hough
{
  struct LinePair;
  
  /**
   * The class Line represents a line by the angle of its normal and its distance to the origin
   * Additionally, each line coontains a value that shows how many votes in hough-transform this line received
   */
  class Line
  {
  public:
    //the angle of the normal
    double alpha;
    
    //the distance of the line to the origin
    int d;
    
    //the number of votes from hough-transform
    int votes;
    
    /**
     * ctor for a line
     * @param alpha the angle of the normal
     * @param d the distance of the line to the origin
     * @param votes the number of votes from hough-transform
     */
    Line(double alpha, int d, int votes);
    
    /**
     * operator to sort lines. A line is smaller than another line if the angle is smaller or if equal the distance
     * is smaller or if equal the votes are less
     * @param other another line
     * @return true if this line is smaller than the other
     */
    bool operator<(const Line& other) const;
    
    /**
     * calculates two cartesian end points for the line
     * @param limitA another line, the intersection makes the first end point
     * @param limitB another line, the intersection makes the second end point
     * @return a pair of end points
     */
    std::pair<base::Vector3d, base::Vector3d> toCartesian(const Line& limitA, const Line& limitB);
    
    /**
     * a line selection algorithm
     * @deprecated
     */
    static std::vector<Line> selectLines(std::vector<Line> lines, std::pair<int,int> basinSize, double spatialResolution, double angleTolerance, double orientation, bool alignLines, bool guessMissing);
    
    /**
     * the line selecting algorithm
     * @param lines all lines found by hough-transform
     * @param basinSize the size (width, height) of the basin
     * @param spatialResolution the resolution of the sonar for calculating from bins to meters
     * @param angleTolerance the lines may be about that off of the basin orientation
     * @param basinOrientation the orientation of the basin with respect to the vehicle (from FOG)
     * @param houghspace the houghspace to create some other lines if needed
     * @param orientation an output parameter that shows the difference of the given basinOrientation and the orientation gotten from the fitting lines
     * @return a vector of the 4 chosen lines
     */
    static std::vector<Line> selectLines2(std::vector<Line> lines, std::pair<int,int> basinSize, double spatialResolution, double angleTolerance, double basinOrientation, Houghspace& houghspace, double*& orientationDrift);
  private:
    /**
     * selects fitting lines by distance, needed for selectLines
     * @deprecated
     */
    static std::vector<Line> selectByDistance(std::vector<Line> lines, std::pair<int,int> basinSize, double spatialResolution);
    
    /**
     * finds a pair of lines that have a given distance to each other
     * @param lines given lines
     * @param distance the desired distance
     * @param houghspace the houghspace to create another line if needed
     * @return the pair of lines
     */
    static std::vector<LinePair> findCorrecpondence(std::vector<Line> lines, int distance, Houghspace& houghspace);
  };
  
  /**
   * This struct represents a pair of lines with their common score
   */
  struct LinePair
  {
    //the Lines
    Line a,b;
    
    //the score
    int score;
    
    /**
     * ctor
     */
    LinePair(Line a, Line b, int score):a(a),b(b),score(score){};
  };
}
#endif // _SONAR_WALL_HOUGH_LINE_HPP_