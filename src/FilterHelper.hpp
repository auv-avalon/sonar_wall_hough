#ifndef _SONAR_WALL_HOUGH_FILTER_HELPER_HPP_
#define _SONAR_WALL_HOUGH_FILTER_HELPER_HPP_

#include <deque>
#include <base/angle.h>
#include <base/samples/sonar_beam.h>
#include "SonarPeak.hpp"

namespace sonar_wall_hough
{
  enum FilterType{minimum, sobelGaussDst, sobelGaussPhi};
  
  class MinWindow
  {
  public:
    MinWindow(int kernelSize);
    ~MinWindow();
    uint8_t pushValue(uint8_t value);
    void clear();
  private:
    int kernelSize;
    std::deque<uint8_t> minDeque;
  };
  
  class BeamFilterDst
  {
  public:
    BeamFilterDst(FilterType type, int kernelSize);
    ~BeamFilterDst();
    base::samples::SonarBeam filter(base::samples::SonarBeam sonarBeam);
    
  private:
    FilterType type;
    int kernelSize;
    MinWindow* minWindow;
  };
  
  class BeamFilterPhi
  {
  public:
    BeamFilterPhi(FilterType type, int kernelSize);
    ~BeamFilterPhi();
    base::samples::SonarBeam filter(base::samples::SonarBeam sonarBeam);
    
  private:
    FilterType type;
    int kernelSize;
    std::deque<base::samples::SonarBeam> lastbeams;
    MinWindow* minWindow;
  };
  
  class Filter
  {
  public:
    Filter(int kernelSize);
    ~Filter();
    std::vector<SonarPeak> filter(base::samples::SonarBeam sonarBeam, int minDistance);
  
  private:
    BeamFilterDst filterDstMin;
    BeamFilterDst filterDstSGDst;
    BeamFilterDst filterDstSGPhi;
    BeamFilterPhi filterPhiMin;
    BeamFilterPhi filterPhiSGDst;
    BeamFilterPhi filterPhiSGPhi;
  };
}
#endif // _SONAR_WALL_HOUGH_FILTER_HELPER_HPP_
