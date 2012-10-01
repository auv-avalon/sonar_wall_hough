#ifndef _SONAR_WALL_HOUGH_FILTER_HELPER_HPP_
#define _SONAR_WALL_HOUGH_FILTER_HELPER_HPP_

#include <deque>
#include <base/angle.h>
#include <base/samples/sonar_beam.h>
#include "SonarPeak.hpp"

namespace sonar_wall_hough
{
  enum FilterType{minimum, sobelGaussDst, sobelGaussPhi};
  
  /**
   * This class is a helper for minimum filtering
   * see http://arxiv.org/abs/cs.DS/0610046 for details
   */
  class MinWindow
  {
  public:
    /**
     * The ctor
     * @param kernelSize is the size of the range whithin that the minimum value is choosen
     */
    MinWindow(int kernelSize);
    
    /**
     * dtor
     */
    ~MinWindow();
    
    /**
     * proceeds the minimum window
     * @param value is the next value
     * @return the next minimum value
     */
    uint8_t pushValue(uint8_t value);
    
    /**
     * only gives the last minimums without feeding new values in
     * @return the next minimum value
     */
    uint8_t pushValue();
    
    /**
     * clears the window
     */
    void clear();
  private:
    //size of the window, set by ctor
    int kernelSize;
    
    //deque for minimum values
    std::deque<uint8_t> minDeque;
    
    //deque for minimum value indices
    std::deque<int> idxDeque;
    
    //last index
    int lastIdx;
  };
  
  /**
   * This class is for performing filtering on a sonar beam whithin its direction
   */
  class BeamFilterDst
  {
  public:
    /**
     * ctor
     * @param type specifies the type of the filter (minimum or sobelGaussDst)
     * @param kernelSize sets the size of the kernel
     */
    BeamFilterDst(FilterType type, int kernelSize);
    
    /**
     * dtor
     */
    ~BeamFilterDst();
    
    /**
     * performs the filtering
     * @param sonarBeam the beam to be filtered
     * @return the filtered beam
     */
    base::samples::SonarBeam filter(base::samples::SonarBeam sonarBeam);
    
  private:
    //the type, set by ctor
    FilterType type;
    
    //the size of the kernel, set by ctor
    int kernelSize;
    
    //A minWindow, only needed for minimum filtering
    MinWindow* minWindow;
  };
  
  /**
   * This class filters a beam in perpendicular direction of the beam
   * For that, multiple successive beams are collected and the entries of the beam in the middle are filtered
   * by considering the apporopriate entries of the neighboured beams
   */
  class BeamFilterPhi
  {
  public:
    /**
     * ctor
     * @param type specifies the type of the filter (minimum or sobelGaussPhi)
     * @param kernelSize sets the size of the kernel
     */
    BeamFilterPhi(FilterType type, int kernelSize);
    
    /**
     * dtor
     */
    ~BeamFilterPhi();
    
    /**
     * performs the filtering
     * @param sonarBeam the next beam to be stored
     * @return the filtered middle beam
     */
    base::samples::SonarBeam filter(base::samples::SonarBeam sonarBeam);
    
  private:
    //the type, set by ctor
    FilterType type;
    
    //the size of the kernel, set by ctor
    int kernelSize;
    
    //the storage for the beams
    std::deque<base::samples::SonarBeam> lastbeams;
    
    //A minWindow, only needed for minimum filtering
    MinWindow* minWindow;
  };
  
  /**
   * This class provides the whole filtering of incoming beams performing minimum, sobel and gauss filtering
   * sonarBeams are provided to the class and it performs all filtering and a local maximum search.
   * Local maxima that exceed a given distance and strength will be called peaks and are returned afterwards.
   * @note the returned peaks will not result from the actual given beam but from the kernelSize/2-th last given beam
   */
  class Filter
  {
  public:
    /**
     * ctor
     * @param kernelSize sets the size of the kernel
     * @param threshold the minimum value of an entry to be called a peak after filtering
     * @param withMinimum if a minimum filtering should be performed or not
     */
    Filter(int kernelSize, uint8_t threshold, bool withMinimum);
    
    /**
     * dtor
     */
    ~Filter();
    
    /**
     * performs the filtering
     * @param sonarBeam the next incoming sonarBeam
     * @param minDistance the minimum distance of an entry to be called a peak after filtering
     * @return a vector of the peaks after filtering
     */
    std::vector<SonarPeak> filter(base::samples::SonarBeam sonarBeam, int minDistance);
  
  private:
    //the threshold, set by ctor
    uint8_t threshold;
    
    //minimum filtering, set by ctor
    bool withMinimum;
    
    //the minimum filter in dst-direction
    BeamFilterDst filterDstMin;
    
    //the sobel-gauss-filter in dst-direction, edges in dst-direction
    BeamFilterDst filterDstSGDst;
    
    //the sobel-gauss-filter in dst-direction, edges in phi-direction
    BeamFilterDst filterDstSGPhi;
    
    //the minimum filter in phi-direction
    BeamFilterPhi filterPhiMin;
    
    //the sobel-gauss-filter in phi-direction, edges in dst-direction
    BeamFilterPhi filterPhiSGDst;
    
    //the sobel-gauss-filter in phi-direction, edges in phi-direction
    BeamFilterPhi filterPhiSGPhi;
  };
}
#endif // _SONAR_WALL_HOUGH_FILTER_HELPER_HPP_
