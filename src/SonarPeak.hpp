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
      
      SonarPeak()
	:alpha()
	,distance()
	,strength()
      {
      }
      
      SonarPeak(base::Angle alpha, boost::uint16_t distance, boost::uint8_t strength)
	:alpha(alpha)
	,distance(distance)
	,strength(strength)
      {
      }
      
      static std::vector<SonarPeak> preprocessSonarBeam(base::samples::SonarBeam sonarBeam)
      {
	int minDistance = 25;
	int sonarThreshold = 20;
	//int minDistBetweenMaxima = 50; //TODO: dynamisch
	
	std::vector<SonarPeak> peaks;
	//beam must have min 3 entries
	if(sonarBeam.beam.size() < 3)
	  return peaks;
	
	for(int i = minDistance; i < (int)sonarBeam.beam.size()-1; i++)
	{
	  if(sonarBeam.beam.at(i) >= sonarThreshold)
	  {
	    if(sonarBeam.beam.at(i) > sonarBeam.beam.at(i-1) && sonarBeam.beam.at(i) > sonarBeam.beam.at(i+1))
	    {
	      //found local maximum above threshold
	      peaks.push_back(SonarPeak(sonarBeam.bearing, i, sonarBeam.beam.at(i)));
	    }
	  }
	}
	//is last entry possible peak?
	if(sonarBeam.beam.back() > sonarThreshold && sonarBeam.beam.back() > sonarBeam.beam.at(sonarBeam.beam.size()-2))
	{
	  peaks.push_back(SonarPeak(sonarBeam.bearing, sonarBeam.beam.size()-1, sonarBeam.beam.back()));
	}
	
	/*
	//if peaks are too close, erase smaller one
	std::vector<SonarPeak>::iterator it = peaks.begin();
	while(it < peaks.end()-1)
	{
	  if((it+1)->distance - it->distance < minDistBetweenMaxima)
	  {
	    if(it->strength > (it+1)->strength)
	    {
	      peaks.erase(it+1);
	    }
	    else
	    {
	      peaks.erase(it);
	      it++;
	    }
	  }
	  else
	  {
	    it++;
	  }
	}
	      return peaks;*/
      }
    };

}
#endif // _SONAR_WALL_HOUGH_SONAR_PEAK_HPP_
