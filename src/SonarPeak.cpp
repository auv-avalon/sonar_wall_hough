#include "Hough.hpp"

namespace sonar_wall_hough
{
  base::samples::SonarBeam SonarPeak::preprevious = base::samples::SonarBeam();
  base::samples::SonarBeam SonarPeak::previous = base::samples::SonarBeam();
  base::samples::SonarBeam SonarPeak::actual = base::samples::SonarBeam();
  base::samples::SonarBeam SonarPeak::next = base::samples::SonarBeam();
  
  base::samples::SonarBeam SonarPeak::preprevious2 = base::samples::SonarBeam();
  base::samples::SonarBeam SonarPeak::previous2 = base::samples::SonarBeam();
  base::samples::SonarBeam SonarPeak::actual2 = base::samples::SonarBeam();
  base::samples::SonarBeam SonarPeak::next2 = base::samples::SonarBeam();
  
  SonarPeak::SonarPeak()
    :alpha()
    ,distance()
    ,strength()
    ,sobelDirection()
  {
  }
  
   SonarPeak::SonarPeak(base::Angle alpha, boost::uint16_t distance, boost::uint8_t strength, base::Angle sobelDirection)
    :alpha(alpha)
    ,distance(distance)
    ,strength(strength)
    ,sobelDirection(sobelDirection)
  {
  }
  
  std::vector<SonarPeak> SonarPeak::preprocessSonarBeam(base::samples::SonarBeam afternext, int closestBin)
  {
    base::samples::SonarBeam filtered;
    
    if(preprevious.beam.empty())
    {
      preprevious = previous;
      previous = actual;
      actual = next;
      next = afternext;
      
      return std::vector<SonarPeak>();
    }
    
    filtered.bearing = actual.bearing;
    filtered.beam.resize(actual.beam.size());
    
    //make filtering (gauss * sobel)
    int resDst, resPhi;
    std::vector<uint8_t>::iterator flt = filtered.beam.begin()+closestBin;
    
    std::vector<uint8_t>::iterator ppr = preprevious.beam.begin()+closestBin;
    std::vector<uint8_t>::iterator prv = previous.beam.begin()+closestBin;
    std::vector<uint8_t>::iterator nxt = next.beam.begin()+closestBin;
    std::vector<uint8_t>::iterator anx = afternext.beam.begin()+closestBin;
    for(std::vector<uint8_t>::iterator act = actual.beam.begin()+closestBin;act < actual.beam.end()-1; flt++, act++, ppr++, prv++, nxt++, anx++)
    {
      
      /*    ppr prv act nxt anx
       * 2 [ 1   2   0  -2  -1 ]
       * 1 [ 4   8   0  -8  -4 ]
       * 0 [ 4  12   0 -12  -4 ] * beam
       *-1 [ 4   8   0  -8  -4 ]
       *-2 [ 1   2   0  -2  -1 ]
       */
      /*resPhi = (-1 * ppr[2])  + (-1 *prv[2]) + ( 1 *nxt[2]) + ( 1 * anx[2])  +
	       (-3 * ppr[2])  + (-3 *prv[2]) + ( 3 *nxt[2]) + ( 3 * anx[2])  +
	       (-4 * ppr[2])  + (-4 *prv[2]) + ( 4 *nxt[2]) + ( 4 * anx[2])  +
	       (-3 * ppr[2])  + (-3 *prv[2]) + ( 3 *nxt[2]) + ( 3 * anx[2])  +
	       (-1 * ppr[2])  + (-1 *prv[2]) + ( 1 *nxt[2]) + ( 1 * anx[2]);
      resPhi /= 16;*/
      //resPhi = ((1*prv[1]) + (-1*nxt[1]) + (2*prv[0]) + (-2*nxt[0]) + (1*prv[-1]) + (-1*nxt[-1])/8);
      //resPhi = sqrt(2) * act[0];
      
      resPhi = std::min<uint8_t>(ppr[2],std::min<uint8_t>(prv[2], std::min<uint8_t>(act[2], std::min<uint8_t>(nxt[2], anx[2]))));
      resPhi = std::min<uint8_t>(resPhi, std::min<uint8_t>(ppr[1],std::min<uint8_t>(prv[1], std::min<uint8_t>(act[1], std::min<uint8_t>(nxt[1], anx[1])))));
      resPhi = std::min<uint8_t>(resPhi, std::min<uint8_t>(ppr[0],std::min<uint8_t>(prv[0], std::min<uint8_t>(act[0], std::min<uint8_t>(nxt[0], anx[0])))));
      resPhi = std::min<uint8_t>(resPhi, std::min<uint8_t>(ppr[-1],std::min<uint8_t>(prv[-1], std::min<uint8_t>(act[-1], std::min<uint8_t>(nxt[-1], anx[-1])))));
      resPhi = std::min<uint8_t>(resPhi, std::min<uint8_t>(ppr[-2],std::min<uint8_t>(prv[-2], std::min<uint8_t>(act[-2], std::min<uint8_t>(nxt[-2], anx[-2])))));
      
      /*uint8_t values[25];
      values[0] = ppr[2]; values[1] = prv[2]; values[2] = act[2]; values[3] = nxt[2]; values[4] = anx[2];
      values[5] = ppr[1]; values[6] = prv[1]; values[7] = act[1]; values[8] = nxt[1]; values[9] = anx[1];
      values[10] = ppr[0]; values[11] = prv[0]; values[12] = act[0]; values[13] = nxt[0]; values[14] = anx[0];
      values[15] = ppr[-1]; values[16] = prv[-1]; values[17] = act[-1]; values[18] = nxt[-1]; values[19] = anx[-1];
      values[20] = ppr[-2]; values[21] = prv[-2]; values[22] = act[-2]; values[23] = nxt[-2]; values[24] = anx[-2];
      std::sort(values+0, values+25);
      resPhi = values[12];*/
		  
      /*    ppr prv act nxt anx
       * 2 [ 1   4   4   4   1 ]
       * 1 [ 2   8  12   8   2 ]
       * 0 [ 0   0   0   0   0 ] * beam
       *-1 [-2  -8 -12  -8  -2 ]
       *-2 [-1  -4  -4  -4  -1 ]
       */
      /*resDst = (-1 * ppr[2])  + (-3 *prv[2]) + (-4 * act[2]) + (-3 *nxt[2]) + (-1 * anx[2])  +
	       (-1 * ppr[2])  + (-3 *prv[2]) + (-4 * act[2]) + (-3 *nxt[2]) + (-1 * anx[2])  +
	       ( 1 * ppr[2])  + ( 3 *prv[2]) + ( 4 * act[2]) + ( 3 *nxt[2]) + ( 1 * anx[2])  +
	       ( 1 * ppr[2])  + ( 3 *prv[2]) + ( 4 * act[2]) + ( 3 *nxt[2]) + ( 1 * anx[2]);
      resDst /= 16;*/
      //resDst = ((1*prv[1]) + (2*act[1]) + (1*nxt[1]) + (-1*prv[-1]) + (-2*act[-1]) + (-1*nxt[-1])/8);
      //resDst = sqrt(2) * act[0];
      resDst = 0;
      uint8_t strength = sqrt(resDst*resDst+resPhi*resPhi);
      //std::cout << (int)resPhi << std::endl;
      if(strength > 0)
      {
	double direction = atan2(resDst,resPhi);
	//std::cout << (int)resPhi << ", ";
	flt[0] = strength;
      }
    }
    
    //std::cout << std::endl;
    /*
    //non-maximum suppression (only in Dst-Direction)
    if(peaks.size() > 2)
    {
      for(std::vector<SonarPeak>::iterator it = peaks.begin()+1; it < peaks.end()-1; it++)
      {
	if(abs((*it).distance-(*(it-1)).distance)<=1 && ((*it).strength < (*(it-1)).strength || (*it).strength < (*(it+1)).strength))
	  peaks.erase(it);
      }
    }
    */
    preprevious = previous;
    previous = actual;
    actual = next;
    next = afternext;
    return preprocessSonarBeam2(filtered, closestBin);
  }
  
std::vector<SonarPeak> SonarPeak::preprocessSonarBeam2(base::samples::SonarBeam afternext2, int closestBin)
{
  std::vector<SonarPeak> peaks;
  
  base::samples::SonarBeam filtered;
    
    if(preprevious2.beam.empty())
    {
      preprevious2 = previous2;
      previous2 = actual2;
      actual2 = next2;
      next2 = afternext2;
      
      return peaks;
    }
  
  //do a sobel + gauss filter
  int resDst, resPhi;
  std::vector<uint8_t>::iterator ppr = preprevious2.beam.begin()+closestBin;
  std::vector<uint8_t>::iterator prv = previous2.beam.begin()+closestBin;
  std::vector<uint8_t>::iterator nxt = next2.beam.begin()+closestBin;
  std::vector<uint8_t>::iterator anx = afternext2.beam.begin()+closestBin;
  for(std::vector<uint8_t>::iterator act = actual2.beam.begin()+closestBin;act < actual2.beam.end()-2; act++, ppr++, prv++, nxt++, anx++)
  {
    resPhi = (1 * ppr[2])  + (2 *prv[2]) + ( -2 *nxt[2]) + ( -1 * anx[2])  +
	       (4 * ppr[1])  + (8 *prv[1]) + ( -8 *nxt[1]) + (-4 * anx[1])  +
	       (6 * ppr[0])  + (12 *prv[0]) + (-12 *nxt[0]) + (-6 * anx[0])  +
	       (4 * ppr[-1])  + (2 *prv[-1]) + ( -2 *nxt[-1]) + ( -4 * anx[-1])  +
	       (1 * ppr[-2])  + (2 *prv[-2]) + (-2 *nxt[-2]) + ( -1 * anx[-2]);
    resPhi /= 16;
    
    resDst = (1 * ppr[2])  + (4 *prv[2]) + (6 * act[2]) + (4 *nxt[2]) + (1 * anx[2])  +
	       (2 * ppr[1])  + (8 *prv[1]) + (12 * act[1]) + (8 *nxt[1]) + (2 * anx[1])  +
	       ( -2 * ppr[-1])  + (-8 *prv[-1]) + (-12* act[-1]) + (-8 *nxt[-1]) + (-2 * anx[-1])  +
	       (-1 * ppr[-2])  + (-4 *prv[-2]) + (-6 * act[-2]) + (-4 *nxt[-2]) + (-1 * anx[-2]);
    resDst /= 16;
    
    uint8_t strength = sqrt(resDst*resDst+resPhi*resPhi);
    //std::cout << (int)resPhi << std::endl;
    if(strength > 50)
    {
      double direction = atan2(resDst,resPhi);
      //std::cout << (int)resPhi << ", ";
      peaks.push_back(SonarPeak(actual2.bearing, (act-actual2.beam.begin()), strength, base::Angle::fromRad(direction)));
    }
    
  }
    
  preprevious2 = previous2;
  previous2 = actual2;
  actual2 = next2;
  next2 = afternext2;
  
  return peaks;
}

  
} //End namespace