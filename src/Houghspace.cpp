#include "Houghspace.hpp"

namespace sonar_wall_hough
{
  
Houghspace::Houghspace(const Config &config)
  : config(config)
  , space(0)
{
  numberOfAngles = ceil(180 / config.sensorAngularResolution / config.anglesPerBin);
  numberOfPosDistances = (config.maxDistance + config.distancesPerBin - 1) / config.distancesPerBin;
  space = new boost::uint8_t [numberOfAngles*(2*numberOfPosDistances+1)];
  space += numberOfPosDistances * numberOfAngles;
}

Houghspace::~Houghspace()
{
  delete[] (space - (numberOfPosDistances * numberOfAngles));
}

/*
boost::uint8_t* Houghspace::at(Line line)
{
  return at(line.alpha, line.d);
}
*/

boost::uint8_t* Houghspace::at(double angle, int dst)
{
  int angleIdx = angle2Idx(angle);
  int dstIdx = dst2Idx(dst);
  return at(angleIdx, dstIdx);
}

boost::uint8_t* Houghspace::at(int angleIdx, int dstIdx)
{
  while(angleIdx < 0)
  {
    angleIdx += numberOfAngles;
    dstIdx *= -1;
  }
  while(angleIdx >= numberOfAngles)
  {
    angleIdx -= numberOfAngles;
    dstIdx *= -1;
  }
  //does not work, modulo in c++ is strange...
  /*
  angleIdx = angleIdx % numberOfAngles;
  if(angleIdx < 0)
  {
    angleIdx += numberOfAngles;
    dstIdx *= -1;
  }

  */
  if(dstIdx < -numberOfPosDistances || dstIdx > numberOfPosDistances)
  {
    return 0;
  }
  else
    return uncheckedAt(angleIdx, dstIdx);
}

boost::uint8_t* Houghspace::uncheckedAt(int angleIdx, int dstIdx)
{
  return space + (dstIdx * numberOfAngles + angleIdx);
}

void Houghspace::clear()
{
  for(int dstIdx = -numberOfPosDistances; dstIdx <= numberOfPosDistances; dstIdx++)
  {
    for(int angleIdx = 0; angleIdx < numberOfAngles; angleIdx++)
    {
      *at(angleIdx, dstIdx) = 0;
    }
  }
}

int Houghspace::getWidth()
{
  return numberOfAngles;
}

int Houghspace::getHeight()
{
  return 2*numberOfPosDistances+1;
}

int Houghspace::getNumberOfAngles()
{
  return numberOfAngles;
}

int Houghspace::getnumberOfPosDistances()
{
  return numberOfPosDistances;
}

}//end namespace