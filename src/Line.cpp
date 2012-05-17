#include "Line.hpp"

namespace sonar_wall_hough
{

Line::Line(double alpha, int d, int votes)
  :alpha(alpha)
  ,d(d)
  ,votes(votes)
{
}

bool Line::operator<(const Line& other) const
{
  if(alpha < other.alpha)
    return true;
  else if(alpha == other.alpha)
  {
    if(d < other.d)
    return true;
    else if(d == other.d)
    {
      if(votes < other.votes)
	return true;
      else
	return false;
    }
    else
      return false;
  }
  else
    return false;
}

std::pair< base::Vector3d, base::Vector3d > Line::toCartesian(const sonar_wall_hough::Line& limitA, const sonar_wall_hough::Line& limitB)
{
  if(alpha == limitA.alpha || alpha == limitB.alpha)
    return std::pair<base::Vector3d, base::Vector3d>(base::Vector3d(0.0,0.0,0.0), base::Vector3d(10.0,0.0,0.0));
  
  double xA = (limitA.d-d*sin(limitA.alpha)/sin(alpha))/(cos(limitA.alpha)-cos(alpha)*sin(limitA.alpha)/sin(alpha));
  double yA = (limitA.d-d*cos(limitA.alpha)/cos(alpha))/(sin(limitA.alpha)-sin(alpha)*cos(limitA.alpha)/cos(alpha));
  
  double xB = (limitB.d-d*sin(limitB.alpha)/sin(alpha))/(cos(limitB.alpha)-cos(alpha)*sin(limitB.alpha)/sin(alpha));
  double yB = (limitB.d-d*cos(limitB.alpha)/cos(alpha))/(sin(limitB.alpha)-sin(alpha)*cos(limitB.alpha)/cos(alpha));
  
  /*
  std::cout << "alpha = " << alpha << ", d = " << d << ", a_alpha = " << limitA.alpha << ", a_d = " << limitA.d << std::endl;
  std::cout << "intersection: x = " << xA << ", y = " << yA << std::endl;
  std::cout << "alpha = " << alpha << ", d = " << d << ", b_alpha = " << limitB.alpha << ", b_d = " << limitB.d << std::endl;
  std::cout << "intersection: x = " << xB << ", y = " << yB << std::endl;
  */
  return std::pair<base::Vector3d, base::Vector3d>(base::Vector3d(xA, yA, 0.0), base::Vector3d(xB, yB, 0.0));
}

std::vector< Line > Line::selectLines(std::vector< Line > lines, std::pair<int,int> basinSize, double spatialResolution, double angleTolerance, double orientation, bool alignLines, bool guessMissing)
{
  if(lines.size() <= 1)
    return lines;
  
  //group parallel lines
  std::vector<std::vector<Line> > parallels;
  for(int i = 0; i < lines.size(); i++)
  {
    std::vector<Line> parallel;
    for(int j = 0; j < lines.size(); j++)
    {
      //all lines have alpha [0-M_PI]
      if(fabs(lines[i].alpha - lines[j].alpha) < angleTolerance)
      {
	//the lines are parallel
	parallel.push_back(lines[j]);
      }
    }
    //sort parallel vector
    std::sort(parallel.begin(), parallel.end());
    
    parallels.push_back(parallel);
  }
  
  lines.clear();
  
  //delete vectors which lines are subset of another vector
  for(int i = 0; i < parallels.size();)
  {
    bool isSubset = false;
    for(int j = 0; j < parallels.size(); j++)
    {
      if(i == j)
	continue;
      
      if(includes(parallels[j].begin(), parallels[j].end(), parallels[i].begin(), parallels[i].end()))
      {
	//i is subset of j
	isSubset = true;
	break;
      }
    }
    if(isSubset)
    {
      parallels.erase(parallels.begin()+i);
    }
    else
    {
      i++;
    }
  }
  //reject parallel lines with wrong distances
  std::vector<double> meanAngle(parallels.size(), 0.0);
  std::vector<int> totalVotes(parallels.size(), 0);
  for(int i = 0; i < parallels.size(); i++)
  {
    parallels[i] = selectByDistance(parallels[i], basinSize, spatialResolution);
    //calculate mean angle and sum of votes for each parallels-vector
    for(int j = 0; j < parallels[i].size(); j++)
    {
      meanAngle[i] += parallels[i].at(j).alpha;
      totalVotes[i] += parallels[i].at(j).votes;
    }
    meanAngle[i] /= parallels[i].size();
  }
  
  //only proceed if there are 2 or more parallels-vectors
  if(parallels.size() < 2)
    return lines;
  
  //select best match of 2 parallel-vectors being perpendicular
  int bestVotes = 0, bestIdxA, bestIdxB;
  for(int a = 0; a < parallels.size(); a++)
  {
    for(int b = a+1; b < parallels.size(); b++)
    {
      if(fabs(fabs(meanAngle[a] - meanAngle[b]) - M_PI/2) < angleTolerance)
      {
	int sumVotes = totalVotes[a]+totalVotes[b];
	if(sumVotes > bestVotes)
	{
	  bestVotes = sumVotes;
	  bestIdxA = a;
	  bestIdxB = b;
	}
      }
    }
  }
  if(bestVotes == 0)
    return lines;
  
  if(alignLines)
  {
    //make parallel lines' angles equal and perpendicular lines 90 degree for easier position estimation
    
    //calculate mean angles
    double meanA = 0.0, meanB = 0.0;
    for(int i = 0; i < parallels[bestIdxA].size(); i++)
    {
      meanA += parallels[bestIdxA].at(i).alpha;
    }
    meanA /= parallels[bestIdxA].size();
    for(int i = 0; i < parallels[bestIdxB].size(); i++)
    {
      meanB += parallels[bestIdxB].at(i).alpha;
    }
    meanB /= parallels[bestIdxB].size();
    
    //make perpendicular
    std::cout << "angles before "<<meanA << " and "<<meanB<<std::endl;
    double diff = fmod(meanA - meanB - M_PI/2, M_PI);
    meanA -= diff/2;
    meanB += diff/2;
    std::cout << "angles are now "<<meanA << " and "<<meanB<<std::endl;
    
    //set new angles
    for(int i = 0; i < parallels[bestIdxA].size(); i++)
    {
      parallels[bestIdxA].at(i).alpha = meanA;
    }
    for(int i = 0; i < parallels[bestIdxB].size(); i++)
    {
      parallels[bestIdxB].at(i).alpha = meanB;
    }
  }
  
  if(guessMissing)
  {
    //guess missing lines
    std::vector<std::vector<Line> > bestPar;
    bestPar.push_back(parallels[bestIdxA]);
    bestPar.push_back(parallels[bestIdxB]);
    for(int i = 0; i < 2; i++)
    {
      std::vector<Line> par = bestPar[i];
      
      if(par.size() < 2)
      {
	//is it more height or width borders of basin? (orientation points into y-direction of basin)
	int neededDistance = 0;
	double diff = fmod(fabs(par[0].alpha - orientation), M_PI);
	if(diff < M_PI/4)
	{
	  neededDistance = basinSize.first / spatialResolution;
	  //std::cout << "More in line with orientation, distance = height: "<<basinSize.first<<" (="<<neededDistance<<")"<<std::endl;
	}
	else
	{
	  neededDistance = basinSize.second / spatialResolution;
	  //std::cout << "More orthogonal to orientation, distance = width: "<<basinSize.second<<" (="<<neededDistance<<")"<<std::endl;
	}
	
	//guess other line (must be on opposite side of auv)
	int dNew = 0;
	if(par[0].d > 0)
	  dNew = par[0].d - neededDistance;
	else
	  dNew = par[0].d + neededDistance;
	
	if(diff < M_PI/4)
	  par.insert(par.begin(), Line(par[0].alpha, dNew, 0));
	else
	  par.push_back(Line(par[0].alpha, dNew, 0));
	
	//std::cout << "added line : alpha = " << par[1].alpha << "[" << par[1].alpha*180.0/M_PI << "°], d = " << par[1].d << ", votes: " << par[1].votes << std::endl;
      }
      lines.insert(lines.end(), par.begin(), par.end());
    }
  }
  
  //print lines
  for(std::vector<Line>::iterator it = lines.begin(); it != lines.end(); it++)
  {
    std::cout << "alpha = " << it->alpha << "[" << it->alpha*180.0/M_PI << "°], d = " << it->d << ", votes: " << it->votes << std::endl;
  }
  
  return lines;
}

std::vector< Line > Line::selectByDistance(std::vector< Line > lines, std::pair<int,int> basinSize, double spatialResolution)
{
  if(lines.size() < 2)
    return lines;
  
  std::vector<double> closestPartner(lines.size(), 1.0);
  double minDist = 1.0;
  int maxVotes = 0;
  for(int i = 0; i < lines.size(); i++)
  {
    if(lines[i].votes > maxVotes)
      maxVotes = lines[i].votes;
    
    for(int j = i+1; j < lines.size(); j++)
    {
      int distance = abs(lines[i].d - lines[j].d);
      for(int k = 0; k < 2; k++)
      {
	int basinSizePart = (k==0?basinSize.first:basinSize.second) / spatialResolution;
	
	double distDiff = abs(distance - basinSizePart) / (double)basinSizePart;
	if(distDiff < 0.1)
	{
	  if(closestPartner[i] > distDiff)
	    closestPartner[i] = distDiff;
	  if(closestPartner[j] > distDiff)
	    closestPartner[j] = distDiff;
	  if(distDiff < minDist)
	    minDist = distDiff;
	  break;
	}
      }
    }
  }
  //only keep lines with lowest distance, or with highest votes if equal
  for(int i = 0; i < lines.size();)
  {
    if(minDist < 1.0)
    {
      if(closestPartner[i] > minDist)
      {
	lines.erase(lines.begin()+i);
	closestPartner.erase(closestPartner.begin()+i);
      }
      else
	i++;	
    }
    else
    {
      if(lines[i].votes < maxVotes)
      {
	lines.erase(lines.begin()+i);
	closestPartner.erase(closestPartner.begin()+i);
      }
      else
	i++;
    }
  }
  
  return lines;
}
  
} //end namespace