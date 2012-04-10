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

std::vector< Line > Line::selectLines(std::vector< Line > lines, std::pair<int,int> basinSize, double angleTolerance, bool alignLines, bool guessMissing)
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
    parallels[i] = selectByDistance(parallels[i], basinSize);
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
    //guess the missing lines to always have 2x2 lines for the output
    if(parallels[bestIdxA].size()+parallels[bestIdxB].size() == 3)
    {
      int fullIdx, halfIdx;
      if(parallels[bestIdxA].size() == 2 && parallels[bestIdxB].size() == 1)
      {
	fullIdx = bestIdxA;
	halfIdx = bestIdxB;
      }
      else if(parallels[bestIdxB].size() == 2 && parallels[bestIdxA].size() == 1)
      {
	fullIdx = bestIdxB;
	halfIdx = bestIdxA;
      }
      int fullDiff = abs(parallels[fullIdx].at(0).d - parallels[fullIdx].at(1).d);
      int halfDiff;
      if(abs(fullDiff - basinSize.first) < abs(fullDiff - basinSize.second))
	halfDiff = basinSize.second;
      else
	halfDiff = basinSize.first;
      
      int newDiff;
      if(parallels[halfIdx].at(0).d < 0)
	newDiff = parallels[halfIdx].at(0).d + halfDiff;
      else
	newDiff = parallels[halfIdx].at(0).d - halfDiff;
      
      parallels[halfIdx].push_back(Line(parallels[halfIdx].at(0).alpha, newDiff, 0));      
    }
  }
  
  lines.insert(lines.end(), parallels[bestIdxA].begin(), parallels[bestIdxA].end());
  lines.insert(lines.end(), parallels[bestIdxB].begin(), parallels[bestIdxB].end());
  
  //print lines
  for(std::vector<Line>::iterator it = lines.begin(); it != lines.end(); it++)
  {
    std::cout << "alpha = " << it->alpha << "[" << it->alpha*180.0/M_PI << "°], d = " << it->d << ", votes: " << it->votes << std::endl;
  }
  
  return lines;
}

std::vector< Line > Line::selectByDistance(std::vector< Line > lines, std::pair<int,int> basinSize)
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
	int basinSizePart = k==0?basinSize.first:basinSize.second;
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