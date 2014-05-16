#include "Line.hpp"

namespace sonar_wall_hough
{

Line::Line(double alpha, int d, int votes)
  :alpha(alpha)
  ,d(d)
  ,votes(votes)
{
}

bool Line::debug = false;

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

std::vector<Line >  Line::selectLines3(std::vector<Line > lines, std::pair< int, int > basinSize, double spatialResolution, double angularTolerance, double basinOrientation, Houghspace& houghspace, double& orientationDrift){
 
  if(debug)
    std::cout << "selecting lines" << std::endl;
  
  std::vector<Line> result;
  
  //find best vertical lines
  std::vector<LinePair> corrVert = findCorrecpondence(lines, basinSize.second/spatialResolution, houghspace);
  
  if(debug)
    std::cout << "Found " << corrVert.size() << " vertical lines" << std::endl;
  
  if(corrVert.size() > 0){        
    
    std::vector<LineQuartet> lines_complete;
    
    for(std::vector<LinePair>::iterator it_pair = corrVert.begin(); it_pair != corrVert.end(); it_pair++){
    
      std::vector<Line> linesHorz;
      for(std::vector<Line>::iterator it = lines.begin(); it != lines.end(); it++){
        
        //is line in 90 degree angle to the vertical line?
        if( std::fabs( (it->alpha - M_PI/2) - it_pair->a.alpha ) < angularTolerance 
          || std::fabs( (it->alpha + M_PI/2) - it_pair->a.alpha ) < angularTolerance){
          linesHorz.push_back(*it);          
        }
        
      }    
    
      std::vector<LinePair> corrHorz = findCorrecpondence(linesHorz, basinSize.first/spatialResolution, houghspace);
      
      for(std::vector<LinePair>::iterator it_horz = corrHorz.begin(); it_horz != corrHorz.end(); it_horz++){
        
        int new_score = (it_pair->score + it_horz->score) / 2.0;
        LineQuartet quartet = LineQuartet(*it_pair, *it_horz, new_score);
        
        lines_complete.push_back(quartet);        
      }
    
    }
    
    if(debug)
      std::cout << "Found " << lines_complete.size() << " line quartets" << std::endl;
    
    //Find best line quartet    
    if(lines_complete.size() > 0){
      
      LineQuartet best = lines_complete.front();
      
      for(std::vector<LineQuartet>::iterator it_quart = lines_complete.begin(); it_quart != lines_complete.end(); it_quart++){
        
        if(it_quart->score > best.score)
          best = *it_quart;        
      }
      
      result.push_back(best.horz.a);
      result.push_back(best.horz.b);
      result.push_back(best.vert.a);
      result.push_back(best.vert.b);
      
      for(std::vector<Line>::iterator it_line = result.begin(); it_line != result.end(); it_line++){
        
        if(it_line->alpha < -M_PI/2.0){
          it_line->alpha += M_PI;
          it_line->d *= -1;
        }
        
        if(it_line->alpha > M_PI/2.0){
          it_line->alpha -= M_PI;
          it_line->d *= -1;
        }
          
      }
      
      
      
      rectifyLines(result, basinOrientation, orientationDrift);
      
    }else{
      return result;
    }    
    
  }  
  
  return result;
}  


std::vector< Line > Line::selectLines2(std::vector< Line > lines, std::pair< int, int > basinSize, double spatialResolution, double angleTolerance, double basinOrientation, Houghspace& houghspace, double& orientationDrift)
{  
  if(debug)
    std::cout << "selecting lines with basin orientation = " << basinOrientation << " and tolerance " << angleTolerance << std::endl;
  
  std::vector<Line> result;
  
  //for each line check if it is along x-axis of basin (horz), y-axis of basin (vert), or none
  std::vector<Line> linesVert, linesHorz;
  for(std::vector<Line>::iterator it = lines.begin(); it != lines.end(); it++)
  {
    //bring line to same angle region as basin
    if((it->alpha - basinOrientation) > M_PI/2) // line angle too big
    {
      while((it->alpha - basinOrientation) > M_PI/2)
      {
	it->alpha -= M_PI;
	it->d *= -1;
      }
    }
    if((basinOrientation - it->alpha) > M_PI/2) // line angle too small
    {
      while((basinOrientation - it->alpha) > M_PI/2)
      {
	it->alpha += M_PI;
	it->d *= -1;
      }
    }
    
    if(fabs(0.0 - (basinOrientation - it->alpha) ) < angleTolerance)
    {
      //only for sauc-e to reject wrong lines where there surely is no wall
      //if(it->d > 0) // dont detect back wall (in direction of open water)
	linesVert.push_back(*it);
    }
    //horizontal lines must have angle about basinOrientation _+_ M_PI/2
    else if(fabs(M_PI/2 - (basinOrientation - it->alpha) ) < angleTolerance)
    {
      it->alpha += M_PI;
      it->d *= -1;
      linesHorz.push_back(*it);
    }
    else if(fabs(-M_PI/2 - (basinOrientation - it->alpha) ) < angleTolerance)
    {
      linesHorz.push_back(*it);
    }
    /*
    //angle difference to vert walls ( = orientation)
    base::Angle diff = base::Angle::fromRad(basinOrientation-it->alpha); //line angles are always between 0 and M_PI, basinOrientation is between -M_PI and M_PI
    if(fabs(0.0 - diff.getRad()) < angleTolerance || fabs(M_PI - diff.getRad()) < angleTolerance || fabs(-M_PI - diff.getRad()) < angleTolerance)
    {
      //line is vertical, align
      if(!(fabs(0.0 - diff.getRad()) < angleTolerance))
      {
	it->d *= -1;
	it->alpha -= M_PI;
      }
      linesVert.push_back(*it);
    }
    else if(fabs(-M_PI/2 - diff.getRad()) < angleTolerance || fabs(M_PI/2 - diff.getRad()) < angleTolerance)
    {
      //line is horizontal, align
      if(fabs(M_PI/2 - diff.getRad()) < angleTolerance)
      {
	it->d *= -1;
	it->alpha -= M_PI;
      }
      
      linesHorz.push_back(*it);
    }
    */
  }
  
  if(debug)
    std::cout << "we have "<<linesVert.size()<<" vertical and " << linesHorz.size()<<" horizontal lines"<<std::endl;
    
  //check if there is a corresponding line or generate such (with votes from houghspace), generate score
  std::vector<LinePair> corrVert = findCorrecpondence(linesVert, basinSize.second/spatialResolution, houghspace);
  std::vector<LinePair> corrHorz = findCorrecpondence(linesHorz, basinSize.first/spatialResolution, houghspace);
  
  if(corrHorz.empty() || corrVert.empty())
  {
    
    if(debug)
      std::cout << "at least one of the horz/vert vectors is empty. returning empty vector" << std::endl;
    
    return result; //returns empty vector meaning "no valid lines found"
  }
  
  //take x and y lines with best score
  
  int bestIdx = 0;
  for(int i = 0; i < corrHorz.size(); i++)
  {
    if(corrHorz[i].score > corrHorz[bestIdx].score)
      bestIdx = i;
  }
  result.push_back(corrHorz[bestIdx].a);
  result.push_back(corrHorz[bestIdx].b);
  
  bestIdx = 0;
  for(int i = 0; i < corrVert.size(); i++)
  {
    if(corrVert[i].score > corrVert[bestIdx].score)
      bestIdx = i;
  }
  result.push_back(corrVert[bestIdx].a);
  result.push_back(corrVert[bestIdx].b);
  
  rectifyLines(result, basinOrientation, orientationDrift);
  
  return result;
}

void Line::rectifyLines(std::vector<Line> &lines, double basinOrientation, double &orientationDrift){

  if(lines.size() < 4)
    return;  
  
  if(debug){
    std::cout << "Rectify lines" << std::endl;
    std::cout << "Horizontal lines: " << std::endl;
    std::cout << "angle: " << lines[0].alpha << " d " << lines[0].d << " votes " << lines[0].votes << std::endl;
    std::cout << "angle: " << lines[1].alpha << " d " << lines[1].d << " votes " << lines[1].votes << std::endl;
    std::cout << "Vertical lines: " << std::endl;
    std::cout << "angle: " << lines[2].alpha << " d " << lines[2].d << " votes " << lines[2].votes << std::endl;
    std::cout << "angle: " << lines[3].alpha << " d " << lines[3].d << " votes " << lines[3].votes << std::endl;
  }
  
  //calculate orientationDrift
  double actualOrientation = 0.0;
  int votes = 0;
  if(lines[0].alpha >= 0.0)
    actualOrientation += (lines[0].alpha - M_PI/2) * lines[0].votes;
  else
    actualOrientation += (lines[0].alpha + M_PI/2) * lines[0].votes;
  votes += lines[0].votes;
  if(lines[1].alpha >= 0.0)
    actualOrientation += (lines[1].alpha - M_PI/2) * lines[1].votes;
  else
    actualOrientation += (lines[1].alpha + M_PI/2) * lines[1].votes;
  votes += lines[1].votes;
  actualOrientation += (lines[2].alpha) * lines[2].votes;
  votes += lines[2].votes;
  actualOrientation += (lines[3].alpha) * lines[3].votes;
  votes += lines[3].votes;
  actualOrientation /= votes;
  orientationDrift = actualOrientation - basinOrientation;
  
  //rectify lines
  //lines[0,1] must be actualOrientation + M_PI/2
  //lines[2,3] must be actualOrientation
  if(fabs(lines[0].alpha - (actualOrientation+M_PI/2)) > M_PI/2)
    lines[0].d *= -1;
  lines[0].alpha = actualOrientation + M_PI/2;
  
  if(fabs(lines[1].alpha - (actualOrientation+M_PI/2)) > M_PI/2)
    lines[1].d *= -1;
  lines[1].alpha = actualOrientation + M_PI/2;
  
  if(fabs(lines[2].alpha - (actualOrientation)) > M_PI/2)
    lines[2].d *= -1;
  lines[2].alpha = actualOrientation;
  
  if(fabs(lines[3].alpha - (actualOrientation)) > M_PI/2)
    lines[3].d *= -1;
  lines[3].alpha = actualOrientation;  
  
}


std::vector<LinePair> Line::findCorrecpondence(std::vector< Line > lines, int distance, Houghspace& houghspace)
{
  //sigma for normal distribution for scoring
  double sigmasq = - (0.01*distance*distance)/(2*log(0.1));
  
  std::vector<LinePair> pairs;
  
  if(debug)
    std::cout << "wanting distance " << distance << std::endl;
  
  for(int i = 0; i < lines.size(); i++)
  {
    
    if(debug)
      std::cout << "looking for best correspondence for line " << lines[i].alpha << "rad/"<<lines[i].d << std::endl;
    
    if(abs(lines[i].d) > distance)
    {
      if(debug)
	std::cout << "Line is too far away. skipping." << std::endl;

      continue;
    }
    
    int score = 0;
    int bestIdx = -1;
    //get best corresponding line
    for(int j = 0; j < lines.size(); j++)
    {
      if(i == j)
	continue;
      
      //lines must lie on opposite sides of origin (i.e. different signs in dst)
      if(lines[i].d<0 == lines[j].d<0)
	continue;
      
      int dstDiff = abs(lines[i].d - lines[j].d) - distance;
      int tscore = (lines[i].votes + lines[j].votes) * exp(-0.5*dstDiff*dstDiff/sigmasq)-0.1;//make score
      
      if(tscore > score)
      {
	bestIdx = j;
	score = tscore;
      }
    }
    //check if new Line with perfect distance would be better
    int pd = lines[i].d > 0 ? lines[i].d - distance : lines[i].d + distance;
    uint8_t* pvotes = houghspace.at(lines[i].alpha, pd);
    int pscore = 0;
    if(pvotes != NULL)
    {
      
      pscore = (*pvotes + lines[i].votes) * exp(-0.5)-0.1; //make score
      if(bestIdx == -1 || pscore >= score)
      {
	bestIdx = -1;
	score = pscore;
      }
    }
    else if(bestIdx == -1) //if invalid in hough space, set to same as lines[i]. This happens when the maxDistance is smaller than the lines dst should be
      score = lines[i].votes;
    
    //push best pair
    if(bestIdx != -1)
    {
      if(debug)
	std::cout << "found " << lines[bestIdx].alpha << "rad/" << lines[bestIdx].d << " -> score = " << score << std::endl;
      
      pairs.push_back(LinePair(lines[i], lines[bestIdx], score));
    }
    else
    {
      if(debug)
	std::cout << "found none, creatin with score = " << score << std::endl;
      
      pairs.push_back(LinePair(lines[i], Line(lines[i].alpha, pd, pscore), score));
    }
    //std::cout << "are they different? " << pairs.back().a.d << " != " << pairs.back().b.d << std::endl;
  }
  
  return pairs;
}

void Line::setDebug(bool d){
  debug = d;
}

} //end namespace