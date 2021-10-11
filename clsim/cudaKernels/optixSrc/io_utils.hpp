/*The MIT License (MIT)

Copyright (c) 2020, Hendrik Schwanekamp hschwanekamp@nvidia.com, Ramona Hohl rhohl@nvidia.com

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGSEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <utility>
#include <cstdlib> 

#include <dataStructCuda.cuh>

void photonsToFile(const std::string& filename, I3CUDASimPhoton *photons, unsigned int nphotons){
  std::cout<< "writing "<< nphotons << " to file "<< filename<<std::endl;
  std::ofstream outputFile; outputFile.open (filename);
  for (unsigned int i = 0; i <  nphotons; i++)
  {  
              outputFile <<photons[i].posAndTime.x << "," << photons[i].posAndTime.y << "," << photons[i].posAndTime.z << std::endl;
  }
  outputFile.close();
}

std::vector<std::string> getNextLineAndSplitIntoTokens(std::istream& str, const char seperator)
{
    std::vector<std::string>   result;
    std::string                line;
    std::getline(str,line);
    std::stringstream          lineStream(line);
    std::string                cell;

    while(std::getline(lineStream,cell, seperator))
    {
        result.push_back(cell);
    }
    // This checks for a trailing comma with no data after it.
    if (!lineStream && cell.empty())
    {
        // If there was a trailing comma then add an empty element.
        result.push_back("");
    }
    return result;
}

void readDOMSFile( std::string &fname, int startCol, int endCol, int nCols,  std::vector<float3> &input)
{ 
  std::ifstream tsvFile(fname);
  const char seperator = ','; // '\t';


  if(!tsvFile.is_open())
    throw std::runtime_error("could not open input file");

  std::cout << "reading  from "<< fname<<std::endl;
  getNextLineAndSplitIntoTokens(tsvFile,seperator);

  while( !tsvFile.eof() )
  {
    auto line = getNextLineAndSplitIntoTokens(tsvFile,seperator);
    if(line.size() < nCols)
      continue;
    //for (int i = startCol; i<endCol; ++i)
      input.push_back(make_float3(std::stof(line[startCol]),  std::stof(line[startCol+1]), std::stof(line[startCol+2]) ));
      
 }
}

void getDOMPos(std::vector<float3>& domPos, float& radius, std::string& domfile){
  readDOMSFile( domfile, 0,3, 5,  domPos);
  radius = 0.1651;
  //int i = 0;
  //std::cout<< " first dom "<< domPos[i].x << "," << domPos[i].y << "," <<domPos[i].z <<std::endl;
 
}

inline float dist( float3 a, float3 b){
    return std::sqrt( (a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) + (a.z-b.z)*(a.z-b.z) );
}

void verifyHitsLocations(float3*  hits, int nhits){
    std::string domfile =  DATA_DIR "doms.tsv";
    //get DOM positions
    std::vector<float3> domPos;
    float radius;
    getDOMPos(domPos, radius,domfile);
 
    const float maxDist = 1000000.0;

    float closestDOM[nhits];
    float avrgdist = 0;
    int found = 0;

    //check if feasible hit position
    for (int ihit =0 ; ihit< nhits; ++ihit)
    {
        closestDOM[ihit] = maxDist;
        int domID =-1;
        for(int idom= 0 ; idom < domPos.size(); ++idom)
        {
            float d = std::abs( dist(domPos[idom], hits[ihit] ) - radius);
            if( d < closestDOM[ihit]){
                closestDOM[ihit] = d;
                domID = idom;

            }
        }

        if(domID !=-1 && closestDOM[ihit]< maxDist){
            ++found;
            avrgdist += closestDOM[ihit];
        //    std::cout<< " domPos "<< domPos[domID].x <<" , "<< domPos[domID].y <<" , "<< domPos[domID].z<<std::endl;
          //  std::cout<< " hitPos "<< hits[ihit].second.x <<" , "<< hits[ihit].y <<" , "<< hits[ihit].z<<std::endl;
         //   std::cout<< " dist "<< closestDOM[ihit]<<std::endl;

        }
    }
    avrgdist /= float(found);
    std::cout <<" found "<< found << " closest DOMS out of "<< nhits <<" hit positions, with avrg distance to closest DOM of "<<avrgdist<<   std::endl;
}
 