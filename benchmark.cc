// -*- C++ -*-
// SpatialDataStructuresBenchmark.cc
// An example to compare the efficiency of different spatial data structures

// #include "CommonDefinitions.h"

// c junk
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

// c++ junk
#include <array>
#include <vector>
#include <list>
#include <string>
#include <algorithm>
#include <chrono>
#include <random>
#include <fstream>

using std::string;
using std::vector;
using std::array;
using std::chrono::high_resolution_clock;
using std::chrono::duration;
using std::chrono::duration_cast;

#include "vorocode.hpp"
#include "voro++.hh"
#include <cstdlib>
#include <set>
#include <assert.h>
#include "malloc_count-0.7/malloc_count.h"

namespace Utilities {

// This function returns a random floating point number between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

void
verifyThatDirectoryExists(const std::string & path) {
  std::ifstream test(path);
  if ((bool)test == false) {
    fprintf(stderr, "Error, cannot find directory at %s.  "
            "Please make it yourself (\"mkdir %s\")\n",
            path.c_str(), path.c_str());
    exit(1);
  }
}

size_t
interpolateNumberLinearlyOnLogScale(const size_t lower,
                                    const size_t upper,
                                    const unsigned int numberOfPoints,
                                    const unsigned int pointIndex) {
  const double percent =
    pointIndex / double(numberOfPoints - 1);
  const double power = std::log10(lower) +
    percent * (std::log10(upper) - std::log10(lower));
  return std::pow(10., power);
}

}

// These files contain the functions which run the different spatial data
//  structures.
// #include "Versions_stlib.h"
// #include "Versions_pcl.h"
// #include "Versions_kdtree2.h"

// This file contains the test scenario generators which form the distributions
//  of points.
// #include "PointScenarioGenerators.h"

// This is a result checking function to make sure we have the right answer.
// void
// checkResult(const vector<unsigned int> & correctResult,
//             const vector<unsigned int> & correctResultDelimiters,
//             const vector<unsigned int> & testResult,
//             const vector<unsigned int> & testResultDelimiters,
//             const std::string & testName) {
//   vector<unsigned int> sortedCopyOfCorrectResult = correctResult;
//   vector<unsigned int> sortedCopyOfTestResult = testResult;
//   char sprintfBuffer[500];
//   if (correctResultDelimiters.size() != testResultDelimiters.size()) {
//     sprintf(sprintfBuffer, "wrong result, "
//             "delimiter sizes don't match: correct is %zu, test is %zu, "
//             "test named " BOLD_ON FG_RED "%s" RESET "\n",
//             correctResultDelimiters.size(), testResultDelimiters.size(),
//             testName.c_str());
//     throw std::runtime_error(sprintfBuffer);
//   }
//   const unsigned int numberOfPoints = correctResultDelimiters.size() - 1;
//   for (unsigned int pointIndex = 0;
//        pointIndex < numberOfPoints; ++pointIndex) {
//     try {
//       const unsigned int beginIndex = correctResultDelimiters[pointIndex];
//       const unsigned int endIndex = correctResultDelimiters[pointIndex + 1];
//       if (testResultDelimiters[pointIndex] != beginIndex) {
//         sprintf(sprintfBuffer,
//                 "beginIndex for point %u doesn't match: correct is %u, "
//                 "test is %u\n",
//                 pointIndex, beginIndex, testResultDelimiters[pointIndex]);
//         throw std::runtime_error(sprintfBuffer);
//       }
//       if (testResultDelimiters[pointIndex + 1] != endIndex) {
//         sprintf(sprintfBuffer,
//                 "endIndex for point %u doesn't match: correct is %u, "
//                 "test is %u\n",
//                 pointIndex, endIndex, testResultDelimiters[pointIndex + 1]);
//         throw std::runtime_error(sprintfBuffer);
//       }

//       // sort the two neighborhoods
//       std::sort(sortedCopyOfCorrectResult.begin() + beginIndex,
//                 sortedCopyOfCorrectResult.begin() + endIndex);
//       std::sort(sortedCopyOfTestResult.begin() + beginIndex,
//                 sortedCopyOfTestResult.begin() + endIndex);
//       // check if the two neighborhoods are the same
//       bool neighborhoodsAreTheSame = true;
//       for (unsigned int index = beginIndex; index < endIndex; ++index) {
//         if (sortedCopyOfCorrectResult[index] != sortedCopyOfTestResult[index]) {
//           neighborhoodsAreTheSame = false;
//           fprintf(stderr, "point %7u neighbor %3u is different: correct %7u "
//                   "test %5u\n",
//                   pointIndex, index - beginIndex,
//                   sortedCopyOfCorrectResult[index],
//                   sortedCopyOfTestResult[index]);
//         }
//       }
//       if (neighborhoodsAreTheSame == false) {
//         sprintf(sprintfBuffer, "wrong result, "
//                 "point %u's neighborhood doesn't match\n",
//                 pointIndex);
//         throw std::runtime_error(sprintfBuffer);
//       }
//     } catch (const std::exception & e) {
//       fprintf(stderr, "wrong result, "
//               "point %u had an error, "
//               "test named " BOLD_ON FG_RED "%s" RESET "\n",
//               pointIndex, testName.c_str());
//       throw;
//     }
//   }
// }

// template <class Function>
std::multiset<long>
runTest(const unsigned int numberOfTrials,
        //const Function function,
        vector<Particle> & particles,
        // const BoundingBox & pointsBoundingBox,
        std::vector<double> searchDist,
        double boxLength,
        vector<unsigned int> * const result,
        vector<unsigned int> * const resultDelimiters,
        double * const initializationTime,
        double * const queryingTime,
        string voroVersion,
        double scale,
        double prec
        ) {

  *initializationTime = std::numeric_limits<double>::max();
  *queryingTime = std::numeric_limits<double>::max();

  bool firstRunDone = false;
  std::multiset<long> vols;

  for (unsigned int trialNumber = 0;
       trialNumber < numberOfTrials; ++trialNumber, firstRunDone=true) {
    double thisTrialsInitializationTime = std::numeric_limits<double>::max();
    double thisTrialsQueryingTime = std::numeric_limits<double>::max();

// Randomly add particles into the container
    if (voroVersion == "voro++") {
      resetNumberOfTimesMallocHasBeenCalled();
      const int n_x=6,n_y=6,n_z=6;

      voro::container con(-boxLength/2, boxLength/2, -boxLength/2, boxLength/2, -boxLength/2, boxLength/2, 
      n_x, n_y, n_z, false, false, false, 8);

      for(unsigned int i=0;i<particles.size();i++) {
          con.put(particles[i].id, particles[i].position.X, particles[i].position.Y, particles[i].position.Z);
      }
      
      voro::voronoicell c;
      double vvol = 0;
      voro::c_loop_all vl(con);
      if (vl.start()) do {
        if (con.compute_cell(c,vl)) {
          double vol = c.volume();
          if (firstRunDone) {
            if (vols.find(vol*prec) == vols.end()) {
              std::cout << "ERROR, new volume: " << vol << std::endl;
            }
          }
          else {
            vols.insert(vol*prec);
          }
          // std::cout << "Cell volume: " << vol << std::endl;
          vvol += vol;
        }
      } while (vl.inc());
      // Sum up the volumes, and check that this matches the container volume
      // double vvol=con.sum_cell_volumes();
      // std::cout << "voro++ volume: " << vvol << std::endl;
      // std::cout << "Malloc called: " << getNumberOfTimesMallocHasBeenCalled() << " times." << std::endl;
    }

    if (voroVersion == "voro--") {
      resetNumberOfTimesMallocHasBeenCalled();
      cellContainer con(std::vector<Particle>(), 1/*(scale != 0) ? (boxLength/pow(particles.size(),1.0/3))*scale : boxLength*/, -boxLength/2, boxLength/2, -boxLength/2, boxLength/2, -boxLength/2, boxLength/2);
      for(unsigned int i=0;i<particles.size();i++) {
          con.put(particles[i].id, particles[i].position.X, particles[i].position.Y, particles[i].position.Z);
      }


      double vvol = 0;
      con.sds.initialize(&con.particles[0], &con.particles[con.particles.size()]);
      
      voronoiCell c;
      for(unsigned int i = 0; i < con.particles.size(); ++i) {
        // std::cout << "Search distance for particle " << i << ": " << searchDist[i] << std::endl;
        con.defaultLength = searchDist[i];
        // std::cout << "Mallocs before cell " << i << ": " << getNumberOfTimesMallocHasBeenCalled() << std::endl;
        con.makeCell(i, c);
        // std::cout << "Mallocs after cell " << i << ": " << getNumberOfTimesMallocHasBeenCalled() << std::endl;
        // std::cout << "Cell " << i << " data: " << std::endl;
        // std::cout << "    edges: " <<  c.edges.capacity() << ", " << c.edges.size() << ", " << c.edges.computeOccupancy() << std::endl;
        // std::cout << "    vertices: " <<  c.vertices.capacity() << ", " << c.vertices.size() << ", " << c.vertices.computeOccupancy() << std::endl;
        double vol = c.volume();
        if (firstRunDone) {
          if (vols.find(vol*prec) == vols.end()) {
            std::cout << "ERROR, new volume: " << vol << std::endl;
          }
        }
        else {
          vols.insert(vol*prec);
        }
        vvol += vol;
      }

      // Sum up the volumes, and check that this matches the container volume
      // double vvol=con.sum_cell_volumes();
      // std::cout << "voro-- volume: " << vvol << std::endl;
      // std::cout << "Malloc called: " << getNumberOfTimesMallocHasBeenCalled() << " times." << std::endl;
      // std::cout << "Diagram memory: " << con.get_memory_usage() << std::endl;
    }

    result->resize(0);
    resultDelimiters->resize(0);
    // Do the test
    // function(particles, neighborSearchDistance,
    //          result, resultDelimiters,
    //          &thisTrialsInitializationTime, &thisTrialsQueryingTime);
    // Take the minimum values from all trials
    *initializationTime =
      std::min(*initializationTime, thisTrialsInitializationTime);
    *queryingTime =
      std::min(*queryingTime, thisTrialsQueryingTime);

  }
  return vols;
      }

// double findMaxDist(vector<Particle> particles) {
//   double maxDist = 0;
//   for(particle1 : particles) {
//     for(particle2 : particles) {
//       if (particle1.position == particle2.position) {
//         double dist = sqrt((particle1.position.X-particle2.position.X)^2 +
//          (particle1.position.Y-particle2.position.Y)^2 +
//          (particle1.position.Z-particle2.position.Z)^2);
//         if dist > maxDist {
//           maxDist = dist;
//         }
//       }
//     }
//   }
//   return maxDist;
// }

int main(int argc, char* argv[]) {
  // ignoreUnusedVariable(argc);
  // ignoreUnusedVariable(argv);

  // ===========================================================================
  // *************************** < Inputs> *************************************
  // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

  const array<double, 2> numberOfPointsRange = {{1e3, 1e8}};
  const unsigned int numberOfDataPoints      = 15;
  const unsigned int numberOfTrialsPerSize   = 3;
  const double neighborSearchDistance        = 1.0;

  unsigned int numberOfPoints = 1000;
  double scale = 1.00001;
  double prec = 1e8;
  if (argc >= 2) {
    numberOfPoints = atoi(argv[1]);
  }
  if (argc >= 3) {
    scale = atof(argv[2]);
  }
  if (argc >= 4) {
    prec = atof(argv[3]);
  }
  // Utilities::interpolateNumberLinearlyOnLogScale(numberOfPointsRange[0],
  //                                                numberOfPointsRange[1],
  //                                                numberOfDataPoints,
  //                                                dataPointIndex);
  

  // typedef PointScenarioGenerators::UniformRandomWithAverageNumberOfNeighbors UniformPointGenerator;
  // typedef PointScenarioGenerators::NonUniformRandomWithAverageNumberOfNeighbors NonUniformPointGenerator;

  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  // *************************** </Inputs> *************************************
  // ===========================================================================

  char sprintfBuffer[500];

  const string prefix = "data/benchmark_";
  const string suffix = "_knuth";

  // Make sure that the data directory exists.
  Utilities::verifyThatDirectoryExists("data");

  // Open output files
  // sprintf(sprintfBuffer, "%s%s_stlib_cellArray_results%s.csv",
  //         prefix.c_str(), UniformPointGenerator::getName().c_str(), suffix.c_str());
  // FILE * stlibCellArrayFileUniform = fopen(sprintfBuffer, "w");
  // sprintf(sprintfBuffer, "%s%s_stlib_octTree_results%s.csv",
  //         prefix.c_str(), UniformPointGenerator::getName().c_str(), suffix.c_str());
  // FILE * stlibOctTreeFileUniform = fopen(sprintfBuffer, "w");
  // sprintf(sprintfBuffer, "%s%s_stlib_kdTree_results%s.csv",
  //         prefix.c_str(), UniformPointGenerator::getName().c_str(), suffix.c_str());
  // FILE * stlibKdTreeFileUniform = fopen(sprintfBuffer, "w");
  // sprintf(sprintfBuffer, "%s%s_pcl_kdTree_results%s.csv",
  //         prefix.c_str(), UniformPointGenerator::getName().c_str(), suffix.c_str());
  // FILE * pclKdTreeFileUniform = fopen(sprintfBuffer, "w");
  // sprintf(sprintfBuffer, "%s%s_pcl_ocTree_results%s.csv",
  //         prefix.c_str(), UniformPointGenerator::getName().c_str(), suffix.c_str());
  // FILE * pclOctreeFileUniform = fopen(sprintfBuffer, "w");
  // sprintf(sprintfBuffer, "%s%s_kdtree2_results%s.csv",
  //         prefix.c_str(), UniformPointGenerator::getName().c_str(), suffix.c_str());
  // FILE * kdtree2FileUniform = fopen(sprintfBuffer, "w");

  // sprintf(sprintfBuffer, "%s%s_stlib_cellArray_results%s.csv",
  //         prefix.c_str(), NonUniformPointGenerator::getName().c_str(), suffix.c_str());
  // FILE * stlibCellArrayFileNonUniform = fopen(sprintfBuffer, "w");
  // // sprintf(sprintfBuffer, "%s%s_stlib_octTree_results%s.csv",
  // //         prefix.c_str(), NonUniformPointGenerator::getName().c_str(), suffix.c_str());
  // // FILE * stlibOctTreeFileNonUniform = fopen(sprintfBuffer, "w");
  // sprintf(sprintfBuffer, "%s%s_stlib_kdTree_results%s.csv",
  //         prefix.c_str(), NonUniformPointGenerator::getName().c_str(), suffix.c_str());
  // FILE * stlibKdTreeFileNonUniform = fopen(sprintfBuffer, "w");
  // sprintf(sprintfBuffer, "%s%s_pcl_kdTree_results%s.csv",
  //         prefix.c_str(), NonUniformPointGenerator::getName().c_str(), suffix.c_str());
  // FILE * pclKdTreeFileNonUniform = fopen(sprintfBuffer, "w");
  // sprintf(sprintfBuffer, "%s%s_pcl_ocTree_results%s.csv",
  //         prefix.c_str(), NonUniformPointGenerator::getName().c_str(), suffix.c_str());
  // FILE * pclOctreeFileNonUniform = fopen(sprintfBuffer, "w");
  // sprintf(sprintfBuffer, "%s%s_kdtree2_results%s.csv",
  //         prefix.c_str(), NonUniformPointGenerator::getName().c_str(), suffix.c_str());
  // FILE * kdtree2FileNonUniform = fopen(sprintfBuffer, "w");

  sprintf(sprintfBuffer, "%s_sanic_results%s.csv",
          prefix.c_str(), suffix.c_str());
  FILE * sanic = fopen(sprintfBuffer, "w");

  srand(1);
  // for(int i = 0;i < 4000000; i++) {
  //   Utilities::rnd();
  //   Utilities::rnd();
  //   Utilities::rnd();
  // }
  // For each size
  for (unsigned int dataPointIndex = 0;
       dataPointIndex < numberOfDataPoints;
       ++dataPointIndex) {

    high_resolution_clock::time_point thisSizesTic =
      high_resolution_clock::now();

    // generate points
    vector<Particle> pointsUniform;
    vector<Particle> pointsNonUniform;

    // Randomly add particles into the container
    for(unsigned int i=0;i<numberOfPoints;i++) {
      double x, y, z;
        x=-neighborSearchDistance/2+Utilities::rnd()*(neighborSearchDistance);
        y=-neighborSearchDistance/2+Utilities::rnd()*(neighborSearchDistance);
        z=-neighborSearchDistance/2+Utilities::rnd()*(neighborSearchDistance);
        pointsUniform.push_back(Particle(i,x,y,z,i));
    }

    // if (numberOfPoints <= 1000) {
      cellContainer con(std::vector<Particle>(), neighborSearchDistance, -neighborSearchDistance/2, neighborSearchDistance/2, -neighborSearchDistance/2, neighborSearchDistance/2, -neighborSearchDistance/2, neighborSearchDistance/2);

      for(unsigned int i=0;i<pointsUniform.size();i++) {
        con.put(pointsUniform[i].id, pointsUniform[i].position.X, pointsUniform[i].position.Y, pointsUniform[i].position.Z);
      }
      con.initialize();
      std::vector<double> searchDist = con.findMaxNeighDist(scale);
      double max = *std::max_element(std::begin(searchDist),std::end(searchDist));
      double min = *std::min_element(searchDist.begin(),searchDist.end());

      std::cout << "Difference between max (" << max << ") and min (" << min << ") search distance: " << max-min << std::endl;
    // }
    // // else {
    // //   searchDist = (neighborSearchDistance/pow(numberOfPoints, 1/3.0))*5;
    // // }
    // std::cout << "Max distance: " << findMaxDist(pointsUniform) << std::endl;

    // const unsigned int averageNumberOfNeighborsPerPoint = 70;
    // UniformPointGenerator::generatePoints(numberOfPoints,
    //                                averageNumberOfNeighborsPerPoint,
    //                                neighborSearchDistance, &pointsUniform);

    // set number of normal distributions to use for non uniform distribution
    // unsigned int numGroups = 3;

    // NonUniformPointGenerator::generatePoints(numberOfPoints,
    //                                averageNumberOfNeighborsPerPoint,
    //                                numGroups,
    //                                neighborSearchDistance, &pointsNonUniform);

    // form the bounding box
    // BoundingBox temp;
    // for (const Point & p : pointsUniform) {
    //   temp += p;
    // }
    // stlib::geom::offset(&temp, 1e-3*(temp.upper[0] - temp.lower[0]));
    // const BoundingBox pointsBoundingBox = temp;

    vector<unsigned int> result;
    vector<unsigned int> resultDelimiters;
    double initializationTime;
    double queryingTime;

    // =========================================================================
    // ************************** < do voro++ > ********************************
    // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    thisSizesTic = high_resolution_clock::now();

    std::multiset<long> voropp = runTest(numberOfTrialsPerSize,
            pointsUniform,
            searchDist,
            neighborSearchDistance,
            &result,
            &resultDelimiters,
            &initializationTime,
            &queryingTime,
            "voro++",
            scale, prec);
    fprintf(sanic, "%10.4e, %10.4e, %10.4e\n",
            double(numberOfPoints), initializationTime, queryingTime);
    std::cout<<"finished voro++"<<std::endl;

    high_resolution_clock::time_point thisSizesToc =
      high_resolution_clock::now();
    double thisSizesElapsedTime =
      duration_cast<duration<double> >(thisSizesToc - thisSizesTic).count();
    printf("finished %8.2e points in %6.2f seconds\n", double(numberOfPoints),
           thisSizesElapsedTime);


    // =========================================================================
    // ************************** < do voro-- > ********************************
    // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    thisSizesTic = high_resolution_clock::now();

    std::multiset<long> voromm = runTest(numberOfTrialsPerSize,
            pointsUniform,
            searchDist,
            neighborSearchDistance,
            &result,
            &resultDelimiters,
            &initializationTime,
            &queryingTime,
            "voro--",
            scale, prec);
    fprintf(sanic, "%10.4e, %10.4e, %10.4e\n",
            double(numberOfPoints), initializationTime, queryingTime);
    std::cout<<"finished voro--"<<std::endl;

    /*high_resolution_clock::time_point*/ thisSizesToc =
      high_resolution_clock::now();
    /*double*/ thisSizesElapsedTime =
      duration_cast<duration<double> >(thisSizesToc - thisSizesTic).count();
    printf("finished %8.2e points in %6.2f seconds\n", double(numberOfPoints),
           thisSizesElapsedTime);

    if (voropp != voromm) {
      std::cout << "SHIT GAIZ NO!!!! Different volumes produced!" << std::endl;
      std::size_t diff = 0;

      std::set<long> seen;

      for (auto x : voromm) {
        if (seen.find(x) != seen.end()) {
          continue;
        }
        seen.insert(x);

        std::size_t mm = voromm.count(x);
        std::size_t pp;
        if (voropp.find(x) == voropp.end()) {
          pp = 0;
        }
        else {
          pp = voropp.count(x);
        }

        if (pp < mm) {
          diff += (mm - pp);
          std::cout << "Different volume: " << (x/prec) << ", counts " << mm << " vs. " << pp << std::endl;
        }
      }

      seen.clear();
      for (auto x : voropp) {
        if (seen.find(x) != seen.end()) {
          continue;
        }
        seen.insert(x);

        std::size_t pp = voropp.count(x);
        std::size_t mm;
        if (voromm.find(x) == voromm.end()) {
          mm = 0;
        }
        else {
          mm = voromm.count(x);
        }

        if (mm < pp) {
          std::cout << "Different volume: " << (x/prec) << ", counts " << mm << " vs. " << pp << std::endl;
        }
      }

      std::cout << "Num different: " << diff << std::endl;
      std::cout << "Percentage different: " << diff/numberOfPoints << std::endl;
    }
    else {
      std::cout << "voro++ and voro-- produced the same volumes to prec " << (1/prec) << std::endl;
    }

//     // =========================================================================
//     // ********************** < do stlib cellArray > ***************************
//     // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

//     // uniform
//     runTest(numberOfTrialsPerSize,
//             pointsUniform,
//             pointsBoundingBox,
//             neighborSearchDistance,
//             &result,
//             &resultDelimiters,
//             &initializationTime,
//             &queryingTime);
//     const vector<unsigned int> correctResultU           = result;
//     const vector<unsigned int> correctResultDelimitersU = resultDelimiters;
//     fprintf(stlibCellArrayFileUniform, "%10.4e, %10.4e, %10.4e\n",
//             double(numberOfPoints), initializationTime, queryingTime);
//     std::cout<<"finished uniform cellarray"<<std::endl;

//     // nonuniform
//     runTest(numberOfTrialsPerSize,
//             findNeighbors_stlibCellArray,
//             pointsNonUniform,
//             pointsBoundingBox,
//             neighborSearchDistance,
//             &result,
//             &resultDelimiters,
//             &initializationTime,
//             &queryingTime);

//     const vector<unsigned int> correctResultNU           = result;
//     const vector<unsigned int> correctResultDelimitersNU = resultDelimiters;
//     fprintf(stlibCellArrayFileNonUniform, "%10.4e, %10.4e, %10.4e\n",
//             double(numberOfPoints), initializationTime, queryingTime);
//     std::cout<<"finished nonuniform cellarray"<<std::endl;

//     // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//     // ********************** < do stlib cellArray > ***************************
//     // =========================================================================

//     // =========================================================================
//     // ********************** < do stlib octTree > ***************************
//     // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

//     // uniform
//     runTest(numberOfTrialsPerSize,
//             findNeighbors_stlibOctTree,
//             pointsUniform,
//             pointsBoundingBox,
//             neighborSearchDistance,
//             &result,
//             &resultDelimiters,
//             &initializationTime,
//             &queryingTime);
//     checkResult(correctResultU,
//                 correctResultDelimitersU,
//                 result,
//                 resultDelimiters,
//                 string("stlib octTree"));
//     fprintf(stlibOctTreeFileUniform, "%10.4e, %10.4e, %10.4e\n",
//             double(numberOfPoints), initializationTime, queryingTime);
//     std::cout<<"finished stlib octree"<<std::endl;

//     // // nonuniform
//     // runTest(numberOfTrialsPerSize,
//     //         findNeighbors_stlibOctTree,
//     //         pointsNonUniform,
//     //         pointsBoundingBox,
//     //         neighborSearchDistance,
//     //         &result,
//     //         &resultDelimiters,
//     //         &initializationTime,
//     //         &queryingTime);
//     // checkResult(correctResultNU,
//     //             correctResultDelimitersNU,
//     //             result,
//     //             resultDelimiters,
//     //             string("stlib octTree"));
//     // fprintf(stlibOctTreeFileNonUniform, "%10.4e, %10.4e, %10.4e\n",
//     //         double(numberOfPoints), initializationTime, queryingTime);

//     // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//     // ********************** < do stlib octTree > ***************************
//     // =========================================================================

//     // =========================================================================
//     // ********************** < do stlib kdTree > ***************************
//     // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

//     // uniform
//     runTest(numberOfTrialsPerSize,
//             findNeighbors_stlibKdTree,
//             pointsUniform,
//             pointsBoundingBox,
//             neighborSearchDistance,
//             &result,
//             &resultDelimiters,
//             &initializationTime,
//             &queryingTime);
//     checkResult(correctResultU,
//                 correctResultDelimitersU,
//                 result,
//                 resultDelimiters,
//                 string("stlib kdTree"));
//     fprintf(stlibKdTreeFileUniform, "%10.4e, %10.4e, %10.4e\n",
//             double(numberOfPoints), initializationTime, queryingTime);

//     // nonuniform
//     runTest(numberOfTrialsPerSize,
//             findNeighbors_stlibKdTree,
//             pointsNonUniform,
//             pointsBoundingBox,
//             neighborSearchDistance,
//             &result,
//             &resultDelimiters,
//             &initializationTime,
//             &queryingTime);
//     checkResult(correctResultNU,
//                 correctResultDelimitersNU,
//                 result,
//                 resultDelimiters,
//                 string("stlib kdTree"));
//     fprintf(stlibKdTreeFileNonUniform, "%10.4e, %10.4e, %10.4e\n",
//             double(numberOfPoints), initializationTime, queryingTime);
//     std::cout<<"finished stlib kdtree"<<std::endl;

//     // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//     // ********************** < do stlib kdTree > ***************************
//     // =========================================================================

//     // =========================================================================
//     // ********************** < do pcl kdTree > ***************************
//     // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

//     // uniform
//     runTest(numberOfTrialsPerSize,
//             findNeighbors_pclKdTree,
//             pointsUniform,
//             pointsBoundingBox,
//             neighborSearchDistance,
//             &result,
//             &resultDelimiters,
//             &initializationTime,
//             &queryingTime);
//     checkResult(correctResultU,
//                 correctResultDelimitersU,
//                 result,
//                 resultDelimiters,
//                 string("pcl kdTree"));
//     fprintf(pclKdTreeFileUniform, "%10.4e, %10.4e, %10.4e\n",
//             double(numberOfPoints), initializationTime, queryingTime);

//     // nonuniform
//     runTest(numberOfTrialsPerSize,
//             findNeighbors_pclKdTree,
//             pointsNonUniform,
//             pointsBoundingBox,
//             neighborSearchDistance,
//             &result,
//             &resultDelimiters,
//             &initializationTime,
//             &queryingTime);
//     checkResult(correctResultNU,
//                 correctResultDelimitersNU,
//                 result,
//                 resultDelimiters,
//                 string("pcl kdTree"));
//     fprintf(pclKdTreeFileNonUniform, "%10.4e, %10.4e, %10.4e\n",
//             double(numberOfPoints), initializationTime, queryingTime);
//     std::cout<<"finished pcl kdtree"<<std::endl;

//     // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//     // ********************** < do pcl kdTree > ***************************
//     // =========================================================================

//     // =========================================================================
//     // ********************** < do pcl octree > ***************************
//     // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

//     // uniform
//     runTest(numberOfTrialsPerSize,
//             findNeighbors_pclOctree,
//             pointsUniform,
//             pointsBoundingBox,
//             neighborSearchDistance,
//             &result,
//             &resultDelimiters,
//             &initializationTime,
//             &queryingTime);
//     checkResult(correctResultU,
//                 correctResultDelimitersU,
//                 result,
//                 resultDelimiters,
//                 string("pcl octree"));
//     fprintf(pclOctreeFileUniform, "%10.4e, %10.4e, %10.4e\n",
//             double(numberOfPoints), initializationTime, queryingTime);

//     // nonuniform
//     runTest(numberOfTrialsPerSize,
//             findNeighbors_pclOctree,
//             pointsNonUniform,
//             pointsBoundingBox,
//             neighborSearchDistance,
//             &result,
//             &resultDelimiters,
//             &initializationTime,
//             &queryingTime);
//     checkResult(correctResultNU,
//                 correctResultDelimitersNU,
//                 result,
//                 resultDelimiters,
//                 string("pcl octree"));
//     fprintf(pclOctreeFileNonUniform, "%10.4e, %10.4e, %10.4e\n",
//             double(numberOfPoints), initializationTime, queryingTime);
//     std::cout<<"finished pcl octree"<<std::endl;

//     // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//     // ********************** < do pcl octree > ***************************
//     // =========================================================================

//     // =========================================================================
//     // ********************** < do kdtree2 > ***************************
//     // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

//     // uniform
//     runTest(numberOfTrialsPerSize,
//             findNeighbors_kdtree2,
//             pointsUniform,
//             pointsBoundingBox,
//             neighborSearchDistance,
//             &result,
//             &resultDelimiters,
//             &initializationTime,
//             &queryingTime);
//     checkResult(correctResultU,
//                 correctResultDelimitersU,
//                 result,
//                 resultDelimiters,
//                 string("kdtree2"));
//     fprintf(kdtree2FileUniform, "%10.4e, %10.4e, %10.4e\n",
//             double(numberOfPoints), initializationTime, queryingTime);

//     // nonuniform
//         runTest(numberOfTrialsPerSize,
//             findNeighbors_kdtree2,
//             pointsNonUniform,
//             pointsBoundingBox,
//             neighborSearchDistance,
//             &result,
//             &resultDelimiters,
//             &initializationTime,
//             &queryingTime);
//     checkResult(correctResultNU,
//                 correctResultDelimitersNU,
//                 result,
//                 resultDelimiters,
//                 string("kdtree2"));
//     fprintf(kdtree2FileNonUniform, "%10.4e, %10.4e, %10.4e\n",
//             double(numberOfPoints), initializationTime, queryingTime);

//     // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//     // ********************** < do kdtree2 > ***************************
//     // =========================================================================

    // const high_resolution_clock::time_point thisSizesToc =
    //   high_resolution_clock::now();
    // const double thisSizesElapsedTime =
    //   duration_cast<duration<double> >(thisSizesToc - thisSizesTic).count();
    // printf("finished %8.2e points in %6.2f seconds\n", double(numberOfPoints),
    //        thisSizesElapsedTime);

//   }


//   fclose(stlibCellArrayFileUniform);
//   fclose(stlibCellArrayFileNonUniform);
//   fclose(stlibOctTreeFileUniform);
//   // fclose(stlibOctTreeFileNonUniform);
//   fclose(stlibKdTreeFileUniform);
//   fclose(stlibKdTreeFileNonUniform);
//   fclose(pclKdTreeFileUniform);
//   fclose(pclKdTreeFileNonUniform);
//   fclose(pclOctreeFileUniform);
//   fclose(pclOctreeFileNonUniform);
//   fclose(kdtree2FileUniform);
//   fclose(kdtree2FileUniform);

//   return 0;
}
fclose(sanic);
return 0;
}
