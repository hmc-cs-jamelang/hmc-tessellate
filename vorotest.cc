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
#include <set>

using std::string;
using std::vector;
using std::array;
using std::chrono::high_resolution_clock;
using std::chrono::duration;
using std::chrono::duration_cast;

#include "vorocode.hpp"
#include "voro++.hh"

// This function returns a random floating point number between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

int main(int argc, char* argv[]){
    // Set the number of particles that are going to be randomly introduced
    int particle=1000;
    if (argc >= 2) {
    particle = atoi(argv[1]);
  }
	std::cout << "Start test" << std::endl;
 // Set up constants for the container geometry
    const double x_min=-0.5,x_max=0.5;
    const double y_min=-0.5,y_max=0.5;
    const double z_min=-0.5,z_max=0.5;
    const double cvol=(x_max-x_min)*(y_max-y_min)*(x_max-x_min);
    // const double RAND_MAX = 1;
    
    // Set up the number of blocks that the container is divided into
    const int n_x=1,n_y=1,n_z=1;
    

    std::cout << "Setting constants" << std::endl;
    
    // Create a container with the geometry given above, and make it
    // non-periodic in each of the three coordinates. Allocate space for
    // eight particles within each computational block
    cellContainer vmcon(std::vector<Particle>(), 2*(x_max-x_min), x_min, x_max, y_min, y_max, z_min, z_max);
    std::cout << "Made voro-- container" << std::endl;
    voro::container vpcon(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,false,false,false,particle);
    std::cout << "Made voro++ container" << std::endl;
    voro::voronoicell_neighbor vpneigh;
    std::cout << "Made voro neighbor" << std::endl;
    std::vector<int> vmneighArray[particle];


    
    // Randomly add particles into the container
    srand(1);
    for(int i=0;i<particle;i++) {
        double x=x_min+rnd()*(x_max-x_min);
        double y=y_min+rnd()*(y_max-y_min);
        double z=z_min+rnd()*(z_max-z_min);
        vmcon.put(i,x,y,z);
        vpcon.put(i,x,y,z);

    }
    std::vector<int> wrongNeighbors;
    std::cout << "Generated particles" << std::endl;
    double vmvol = vmcon.sum_cell_volumes();

    std::vector<int> neighborsvm;
    std::vector<double> volumesvm;
    std::vector<double> volumesvp;
    std::cout << "Finding voro-- neighbors" << std::endl;
    voronoiCell c;
    for (int i = 0; i<particle; ++i) {
    	c = vmcon.makeCell(i);
        c.neighbors(neighborsvm);
        volumesvm.push_back(c.volume());
    	vmneighArray[i] = neighborsvm;
    	neighborsvm.resize(0);
    }
    voro::c_loop_all cl(vpcon);
    std::cout << "Made c loop all" << std::endl;

	if(cl.start()) {
		std::cout << "cl starting" << std::endl;
		do {
            std::set<int> vpset;
            std::set<int> vmset;
            std::vector<int> wrong;
            for(int i = 0; i < vmneighArray[cl.pid()].size(); ++i) {
                vmset.insert(vmneighArray[cl.pid()][i]);
            }
			std::cout << "In the loop" << std::endl;
			if(vpcon.compute_cell(vpneigh,cl)) {
				std::cout << "Cell computed" << std::endl;

                volumesvp.push_back(vpneigh.volume());

 	        	// Gather information about the computed Voronoi cell
                std::vector<int> temp;
 	        	vpneigh.neighbors(temp);
                for(int i = 0; i < temp.size(); ++i){
                    if(temp[i] >= 0) {
                        vpset.insert(temp[i]);
                    }
                }

 	        	std::cout << "Number of voro++ neighbors: " << vpset.size() << std::endl;
 	        	std::cout << "Number of voro-- neighbors: " << vmset.size() << std::endl;
                if(vpset.size() != vmset.size()) {
                    std::cout << "Wrong number of neighbors for particle: " << cl.pid() << std::endl;
                    wrongNeighbors.push_back(cl.pid());
                }

                for(std::set<int>::iterator it = vmset.begin(); it != vmset.end(); ++ it) {
                    if(vpset.find(*it) == vpset.end()){
                        wrong.push_back(*it);
                        std::cout << "Wrong neighbor for particle " << cl.pid() << ": " << *it << std::endl; 
                    }
                }

 	        	std::cout << "Printing voro++ neighbors for: " << cl.pid() << std::endl;
 	        	for (std::set<int>::iterator a = vpset.begin(); a != vpset.end(); ++a) {
 	        		std::cout << *a << ", ";
 	        	}
 	        	std::cout << std::endl;


 	        	std::cout << "Printing voro-- neighbors for: " << cl.pid() << std::endl;
                for (std::set<int>::iterator a = vmset.begin(); a != vmset.end(); ++a) {
 	        		std::cout << *a << ", ";
 	        	}
 	        	std::cout << std::endl;
			}

			std::cout << "End of loop" << std::endl;
 		} while (cl.inc());
	}
    // std::cout << "Getting the neighbors for the first cell at: " << vmcon.cells[0]->particle.position.X << " " << vmcon.cells[0]->particle.position.Y << " " << vmcon.cells[0]->particle.position.Z  << std::endl;


    // vmcon.cells[0]->neighbors(neighborsvm);
    // if (neighborsvm.size() != vpset.size()){
    // 	std::cout << "Wrong number of neighbors" << std::endl;
    // }
    // for (unsigned int i = 0; i < neighborsvm.size(); i++) {
    //     for (unsigned int j = 0; j < vpset.size(); j++) {
    //     	if (neighborsvm[i] == vpset[j]) {
    //     		std::cout << "Same neighbor: " << neighborsvm[i] << std::endl;
    //     	}
    //     }
    //     std::cout << "Next neighbor check" << std::endl;
    // }

    // std::cout << "Max neighbor dist: " << vmcon.findMaxNeighDist() << std::endl;

    std::cout << "Number of wrong neighbor list lengths: " << wrongNeighbors.size() << std::endl;
    for(auto stuff : wrongNeighbors) {
        std::cout << stuff << " ";
    }
    std::cout << std::endl;

    std::cout << "Wrong volumes: " << std::endl;
    for(int i = 0; i < volumesvp.size(); ++i) {
        if(std::abs(volumesvp[i] - volumesvm[i]) > 1e-15) {
            printf("Particle: %d, Voro++ volume: %.17g, Voro-- volume: %.17g\n", i,volumesvp[i],volumesvm[i]);
        }
    }
    
    std::cout << "End full diagram test." << std::endl;
}