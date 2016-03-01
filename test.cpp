#include <cstdlib>
#include <vector>
#include "VoronoiDiagram.hpp"

const int SCALE = 10;
const int NUM_PARTICLES = 10;

using namespace std;
using namespace voro;

double random_double()
{
	return ((double) rand() / RAND_MAX) * SCALE;
}

bool add_n_particles(Diagram& d, int n)
{
	for (int id = 0; id < n; ++id) {
		d.addParticle(random_double(), random_double(), random_double(), id);
	}

	return true;
}

bool neighbor_test_nogroups()
{
	int n = NUM_PARTICLES;

	Diagram d = Diagram(0, SCALE, 0, SCALE, 0, SCALE);
	add_n_particles(d, n);

	vector<int> neighborList;
	Cell c;
	for (int i = 0; i < n; ++i) {
		neighborList.clear();

		c = d.getCell(i);
		c.computeNeighbors(neighborList);

		cout << "Neighbor list of particle " << i << ": ";
		for (auto neighbor : neighborList) {
			cout << neighbor << " ";
		}
		cout << endl;
	}

	return true;
}

int main()
{
	neighbor_test_nogroups();
	cout << "Neighbor test complete" << endl;
	
	return 0;
}
