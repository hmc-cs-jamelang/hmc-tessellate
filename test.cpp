/*
 * test.cpp
 *
 * Tests for the hmc-tessellate interface.
 *
 * Does NOT test for correct results, simply that the interface works.
 */

#include <cstdlib>
#include <vector>
#include "VoronoiDiagram.hpp"

using namespace std;
using namespace voro;

// Determines the size of each dimension in each diagram
const int SCALE = 10;
// The number of particles to test in each diagram
const int NUM_PARTICLES = 10;

#define NOGROUP_HEADER(n,d)						\
	n = NUM_PARTICLES;							\
	d = Diagram(0, SCALE, 0, SCALE, 0, SCALE);	\
	add_n_particles(d, n);

/*
 * random_double
 *
 * Helper function for tests. Generate a random double between 0 and SCALE.
 */
double random_double()
{
	return ((double) rand() / RAND_MAX) * SCALE;
}

/*
 * add_n_particles
 *
 * Helper function for tests. Inserts n particles in diagram d.
 */
bool add_n_particles(Diagram& d, int n)
{
	for (int id = 0; id < n; ++id) {
		d.addParticle(random_double(), random_double(), random_double(), id);
	}

	return true;
}

/*
 * neighbor_test_nogroups
 *
 * Check that each cell can find its neighbors without using groups.
 */
bool neighbor_test_nogroups()
{
	int n;
	Diagram d;
	NOGROUP_HEADER(n, d);

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

/*
 * volume_test_nogroups
 *
 * Check that each cell can find its volume when there are no groups.
 */
bool volume_test_nogroups()
{
	int n;
	Diagram d;
	NOGROUP_HEADER(n,d);

	Cell c;
	for (int i = 0; i < n; ++i) {
		c = d.getCell(i);

		cout << "Volume of cell " << i << ": " << c.computeVolume() << endl;
	}

	return true;
}

/*
 * vertices_test_nogroups
 *
 * Check that each cell can find its vertices when there are no groups.
 */
bool vertices_test_nogroups()
{
	int n;
	Diagram d;
	NOGROUP_HEADER(n,d);

	vector<double> verticesList;
	Cell c;
	for (int i = 0; i < n; ++i) {
		verticesList.clear();

		c = d.getCell(i);
		c.computeVertices(verticesList);

		cout << "Vertices of cell " << i << ": ";
		for (auto vertex : verticesList) {
			cout << vertex << " ";
		}
		cout << endl;
	}

	return true;
}

/*
 * face_areas_test_nogroups
 *
 * Check that each cell can find its face areas when there are no groups.
 */
bool face_areas_test_nogroups()
{
	int n;
	Diagram d;
	NOGROUP_HEADER(n,d);

	vector<double> areaList;
	Cell c;
	for (int i = 0; i < n; ++i) {
		areaList.clear();

		c = d.getCell(i);
		c.computeFaceAreas(areaList);

		cout << "Face areas of cell " << i << ": ";
		for (auto area : areaList) {
			cout << area << " ";
		}
		cout << endl;
	}

	return true;
}

/*
 * fake_cell_test_nogroups
 *
 * Check that a fake cell works correctly when there are no groups.
 */
bool fake_cell_test_nogroups()
{
	int n;
	Diagram d;
	NOGROUP_HEADER(n,d);

	double x = random_double();
	double y = random_double();
	double z = random_double();

	Cell c = d.getCell(x, y, z);

	vector<int> neighborList;
	vector<double> verticesList;
	vector<double> areaList;

	double volume = c.computeVolume();
	c.computeNeighbors(neighborList);
	c.computeVertices(verticesList);
	c.computeFaceAreas(areaList);

	cout << "Computed fake particle at (" << x << ", " << y << ", " << z << ")" << endl;
	cout << "\tVolume: " << volume << endl;
	
	cout << "\tNeighbors: ";
	for (auto neighbor : neighborList) {
		cout << neighbor << " ";
	}
	cout << endl;

	cout << "\tVertices: ";
	for (auto vertex : verticesList) {
		cout << vertex << " ";
	}
	cout << endl;

	cout << "\tFace areas: ";
	for (auto area : areaList) {
		cout << area << " ";
	}
	cout << endl;

	return true;
}

/*
 * main
 *
 * Runs the tests.
 */
int main()
{
	neighbor_test_nogroups();
	cout << "Neighbor test complete" << endl;

	volume_test_nogroups();
	cout << "Volume test complete" << endl;

	vertices_test_nogroups();
	cout << "Vertices test complete" << endl;

	face_areas_test_nogroups();
	cout << "Face areas test complete" << endl;

	fake_cell_test_nogroups();
	cout << "Fake cell test complete" << endl;
	
	return 0;
}
