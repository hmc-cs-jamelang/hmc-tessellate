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

// Standard function beginning for non-group tests
#define NOGROUP_HEADER(n,d)						\
	n = NUM_PARTICLES;							\
	d = Diagram(0, SCALE, 0, SCALE, 0, SCALE);	\
	add_n_particles(d, n);

// Standard function beginning for group tests
#define GROUP_HEADER(n,d,g1,g2,t)				\
	n = NUM_PARTICLES;							\
	d = Diagram(0, SCALE, 0, SCALE, 0, SCALE);	\
	g1 = 0;										\
	g2 = 1;										\
	add_n_particles_group(d, n/2, g1);			\
	add_n_particles_group(d, n-n/2, g2);		\
	t = d.targetGroups(g1);

#define TEST_PASS_OR_FAIL(name,fn)						\
	if (fn) {											\
		cout << name << " test complete" << endl;		\
	}													\
	else {												\
		cout << name << " test failed" << endl;			\
		return false;									\
	}

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
void add_n_particles(Diagram& d, int n)
{
	for (int id = 0; id < n; ++id) {
		d.addParticle(random_double(), random_double(), random_double(), id);
	}
}

/*
 * add_n_particles_group
 *
 * Helper function for tests. Inserts n particles in diagram d assigned to group g.
 */
void add_n_particles_group(Diagram& d, int n, int g)
{
	for (int id = 0; id < n; ++id) {
		d.addParticle(random_double(), random_double(), random_double(), id, g);
	}
}

/*
 * neighbor_test_nogroups
 *
 * Check that each cell can find its neighbors when there are no groups.
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
 * neighbor_test_targetgroups
 *
 * Check that each cell can find its neighbors when using target groups.
 */
bool neighbor_test_targetgroups()
{
	int red, blue, n;
	Diagram d;
	TargetGroup t;
	GROUP_HEADER(n,d,red,blue,t);

	vector<int> neighborList;
	Cell c;
	for (int i = 0; i < n; ++i) {
		neighborList.clear();

		c = d.getCell(i, t);
		c.computeNeighbors(neighborList);

		cout << "Neighbors of particle " << i << " in red group: ";
		for (auto neighbor : neighborList) {
			cout << neighbor << " ";
			if (neighbor >= n/2) {
				cout << "INVALID NEIGHBOR FOUND" << endl;
				return false;
			}
		}
		cout << endl;
	}

	cout << "Comparing neighbors with and without groups: ";
	t = d.targetGroups(red, blue);
	vector<int> neighborListNoGroups;
	Cell cn;
	for (int i = 0; i < n; ++i) {
		neighborList.clear();
		neighborListNoGroups.clear();

		c = d.getCell(i, t);
		cn = d.getCell(i);
		c.computeNeighbors(neighborList);
		cn.computeNeighbors(neighborListNoGroups);

		if (neighborList != neighborListNoGroups) {
			cout << "FAILED" << endl;
			return false;
		}
	}
	cout << "COMPLETE" << endl;

	return true;
}

/*
 * run_nogroup_tests
 *
 * Run the tests that do not use groups.
 */
bool run_nogroup_tests()
{
	TEST_PASS_OR_FAIL("Neighbor",neighbor_test_nogroups());
	TEST_PASS_OR_FAIL("Volume",volume_test_nogroups());
	TEST_PASS_OR_FAIL("Vertices",vertices_test_nogroups());
	TEST_PASS_OR_FAIL("Face areas",face_areas_test_nogroups());
	TEST_PASS_OR_FAIL("Fake cell",fake_cell_test_nogroups());

	return true;
}

/*
 * run_targetgroup_tests
 *
 * Run the tests that use target groups.
 */
bool run_targetgroup_tests()
{
	TEST_PASS_OR_FAIL("Neighbor",neighbor_test_targetgroups());

	return true;
}

/*
 * main
 *
 * Runs the tests.
 */
int main()
{
	TEST_PASS_OR_FAIL("No group",run_nogroup_tests());
	TEST_PASS_OR_FAIL("Target group",run_targetgroup_tests());
	
	return 0;
}
