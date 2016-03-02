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

// Macro for printing the elements of a vector, which we do frequently
#define ITERATE_AND_PRINT(vect)					\
	for (auto& element : vect) {				\
		cout << element << " ";					\
	}											\
	cout << endl;

// Macro for comparing two values. Print and return a failure state if they
// are different.
#define FAIL_IF_DIFFERENT(fst,snd) \
	if (fst != snd) {			   \
		cout << "FAILED" << endl;  \
		return false;			   \
	}

// Macro for checking whether Diagram function fn affects groups and non-groups
// with the same particles in the same way.
#define COMPARE_TARGET_NOGROUP(groupVector,nogroupVector,fn)		\
	t = d.targetGroups(red, blue);									\
	Cell cn;														\
	for (int i = 0; i < n; ++i) {									\
		groupVector.clear();										\
		nogroupVector.clear();										\
																	\
		c = d.getCell(i, t);										\
		cn = d.getCell(i);											\
		c.fn(groupVector);											\
		cn.fn(nogroupVector);										\
																	\
		if (groupVector != nogroupVector) {							\
			cout << "FAILED" << endl;								\
			return false;											\
		}															\
	}																\
	cout << "COMPLETE" << endl;

// Macro for running a test and stating whether it passed or failed.
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
		ITERATE_AND_PRINT(neighborList);
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
		ITERATE_AND_PRINT(verticesList);
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
		ITERATE_AND_PRINT(areaList);
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
	ITERATE_AND_PRINT(neighborList);

	cout << "\tVertices: ";
	ITERATE_AND_PRINT(verticesList);

	cout << "\tFace areas: ";
	ITERATE_AND_PRINT(areaList);

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

	vector<int> neighborListNoGroups;
	cout << "Comparing neighbors with and without groups: ";
	COMPARE_TARGET_NOGROUP(neighborList,neighborListNoGroups,computeNeighbors);

	return true;
}

/*
 * vertices_test_targetgroups
 *
 * Check that each cell can find its vertices when using target groups.
 */
bool vertices_test_targetgroups()
{
	int red, blue, n;
	Diagram d;
	TargetGroup t;
	GROUP_HEADER(n,d,red,blue,t);

	vector<double> verticesList;
	Cell c;
	for (int i = 0; i < n; ++i) {
		verticesList.clear();

		c = d.getCell(i, t);
		c.computeVertices(verticesList);

		cout << "Vertices of particle " << i << " in red group: ";
		ITERATE_AND_PRINT(verticesList);
	}

	vector<double> verticesListNoGroups;
	cout << "Comparing vertices with and without groups: ";
	COMPARE_TARGET_NOGROUP(verticesList,verticesListNoGroups,computeVertices);

	return true;
}

/*
 * face_area_test_targetgroups
 *
 * Check that each cell can find its face areas when using target groups.
 */
bool face_area_test_targetgroups()
{
	int red, blue, n;
	Diagram d;
	TargetGroup t;
	GROUP_HEADER(n,d,red,blue,t);

	vector<double> areaList;
	Cell c;
	for (int i = 0; i < n; ++i) {
		areaList.clear();

		c = d.getCell(i, t);
		c.computeFaceAreas(areaList);

		cout << "Face areas of particle " << i << " in red group: ";
		ITERATE_AND_PRINT(areaList);
	}

	vector<double> areaListNoGroups;
	cout << "Comparing face areas with and without groups: ";
	COMPARE_TARGET_NOGROUP(areaList,areaListNoGroups,computeFaceAreas);

	return true;
}

/*
 * volume_test_targetgroups
 *
 * Check that each cell can find its volume when using target groups.
 */
bool volume_test_targetgroups()
{
	int red, blue, n;
	Diagram d;
	TargetGroup t;
	GROUP_HEADER(n,d,red,blue,t);

	double volume;
	Cell c;
	for (int i = 0; i < n; ++i) {
		c = d.getCell(i, t);
		volume = c.computeVolume();

		cout << "Volume of particle " << i << " in red group: " << volume << endl;
	}

	cout << "Comparing volumes with and without groups: ";
	t = d.targetGroups(red, blue);
	Cell cn;
	for (int i = 0; i < n; ++i) {
		c = d.getCell(i, t);
		cn = d.getCell(i);

		FAIL_IF_DIFFERENT(c.computeVolume(),cn.computeVolume());
	}
	cout << "COMPLETE" << endl;

	return true;
}

/*
 * fake_cell_test_targetgroups
 *
 * Check that a fake cell works correctly when using target groups.
 */
bool fake_cell_test_targetgroups()
{
	int red, blue, n;
	Diagram d;
	TargetGroup t;
	GROUP_HEADER(n,d,red,blue,t);

	double x = random_double();
	double y = random_double();
	double z = random_double();

	Cell c = d.getCell(x, y, z, t);

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
	ITERATE_AND_PRINT(neighborList);

	cout << "\tVertices: ";
	ITERATE_AND_PRINT(verticesList);

	cout << "\tFace areas: ";
	ITERATE_AND_PRINT(areaList);

	cout << "\tComparing fake cells with and without groups:" << endl;
	t = d.targetGroups(red, blue);
	c = d.getCell(x, y, z, t);
	Cell cn = d.getCell(x, y, z);

	neighborList.clear();
	verticesList.clear();
	areaList.clear();

	volume = c.computeVolume();
	c.computeNeighbors(neighborList);
	c.computeVertices(verticesList);
	c.computeFaceAreas(areaList);

	vector<int> neighborListNG;
	vector<double> verticesListNG;
	vector<double> areaListNG;

	double volumeNG = c.computeVolume();
	cn.computeNeighbors(neighborListNG);
	cn.computeVertices(verticesListNG);
	cn.computeFaceAreas(areaListNG);

	cout << "\t\tComparing volumes: ";
	FAIL_IF_DIFFERENT(volume,volumeNG);
	cout << "COMPLETE" << endl;

	cout << "\t\tComparing neighbors: ";
	FAIL_IF_DIFFERENT(neighborList,neighborListNG);
	cout << "COMPLETE" << endl;

	cout << "\t\tComparing vertices: ";
	FAIL_IF_DIFFERENT(verticesList,verticesListNG);
	cout << "COMPLETE" << endl;

	cout << "\t\tComparing face areas: ";
	FAIL_IF_DIFFERENT(areaList,areaListNG);
	cout << "COMPLETE" << endl;

	return true;
}

/*
 * sourcegroup_test
 *
 * Check that source groups work correctly. Since this should just return a vector
 * with specific indices, the test is minimal.
 */
bool sourcegroup_test()
{
	int red, blue, n;
	Diagram d;
	TargetGroup t;
	GROUP_HEADER(n,d,red,blue,t);

	vector<int> redVect;
	vector<int> blueVect;
	vector<int> all;
	for (int i = 0; i < n/2; ++i) {
		redVect.push_back(i);
		all.push_back(i);
	}
	for (int i = n/2; i < n; ++i) {
		blueVect.push_back(i);
		all.push_back(i);
	}

	cout << "Checking correctness of red source group: ";
	FAIL_IF_DIFFERENT(d.sourceGroups(red),redVect);
	cout << "COMPLETE" << endl;

	cout << "Checking correctness of blue source group: ";
	FAIL_IF_DIFFERENT(d.sourceGroups(blue),blueVect);
	cout << "COMPLETE" << endl;

	cout << "Checking correctness of red and blue source group: ";
	FAIL_IF_DIFFERENT(d.sourceGroups(red, blue),all);
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
	TEST_PASS_OR_FAIL("Vertices",vertices_test_targetgroups());
	TEST_PASS_OR_FAIL("Face areas",face_area_test_targetgroups());
	TEST_PASS_OR_FAIL("Volume",volume_test_targetgroups());
	TEST_PASS_OR_FAIL("Fake cell",fake_cell_test_targetgroups());

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
	TEST_PASS_OR_FAIL("Source group",sourcegroup_test());
	
	return 0;
}
