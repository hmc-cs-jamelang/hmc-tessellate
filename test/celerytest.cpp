/*
 * A set of tests for the Celery class, an implementation of
 * a cell array.
 */

#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>
#include "../src/cellarray.hpp"

using namespace std;
using namespace hmc::spatial;

// Set up random number generation
uniform_real_distribution<double> rand_double(0, 1);
default_random_engine generator;

// Run one test, stop testing on failure.
#define RUN_ONE_TEST(test)					   \
	if (test) {								   \
		cout << "*** TEST PASSED ***" << endl; \
	}										   \
	else {									   \
		cout << "*** TEST FAILED ***" << endl; \
		return 0;							   \
	}

// A simple point type that can be accepted by the cell array
typedef struct TestPoint {
	double x, y, z;

	TestPoint(double a, double b, double c)
		: x(a), y(b), z(c)
	{}

	friend ostream& operator<<(ostream& os, const TestPoint& pt)
	{
		os << "(" << pt.x << ", " << pt.y << ", " << pt.z << ")";
		return os;
	}
} TestPoint;

/*
 * rand_in_range
 *
 * Returns a random double within the specified range.
 */
double rand_in_range(double min, double max)
{
	double normalized = rand_double(generator);
	return min + normalized * (max - min);
}

/*
 * generate_random_points
 *
 * Creates numPoints randomly assigned points within the specified
 * domain.
 */
void generate_random_points(double xmin, double xmax,
							double ymin, double ymax,
							double zmin, double zmax,
							size_t numPoints,
							vector<TestPoint>& points)
{
	for (size_t i = 0; i < numPoints; ++i) {
		double x, y, z;

		x = rand_in_range(xmin, xmax);
		y = rand_in_range(ymin, ymax);
		z = rand_in_range(zmin, zmax);

		points.push_back(TestPoint(x,y,z));
	}
}

/*
 * print_container
 *
 * Print the contents of a container on a single line.
 */
template<typename T>
void print_container(T begin, T end)
{
	for (T i = begin; i != end; ++i) {
		cout << *i << ", ";
	}
}

/*
 * check_cells_delimiters
 *
 * Check that the vector of cell information matches the vector of
 * delimiter information.
 */
bool check_cells_delimiters(const vector<size_t> & cells,
							const vector<size_t> & delimiters)
{
	size_t delimiters_size = delimiters.size();
	size_t cells_size = cells.size();

	// The first cell is always cell zero, and always starts at the
	// first index
	if (delimiters[0] != 0) {
		cout << "DELIMITER ERROR: first delimiter is not 0" << endl;
		return false;
	}

	for (size_t i = 1; i < delimiters_size; ++i) {
		// Ensure that the delimiters do not set the beginning of a cell
		// after some points in that cell
		if (i <= cells[delimiters[i]-1] && delimiters[i] != delimiters[i-1]) {
			cout << "DELIMITER ERROR: delimiter " << i <<
				" begins after cell begins" << endl;
			return false;
		}

		// Once there are no more points, the rest of the groups should be
		// "one-off-the-end"
		if (delimiters[i] == cells_size) {
			for (size_t j = i + 1; j < delimiters_size; ++j) {
				if (delimiters[j] != cells_size) {
					cout << "DELIMITER ERROR: delimiter " << j <<
						" should be " << cells_size << endl;
					return false;
				}
			}

			return true;
		}

		// If there are more points, check that beginning point is either
		// in the correct cell, or the beginning of the cell is the same as
		// the end. This case should NEVER occur for the last entry in the
		// delimiters vector, since it must be "one-off-the-end" and
		// therefore caught by the previous case, so we don't need to worry
		// about indexing at i+1
		else if (i != cells[delimiters[i]] && delimiters[i] != delimiters[i+1]) {
			cout << "DELIMITER ERROR: delimiter " << i <<
				" gives an index to cell " << delimiters[i] << endl;
			return false;
		}
	}

	// The final delimiter should ALWAYS be "one-off-the-end", so this test
	// should fail we do not find such a delimiter
	cout << "DELIMITER ERROR: the last delimiter is not one-off-the-end" << endl;
	return false;
}

/*
 * insert_n_test
 *
 * Test insertion of n points with the given bounding box.
 */
bool insert_n_test(double xmin, double xmax, double ymin, double ymax,
				   double zmin, double zmax, size_t n)
{
	vector<TestPoint> pts;
	generate_random_points(xmin, xmax, ymin, ymax, zmin, zmax, n, pts);

	Celery<TestPoint> sds(xmin, xmax, ymin, ymax, zmin, zmax, pts.begin(), pts.end());

	vector<TestPoint> cellarray = sds.getPoints();
	if (pts.size() != cellarray.size()) {
		return false;
	}

	vector<size_t> cells = sds.getCells();
	vector<size_t> delimiters = sds.getDelimiters();

	cout << "Number of cells per dimension: " << sds.getNumCellsDim() << endl;
	cout << "Number of cells: " << (sds.getDelimiters().size() - 1) << endl;
	cout << "x inverse cell size: " << sds.getCellSizeInvX() << endl;
	cout << "y inverse cell size: " << sds.getCellSizeInvY() << endl;
	cout << "z inverse cell size: " << sds.getCellSizeInvZ() << endl;

	// cout << "Points: ";
	// print_container(cellarray.begin(), cellarray.end());
	// cout << endl;

	// cout << "Cells: ";
	// print_container(cells.begin(), cells.end());
	// cout << endl;

	// cout << "Delimiters: ";
	// print_container(delimiters.begin(), delimiters.end());
	// cout << endl;

	if (!check_cells_delimiters(cells, delimiters)) {
		cout << "TEST FAILED: cells and delimiters do not match" << endl;
		return false;
	}

	return true;
}

/*
 * empty_cells_test
 *
 * Test that particles are inserted correctly when some cells are missing.
 * Must have dmin <= pmin < pmax <= dmax
 */
bool empty_cells_test(double dmin, double dmax, double pmin, double pmax, std::size_t n)
{
	vector<TestPoint> pts;
	generate_random_points(pmin, pmax, pmin, pmax, pmin, pmax, n, pts);

	Celery<TestPoint> sds(dmin, dmax, dmin, dmax, dmin, dmax, pts.begin(), pts.end());

	vector<TestPoint> cellarray = sds.getPoints();
	if (pts.size() != cellarray.size()) {
		return false;
	}

	vector<size_t> cells = sds.getCells();
	vector<size_t> delimiters = sds.getDelimiters();

	cout << "Number of cells per dimension: " << sds.getNumCellsDim() << endl;
	cout << "x inverse cell size: " << sds.getCellSizeInvX() << endl;
	cout << "y inverse cell size: " << sds.getCellSizeInvY() << endl;
	cout << "z inverse cell size: " << sds.getCellSizeInvZ() << endl;

	// cout << "Points: ";
	// print_container(cellarray.begin(), cellarray.end());
	// cout << endl;

	// cout << "Cells: ";
	// print_container(cells.begin(), cells.end());
	// cout << endl;

	// cout << "Delimiters: ";
	// print_container(delimiters.begin(), delimiters.end());
	// cout << endl;

	if (!check_cells_delimiters(cells, delimiters)) {
		cout << "TEST FAILED: cells and delimiters do not match" << endl;
		return false;
	}

	return true;
}

/*
 * initialize_test
 *
 * Test that the initialize method of creating a cell array works correctly.
 */
bool initialize_test(double xmin, double xmax, double ymin, double ymax,
				   double zmin, double zmax, size_t n)
{
	Celery<TestPoint> sds;
	std::vector<TestPoint> pts;

	generate_random_points(xmin, xmax, ymin, ymax, zmin, zmax, n, pts);
	sds.initialize(pts.begin(), pts.end());

	vector<TestPoint> cellarray = sds.getPoints();
	if (pts.size() != cellarray.size()) {
		return false;
	}

	vector<size_t> cells = sds.getCells();
	vector<size_t> delimiters = sds.getDelimiters();

	cout << "Number of cells per dimension: " << sds.getNumCellsDim() << endl;
	cout << "x inverse cell size: " << sds.getCellSizeInvX() << endl;
	cout << "y inverse cell size: " << sds.getCellSizeInvY() << endl;
	cout << "z inverse cell size: " << sds.getCellSizeInvZ() << endl;
	cout << "Lower x boundary: " << sds.getXMin() << endl;
	cout << "Upper x boundary: " << sds.getXMax() << endl;
	cout << "Lower y boundary: " << sds.getYMin() << endl;
	cout << "Upper y boundary: " << sds.getYMax() << endl;
	cout << "Lower z boundary: " << sds.getZMin() << endl;
	cout << "Upper z boundary: " << sds.getZMax() << endl;

	// cout << "Points: ";
	// print_container(cellarray.begin(), cellarray.end());
	// cout << endl;

	// cout << "Cells: ";
	// print_container(cells.begin(), cells.end());
	// cout << endl;

	// cout << "Delimiters: ";
	// print_container(delimiters.begin(), delimiters.end());
	// cout << endl;

	if (!check_cells_delimiters(cells, delimiters)) {
		cout << "TEST FAILED: cells and delimiters do not match" << endl;
		return false;
	}

	return true;
}

int main()
{
	cout << "*** TESTING INSERTION OF A SINGLE POINT ***" << endl;
	RUN_ONE_TEST(insert_n_test(0, 1, 0, 1, 0, 1, 1));

	cout << "*** TESTING INSERTION OF FIVE POINTS ***" << endl;
	RUN_ONE_TEST(insert_n_test(0, 1, 0, 1, 0, 1, 5));

	cout << "*** TESTING INSERTION OF A HUNDRED POINTS ***" << endl;
	RUN_ONE_TEST(insert_n_test(0, 1, 0, 1, 0, 1, 100));

	cout << "*** TESTING INSERTION OF A HUNDRED POINTS WITH ALTERNATIVE BOUNDING BOX ***" << endl;
	RUN_ONE_TEST(insert_n_test(50, 51, 50, 51, 50, 51, 100));

	cout << "*** TESTING INSERTION OF A HUNDRED POINTS WITH NEGATIVE BOUNDING BOX ***" << endl;
	RUN_ONE_TEST(insert_n_test(-10, -9, -10, -9, -10, -9, 100));

	// Million
	cout << "*** TESTING INSERTION OF A MILLION POINTS ***" << endl;
	RUN_ONE_TEST(insert_n_test(0, 1, 0, 1, 0, 1, 1000000));

	cout << "*** TESTING INSERTION OF A MILLION POINTS WITH ALTERNATIVE BOUNDING BOX ***" << endl;
	RUN_ONE_TEST(insert_n_test(100, 101, 100, 101, 100, 101, 1000000));

	cout << "*** TESTING INSERTION OF A MILLION POINTS WITH NEGATIVE BOUNDING BOX ***" << endl;
	RUN_ONE_TEST(insert_n_test(-100, -99, -100, -99, -100, -99, 1000000));

	// 10 million
	// cout << "*** TESTING INSERTION OF TEN MILLION POINTS ***" << endl;
	// RUN_ONE_TEST(insert_n_test(0, 1, 0, 1, 0, 1, 10000000));

	// cout << "*** TESTING INSERTION OF TEN MILLION POINTS WITH ALTERNATIVE BOUNDING BOX ***" << endl;
	// RUN_ONE_TEST(insert_n_test(100, 101, 100, 101, 100, 101, 10000000));

	// cout << "*** TESTING INSERTION OF TEN MILLION POINTS WITH NEGATIVE BOUNDING BOX ***" << endl;
	// RUN_ONE_TEST(insert_n_test(-100, -99, -100, -99, -100, -99, 10000000));

	// Test partially filled cells
	cout << "*** TESTING INSERTION INTO FIRST HALF OF BOUNDING BOX ***" << endl;
	RUN_ONE_TEST(empty_cells_test(0, 2, 0, 1, 1000000));

	cout << "*** TESTING INSERTION INTO SECOND HALF OF BOUNDING BOX ***" << endl;
	RUN_ONE_TEST(empty_cells_test(0, 2, 1, 2, 1000000));

	// Test the initialize method of creating a cell array
	cout << "*** TESTING THE INITIALIZATION METHOD OF INSERTING POINTS ***" << endl;
	RUN_ONE_TEST(initialize_test(0, 1, 0, 1, 0, 1, 1000000));

	return 0;
}
