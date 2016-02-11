// Voro++ replacement header
//
// Sandia 2015-16 Clinic Team
// sandia2015l@cs.hmc.edu

#ifndef VORO_CELL
#define VORO_CELL

#include <stack>
#include <vector>
#include <iostream>
#include <string>
#include <math.h>
#include <stlib/geom/orq/KDTree.h>
#include <stlib/geom/orq/CellArrayNeighbors.h>

// Taken from http://www.flipcode.com/archives/Faster_Vector_Math_Using_Templates.shtml
typedef struct Vector3
{
	double X, Y, Z;
	Vector3(void){};//hack{}
	Vector3(const double x, const double y, const double z);
	Vector3 operator + (const Vector3& A) const;
	Vector3 operator + (const double A) const;
	Vector3 operator * (const double A) const;
	Vector3 operator - (const Vector3& A) const;
	Vector3 operator / (const double A) const;
	bool operator != (const Vector3& A) const;
	bool operator == (const Vector3& A) const;
	double Dot(const Vector3& A) const;
	Vector3 Cross(const Vector3& A) const;
	double distanceTo(const Vector3& other) const;
} Vector3;



typedef struct FaceVertex {
	std::vector<struct HalfEdge*> edges;
} FaceVertex;

typedef struct Vertex {
	struct Vector3 position;
	bool deleteFlag;
	bool seen;

	Vertex(void);
	Vertex(double x, double y, double z);
	Vertex(Vector3 pos);
} Vertex;

typedef struct Particle {
	Vector3 position;
	int id;
	int group;
	Particle(int Id, double x, double y, double z);
	Particle(int Id, Vector3 pos);
	Particle(void){};//hack{}

	// index should not be an int, nor should the default group be 0
	Particle(double x, double y, double z, int index, int id, int group = 0)
	{
		// yeah no
	}

} Particle;

typedef struct HalfEdge {
	HalfEdge( void );
	HalfEdge(Vertex* vertex);
	HalfEdge(Vertex* vertex, Particle* neighbor);

	// target, next
	HalfEdge(Vertex* vertex, HalfEdge* edge2, Particle* neighbor);

	// target, flip, next
	HalfEdge(Vertex* vertex, HalfEdge* edge1, HalfEdge* edge2,
			 Particle* neighbor);

	// We should eventually remove this from our struct, or turn it into a
	// class, or find some better way of doing this.
	HalfEdge* getPrev();

	Vertex* target;
	struct HalfEdge* flip;
	struct HalfEdge* next;
	bool deleteFlag;
	bool seen;

	FaceVertex* face;

	// The neighbor particle whose cutting plane created the face
	// associated with this HalfEdge
	const Particle* creator;
} HalfEdge;



enum side {inside, outside, incident};

class voronoiCell{
public:
	voronoiCell()
	{
		// do nothing
	}

	void clear()
	{
		// Probably should reset some stuff
	}

	void initialize(Vector3 position, double xmin, double xmax, double ymin, double ymax,
					double zmin, double zmax)
	{
		// Should, you know, initialize stuff
	}

	voronoiCell(std::string shape, double length, Particle particle, double maxRadius,
				double x_min, double x_max, double y_min, double y_max,
				double z_min, double z_max);
	~voronoiCell(){};//hack{}
	
	struct Particle particle;
	struct Particle neighborParticle;
	double maxRadius;

	void cutCell(const Particle& neighbor){};//hack{}

	// CustomInterface IO

	// why was this virtual?
	void neighbors(std::vector<int> &v){};//hack{}

	void faceAreas(std::vector<double> &v);

	double volume(){};//hack{}
	void vertices(std::vector<double> &v){};//hack{}
	void face_vertices(std::vector<int> &v);
	void face_areas(std::vector<double> &v){};//hack{}

	void drawGnuplot(double dispX, double dispY, double dispZ, FILE* fp, HalfEdge* edge);
	// only public for debugging purposes
	HalfEdge* firstEdge;

	std::vector<FaceVertex*> faceVertices;

private:
	// Finds side of plane that point is on
	side planeSide(Vertex* vertex);

	// I love Dani-chan
	// Finds a new edge to be first edge (one that wont be cut off)
	HalfEdge* maintainFirstEdge(HalfEdge* edge);

	// Adds a vertex in the middle of an edge
	void splitEdge(HalfEdge* edge, Vertex* newVertex);

	// Location of intersection of plane with edge
	Vector3 planeEdgeIntersect(HalfEdge* edge);

	// Creates a halfedge and its flip
	HalfEdge* addEdgePair(Vertex* vertex1, Vertex* vertex2);

	// Makes a one-sided face on a vertex.
	// Used as an intermediate stage.
	HalfEdge* makeSelfLoopOnVertex(Vertex* newVertex);

	// Sets two edges to be each other's flips
	void setFlip(HalfEdge* forward, HalfEdge* back);

	// Half of the process of selfLoopOnVertex
	// Makes one halfedge of the face (but not its flip)
	HalfEdge* makeOneEdgeFace(Vertex* vertex);

	// orig is pointing outside the plane.
	// This goes until it finds one that crosses the plane.
	HalfEdge* findNextIncidentEdge(HalfEdge* orig);

	// Distance from plane
	double planeDist(Vertex* vertex);

	// Returns whether or not an incident edge exists,
	// and if so returns that edge in returnEdge
	bool findSomeIncidentEdge(HalfEdge* &returnEdge);

	// Resets the seen flag on all edges
	void resetEdges(HalfEdge* edge);
	// Resets the seen flag on all edges and vertices
	void resetEdgesAndVertices(HalfEdge* edge);

	// Deprecated
	// Volume uses it though...
	void reset(HalfEdge* edge);

	// Attempted but failed, supposed to be resetEdgesAndVertices
	// but more robust. Doesn't work.
	void seenSearch(std::stack<HalfEdge*>* seenStackEdge,
					std::stack<Vertex*>* seenStackVertex,
					HalfEdge* edge);

	// Delete things, starting with an edge-to-be-deleted
	void cleanUp(HalfEdge* edge);
	// Gets the connected component to delete
	void deleteSearch(std::stack<HalfEdge*>* deleteStackEdge,
					  std::stack<Vertex*>* deleteStackVertex,
					  HalfEdge* edge);

	// IO helper functions
	void getVertex(HalfEdge* testEdge, std::vector<Vertex*> &vertices);
	void getFaceVertex(HalfEdge* testEdge);
	double tetVolume(Vertex* vertex1, Vertex* vertex2,
					 Vertex* vertex3, Vertex* vertex4);

	double triArea(Vertex* vertex1, Vertex* vertex2, Vertex* vertex3);

	// destructor helper
	void getEdgeAndVertex(HalfEdge* testEdge, std::stack<HalfEdge*> &edgeStack,
						  std::stack<Vertex*> &vertexStack);

	const double tolerance = 1e-11;
};

// Pulled from http://www.cacr.caltech.edu/~sean/projects/stlib/html/geom/classstlib_1_1geom_1_1CellArrayNeighbors.html#a94f9656a64eb93f12859b17306ba8dd0
// stlib documentation on CellArrayNeighbors
struct Location :
	public std::unary_function<Particle*, std::array<double,3> > {
	result_type
	operator()(argument_type r) {
		result_type location = {{r->position.X, r->position.Y, r->position.Z}};
		return location;
	}
};

class cellContainer{
public:
	stlib::geom::CellArrayNeighbors<double, 3,
									Particle*,
									Location > sds;
	struct std::vector<Particle> particles;
	struct std::vector<voronoiCell*> cells;
	double defaultLength;
	bool calculated = false;
	double x_min, x_max, y_min, y_max, z_min, z_max;

	voronoiCell* makeCell(Particle particle);
	cellContainer(std::vector<Particle> parts, double defaultLen);
	cellContainer(std::vector<Particle> parts, double defaultLen,
				  double x_min, double x_max, double y_min, double y_max,
				  double z_min, double z_max);
	~cellContainer();
	double sum_cell_volumes();
	void put(int id, double x, double y, double z);
	double findMaxNeighDist();   

};

//#include "vorocode-private.hpp"

#endif
