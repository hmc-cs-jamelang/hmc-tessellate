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
	inline Vector3(void);
	inline Vector3(const double x, const double y, const double z);
	inline Vector3 operator + (const Vector3& A) const;
	inline Vector3 operator + (const double A) const;
	inline Vector3 operator * (const double A) const;
	inline Vector3 operator - (const Vector3& A) const;
	inline Vector3 operator / (const double A) const;
	inline bool operator != (const Vector3& A) const;
	inline bool operator == (const Vector3& A) const;
	inline double Dot(const Vector3& A) const;
	inline Vector3 Cross(const Vector3& A) const;
	inline double distanceTo(const Vector3& other) const;
} Vector3;



typedef struct FaceVertex {
	std::vector<struct HalfEdge*> edges;
} FaceVertex;

typedef struct Vertex {
	struct Vector3 position;
	bool deleteFlag;
	bool seen;

	inline Vertex(void);
	inline Vertex(double x, double y, double z);
	inline Vertex(Vector3 pos);
} Vertex;

typedef struct Particle {
	Vector3 position;
	int id;
	int group;
	inline Particle(int Id, double x, double y, double z);
	inline Particle(int Id, Vector3 pos);
	inline Particle(void);

	// index should not be an int, nor should the default group be 0
	inline Particle(double x, double y, double z, int index, int id, int group = 0)
	{
		// yeah no
	}

} Particle;

typedef struct HalfEdge {
	inline HalfEdge( void );
	inline HalfEdge(Vertex* vertex);
	inline HalfEdge(Vertex* vertex, Particle* neighbor);

	// target, next
	inline HalfEdge(Vertex* vertex, HalfEdge* edge2, Particle* neighbor);

	// target, flip, next
	inline HalfEdge(Vertex* vertex, HalfEdge* edge1, HalfEdge* edge2,
			 Particle* neighbor);

	// We should eventually remove this from our struct, or turn it into a
	// class, or find some better way of doing this.
	inline HalfEdge* getPrev();

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
	inline voronoiCell();

	// inline void clear();

	// inline void initialize(Vector3 position, double xmin, double xmax, double ymin, double ymax,
	// 					   double zmin, double zmax);

	inline void initialize(std::string shape, double length, Particle seedParticle, 
					  double maxRadius, double x_min, double x_max,
					  double y_min, double y_max, double z_min, double z_max);

	inline voronoiCell(std::string shape, double length, Particle particle, double maxRadius,
				double x_min, double x_max, double y_min, double y_max,
				double z_min, double z_max);
	inline ~voronoiCell();
	
	struct Particle particle;
	struct Particle neighborParticle;
	double maxRadius;

	inline void cutCell(const Particle& neighbor);

	// CustomInterface IO

	// why was this virtual?
	inline void neighbors(std::vector<int> &v);

	inline void faceAreas(std::vector<double> &v);

	inline double volume();
	inline void vertices(std::vector<double> &v);
	inline void face_vertices(std::vector<int> &v);
	inline void face_areas(std::vector<double> &v);

	inline void drawGnuplot(double dispX, double dispY, double dispZ, FILE* fp, HalfEdge* edge);
	// only public for debugging purposes
	HalfEdge* firstEdge;

	std::vector<FaceVertex*> faceVertices;

private:
	// Finds side of plane that point is on
	inline side planeSide(Vertex* vertex);

	// I love Dani-chan
	// Finds a new edge to be first edge (one that wont be cut off)
	inline HalfEdge* maintainFirstEdge(HalfEdge* edge);

	// Adds a vertex in the middle of an edge
	inline void splitEdge(HalfEdge* edge, Vertex* newVertex);

	// Location of intersection of plane with edge
	inline Vector3 planeEdgeIntersect(HalfEdge* edge);

	// Creates a halfedge and its flip
	inline HalfEdge* addEdgePair(Vertex* vertex1, Vertex* vertex2);

	// Makes a one-sided face on a vertex.
	// Used as an intermediate stage.
	inline HalfEdge* makeSelfLoopOnVertex(Vertex* newVertex);

	// Sets two edges to be each other's flips
	inline void setFlip(HalfEdge* forward, HalfEdge* back);

	// Half of the process of selfLoopOnVertex
	// Makes one halfedge of the face (but not its flip)
	inline HalfEdge* makeOneEdgeFace(Vertex* vertex);

	// orig is pointing outside the plane.
	// This goes until it finds one that crosses the plane.
	inline HalfEdge* findNextIncidentEdge(HalfEdge* orig);

	// Distance from plane
	inline double planeDist(Vertex* vertex);

	// Returns whether or not an incident edge exists,
	// and if so returns that edge in returnEdge
	inline bool findSomeIncidentEdge(HalfEdge* &returnEdge);

	// Resets the seen flag on all edges
	inline void resetEdges(HalfEdge* edge);
	// Resets the seen flag on all edges and vertices
	inline void resetEdgesAndVertices(HalfEdge* edge);

	// Deprecated
	// Volume uses it though...
	inline void reset(HalfEdge* edge);

	// Attempted but failed, supposed to be resetEdgesAndVertices
	// but more robust. Doesn't work.
	inline void seenSearch(std::stack<HalfEdge*>* seenStackEdge,
					std::stack<Vertex*>* seenStackVertex,
					HalfEdge* edge);

	// Delete things, starting with an edge-to-be-deleted
	inline void cleanUp(HalfEdge* edge);
	// Gets the connected component to delete
	inline void deleteSearch(std::stack<HalfEdge*>* deleteStackEdge,
					  std::stack<Vertex*>* deleteStackVertex,
					  HalfEdge* edge);

	// IO helper functions
	inline void getVertex(HalfEdge* testEdge, std::vector<Vertex*> &vertices);
	inline void getFaceVertex(HalfEdge* testEdge);
	inline double tetVolume(Vertex* vertex1, Vertex* vertex2,
					 Vertex* vertex3, Vertex* vertex4);

	inline double triArea(Vertex* vertex1, Vertex* vertex2, Vertex* vertex3);

	// destructor helper
	inline void getEdgeAndVertex(HalfEdge* testEdge, std::stack<HalfEdge*> &edgeStack,
						  std::stack<Vertex*> &vertexStack);

	const double tolerance = 1e-11;
};

// Pulled from http://www.cacr.caltech.edu/~sean/projects/stlib/html/geom/classstlib_1_1geom_1_1CellArrayNeighbors.html#a94f9656a64eb93f12859b17306ba8dd0
// stlib documentation on CellArrayNeighbors
struct Location :
	public std::unary_function<Particle*, std::array<double,3> > {
	result_type
	inline operator()(argument_type r) {
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

	inline voronoiCell* makeCell(Particle particle);
	inline cellContainer(std::vector<Particle> parts, double defaultLen);
	inline cellContainer(std::vector<Particle> parts, double defaultLen,
				  double x_min, double x_max, double y_min, double y_max,
				  double z_min, double z_max);
	inline ~cellContainer();
	inline double sum_cell_volumes();
	inline void put(int id, double x, double y, double z);
	inline double findMaxNeighDist();   

};

#include "vorocode-private.hpp"

#endif
