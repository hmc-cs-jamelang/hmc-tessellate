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
#include "structpool.hpp"

struct HalfEdge;
struct Vertex;

typedef StructPool<HalfEdge>::Index EdgeIndex;
typedef StructPool<Vertex>::Index VertexIndex;
// typedef StructPool<voronoiCell>::Index CellIndex;
// typedef StructPool<FaceVertex>::Index FaceVertexIndex;

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
	inline Vector3 normalized() const;
	inline double Dot(const Vector3& A) const;
	inline Vector3 Cross(const Vector3& A) const;
	inline double distanceTo(const Vector3& other) const;
} Vector3;



typedef struct FaceVertex {
	// std::vector<EdgeIndex> edges;
	EdgeIndex firstEdge;
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

	// Should be removed later hopefully
	size_t index;

	inline Particle(int Id, double x, double y, double z, size_t index);
	inline Particle(int Id, Vector3 pos, size_t index);
	inline Particle(void);

	// index should not be an int, nor should the default group be 0
	inline Particle(double x, double y, double z, int index, int id, int group = 0)
	{
		// yeah no
	}

	inline x()
	{
		position.X;
	}

	inline y()
	{
		position.Y;
	}

	inline z()
	{
		position.Z;
	}

} Particle;

typedef struct HalfEdge {
	inline HalfEdge( void );
	inline HalfEdge(VertexIndex vertex);
	inline HalfEdge(VertexIndex vertex, size_t neighbor);

	// target, next
	inline HalfEdge(VertexIndex vertex, EdgeIndex edge2, size_t neighbor);

	// target, flip, next
	inline HalfEdge(VertexIndex vertex, EdgeIndex edge1, EdgeIndex edge2,
			 size_t neighbor);

	// We should eventually remove this from our struct, or turn it into a
	// class, or find some better way of doing this.
	inline EdgeIndex getPrev();

	VertexIndex target;
	EdgeIndex flip;
	EdgeIndex next;
	bool deleteFlag;
	bool seen;

	// FaceVertex* face;

	// The neighbor particle whose cutting plane created the face
	// associated with this HalfEdge
	size_t creator;
} HalfEdge;



enum side {inside, outside, incident};

class voronoiCell{
public:
	inline voronoiCell()
	{
		// do nothing
	}

	inline void clear()
	{
		// Probably should reset some stuff
	}

	// inline voronoiCell& operator=(voronoiCell& rhs) = default;

	inline void initialize(Vector3 position, double xmin, double xmax, double ymin, double ymax,
					double zmin, double zmax)
	{
		// Should, you know, initialize stuff
	}

	inline void reconstruct(const size_t& particleIndex, const Particle& particle, const double& maxRadius,
				const double& x_min, const double& x_max, const double& y_min, const double& y_max,
				const double& z_min, const double& z_max);
	// inline ~voronoiCell();

	struct Particle particle;
	struct Particle neighborParticle;
	double maxRadius;

	inline void cutCell(const Particle& neighbor, size_t index);

	// CustomInterface IO

	// why was this virtual?
	inline void neighbors(std::vector<int> &v);

	inline void faceAreas(std::vector<double> &v);

	inline double volume();
	inline void computeVertices(std::vector<double> &v);
	inline void face_vertices(std::vector<int> &v);
	inline void face_areas(std::vector<double> &v);

	inline void drawGnuplot(double dispX, double dispY, double dispZ, FILE* fp, EdgeIndex edge);
	// only public for debugging purposes
	EdgeIndex firstEdge;

	StructPool<HalfEdge> edges;
	StructPool<Vertex> vertices;
	std::vector<FaceVertex> faceVertices;
	std::vector<Particle*> NeighborParticles;

	inline std::size_t get_memory_usage();

private:

	std::vector<EdgeIndex> deleteStackEdge;

    std::vector<VertexIndex> deleteStackVertex;

    std::vector<EdgeIndex> seenStackEdge;

    std::vector<VertexIndex> seenStackVertex;

    std::vector<EdgeIndex> outsideEdges;

    std::vector<EdgeIndex> testEdges;



	// Finds side of plane that point is on
	inline side planeSide(VertexIndex vertex);

	// I love Dani-chan
	// Finds a new edge to be first edge (one that wont be cut off)
	inline EdgeIndex maintainFirstEdge(EdgeIndex edge);

	// Adds a vertex in the middle of an edge
	inline void splitEdge(EdgeIndex edge, VertexIndex newVertex);

	// Location of intersection of plane with edge
	inline Vector3 planeEdgeIntersect(EdgeIndex edge);

	// Creates a halfedge and its flip
	inline EdgeIndex addEdgePair(VertexIndex vertex1, VertexIndex vertex2);

	// Makes a one-sided face on a vertex.
	// Used as an intermediate stage.
	inline EdgeIndex makeSelfLoopOnVertex(VertexIndex newVertex);

	// Sets two edges to be each other's flips
	inline void setFlip(EdgeIndex forward, EdgeIndex back);

	// Half of the process of selfLoopOnVertex
	// Makes one halfedge of the face (but not its flip)
	inline EdgeIndex makeOneEdgeFace(VertexIndex vertex);

	// orig is pointing outside the plane.
	// This goes until it finds one that crosses the plane.
	inline EdgeIndex findNextIncidentEdge(EdgeIndex orig);

	// Distance from plane
	inline double planeDist(VertexIndex vertex);

	// Returns whether or not an incident edge exists,
	// and if so returns that edge in returnEdge
	inline bool findSomeIncidentEdge(EdgeIndex &returnEdge);

	// Resets the seen flag on all edges
	inline void resetEdges(EdgeIndex edge);
	// Resets the seen flag on all edges and vertices
	inline void resetEdgesAndVertices(EdgeIndex edge);

	// Deprecated
	// Volume uses it though...
	inline void reset(EdgeIndex edge);

	// Attempted but failed, supposed to be resetEdgesAndVertices
	// but more robust. Doesn't work. Should probably get rid of
	// due to memory pools
	inline void seenSearch(EdgeIndex edge);

	// Delete things, starting with an edge-to-be-deleted
	inline void cleanUp(EdgeIndex edge);

	// Gets the connected component to delete.
	// Can probably do this easier with memory pools
	inline void deleteSearch(EdgeIndex edge);

	// IO helper functions
	inline void getVertex(EdgeIndex testEdge, std::vector<VertexIndex> &vertices);
	inline void getFaceVertex(EdgeIndex testEdge);
	inline double tetVolume(VertexIndex vertex1, VertexIndex vertex2,
					 VertexIndex vertex3, VertexIndex vertex4);

	inline double triArea(VertexIndex vertex1, VertexIndex vertex2, VertexIndex vertex3);

	// destructor helper
	inline void getEdgeAndVertex(EdgeIndex testEdge, std::stack<EdgeIndex> &edgeStack,
						  std::stack<VertexIndex> &vertexStack);

	static constexpr double tolerance = 1e-11;
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
	double defaultLength;
	bool calculated = false;
	double x_min, x_max, y_min, y_max, z_min, z_max;
	// struct std::vector<voronoiCell*> cells;
	std::vector<Particle*> NeighborParticles;


	// inline void makeCell(size_t particleIndex, voronoiCell& cell);
	inline void makeCell(size_t particleIndex, voronoiCell& cell, double searchDist) const;
	inline cellContainer(std::vector<Particle> parts, double defaultLen);
	inline cellContainer(std::vector<Particle> parts, double defaultLen,
				  double x_min, double x_max, double y_min, double y_max,
				  double z_min, double z_max);
	inline ~cellContainer();
	inline void initialize();
	inline double sum_cell_volumes();
	inline void put(int id, double x, double y, double z);
	inline std::vector<double> findMaxNeighDist(double scale);
	// inline std::size_t get_memory_usage();

};

#include "vorocode-private.hpp"

#endif
