// added just for one use of abs; let me know if there is a better way
// to do this
#pragma once
#ifndef voro_priv
#define voro_priv
#include <cmath>
#include <cstdlib>
#include <stdio.h>
#include <limits>



// vector3 methods

Vector3::Vector3(void) {}
Vector3::Vector3(const double x, const double y, const double z)
{ X = x; Y = y; Z = z; }

Vector3 Vector3::operator + (const Vector3& A) const
{ return Vector3( X + A.X, Y + A.Y, Z + A.Z ); }

Vector3 Vector3::operator + (const double A) const
{ return Vector3( X + A, Y + A, Z + A ); }

Vector3 Vector3::operator * (const double A) const
{ return Vector3( X * A, Y * A, Z * A ); }

Vector3 Vector3::operator - (const Vector3& A) const
{ return Vector3( X - A.X, Y - A.Y, Z - A.Z ); }

Vector3 Vector3::operator / (const double A) const
{ return Vector3( X / A, Y / A, Z / A ); }

bool Vector3::operator != (const Vector3& A) const
{ return (X != A.X) || (Y != A.Y) || (Z != A.Z); }

bool Vector3::operator == (const Vector3& A) const
{ return (X == A.X) && (Y == A.Y) && (Z == A.Z); }

double Vector3::Dot(const Vector3& A) const
{ return A.X*X + A.Y*Y + A.Z*Z; }

Vector3 Vector3::Cross(const Vector3& A) const
{ return Vector3( Y * A.Z - Z * A.Y, Z * A.X - X * A.Z, X * A.Y - Y * A.X ); }

double Vector3::distanceTo(const Vector3& other) const
{ return sqrt(pow(X-other.X, 2) + pow(Y-other.Y, 2) + pow(Z-other.Z, 2)); }

// vertex constructors
Vertex::Vertex(void) {
    deleteFlag = false;
    seen = false;
}

Vertex::Vertex(double x, double y, double z) {
    position = Vector3(x, y, z);
    deleteFlag = false;
    seen = false;
}

Vertex::Vertex(Vector3 pos) {
    position = pos;
    deleteFlag = false;
    seen = false;
}

Particle::Particle(int Id, double x, double y, double z) {
    position = Vector3(x, y, z);
    id = Id;
}

Particle::Particle(int Id, Vector3 pos) {
    position = pos;
    id = Id;
}

Particle::Particle( void ) {
    id = -1;
}


// halfedge constructors

HalfEdge::HalfEdge( void ) {
    deleteFlag = false;
    seen = false;
}

HalfEdge::HalfEdge(VertexIndex vertex) {
    target = vertex;
    deleteFlag = false;
    seen = false;
}

HalfEdge::HalfEdge(VertexIndex vertex, Particle* neighbor) {
    target = vertex;
    deleteFlag = false;
    seen = false;
    creator = neighbor;
}

// This constructor builds a HalfEdge with a specified target and next.
HalfEdge::HalfEdge(VertexIndex vertex, EdgeIndex edge2, Particle* neighbor) {
    target = vertex;
    next = edge2;
    deleteFlag = false;
    seen = false;
    creator = neighbor;
}

HalfEdge::HalfEdge(VertexIndex vertex, EdgeIndex edge1, EdgeIndex edge2, Particle* neighbor) {
    target = vertex;
    flip = edge1;
    next = edge2;
    deleteFlag = false;
    seen = false;
    creator = neighbor;
}

// // Intended to get the previous edge, should be used rarely.
// EdgeIndex HalfEdge::getPrev() {
//     EdgeIndex other = this;

//     //TODO: Optimize
//     // HalfEdge* other = next->next;
//     while (other->next != this) {
//         other = other->next;
//     }
//     return other;
// }

// HalfEdge* firstEdge = getFirstEdge();

voronoiCell::voronoiCell(std::string shape, double length, Particle seedParticle,
                         double maxRadius, double x_min, double x_max,
                         double y_min, double y_max, double z_min, double z_max) {
    particle = seedParticle;
    maxRadius = maxRadius;
    faceVertices = std::vector<FaceVertex>();

    if (shape == "cube") {
        // Create half-edge arrays corresponding to the faces of the cube.
        EdgeIndex plusZ[4], minusZ[4], plusY[4], minusY[4], plusX[4], minusX[4];

        // Initialize the vertices of the cube manually.
        //
        // Vertices are labeled by octant (with respect to the central
        // particle) as follows:
        // First octant: positive x,y,z
        // Second octant: positive y,z, negative x
        // Third octant: positive z, negative x,y
        // Fourth octant: positive x,z, negative y
        // Fifth octant: positive x,y, negative z
        // Sixth octant: positive y, negative x,z
        // Seventh octant: negative x,y,z
        // Eighth octant: positive x, negative y,z

        VertexIndex oct1 = vertices.create(Vertex(x_max,
                                  y_max,
                                  z_max));

        VertexIndex oct2 = vertices.create(Vertex(x_min,
                                  y_max,
                                  z_max));

        VertexIndex oct3 = vertices.create(Vertex(x_min,
                                  y_min,
                                  z_max));

        VertexIndex oct4 = vertices.create(Vertex(x_max,
                                  y_min,
                                  z_max));

        VertexIndex oct5 = vertices.create(Vertex(x_max,
                                  y_max,
                                  z_min));

        VertexIndex oct6 = vertices.create(Vertex(x_min,
                                  y_max,
                                  z_min));

        VertexIndex oct7 = vertices.create(Vertex(x_min,
                                  y_min,
                                  z_min));

        VertexIndex oct8 = vertices.create(Vertex(x_max,
                                  y_min,
                                  z_min));

        // We will construct each face by beginning with the edge whose target
        // is in the lowest-numbered quadrant, then continue around the face
        // by iteratively creating the prev of the most recently created edge.
        //
        // Traversing the nexts around a face leads if you are looking
        // directly at it should lead to clockwise traversal of the face.

        // Note also our choice of default initialization for the particle
        // associated with each HalfEdge; since no planes have cut this cell
        // yet, we associate it with the particle at the center of this cell.
        plusZ[0] = edges.create(HalfEdge(oct1, &particle));             // O2 to O1

        plusZ[1] = edges.create(HalfEdge(oct2, plusZ[0], &particle));   // O3 to O2

        plusZ[2] = edges.create(HalfEdge(oct3, plusZ[1], &particle));   // O4 to O3

        plusZ[3] = edges.create(HalfEdge(oct4, plusZ[2], &particle));   // O1 to O4

        edges[plusZ[0]].next = plusZ[3];

        minusZ[0] = edges.create(HalfEdge(oct5, &particle));            // O8 to O5

        minusZ[1] = edges.create(HalfEdge(oct8, minusZ[0], &particle)); // O7 to O8

        minusZ[2] = edges.create(HalfEdge(oct7, minusZ[1], &particle)); // O6 to O7

        minusZ[3] = edges.create(HalfEdge(oct6, minusZ[2], &particle)); // O5 to O6

        edges[minusZ[0]].next = minusZ[3];

        plusY[0] = edges.create(HalfEdge(oct1, &particle));             // O5 to O1

        plusY[1] = edges.create(HalfEdge(oct5, plusY[0], &particle));   // O6 to O5

        plusY[2] = edges.create(HalfEdge(oct6, plusY[1], &particle));   // O2 to O6

        plusY[3] = edges.create(HalfEdge(oct2, plusY[2], &particle));   // O1 to O2

        edges[plusY[0]].next = plusY[3];

        minusY[0] = edges.create(HalfEdge(oct3, &particle));            // O7 to O3

        minusY[1] = edges.create(HalfEdge(oct7, minusY[0], &particle)); // O8 to O7

        minusY[2] = edges.create(HalfEdge(oct8, minusY[1], &particle)); // O4 to O8

        minusY[3] = edges.create(HalfEdge(oct4, minusY[2], &particle)); // O3 to O4

        edges[minusY[0]].next = minusY[3];

        plusX[0] = edges.create(HalfEdge(oct1, &particle));             // O4 to O1

        plusX[1] = edges.create(HalfEdge(oct4, plusX[0], &particle));   // O8 to O4

        plusX[2] = edges.create(HalfEdge(oct8, plusX[1], &particle));   // O5 to O8

        plusX[3] = edges.create(HalfEdge(oct5, plusX[2], &particle));   // O1 to O5

        edges[plusX[0]].next = plusX[3];

        minusX[0] = edges.create(HalfEdge(oct2, &particle));            // O6 to O2

        minusX[1] = edges.create(HalfEdge(oct6, minusX[0], &particle)); // O7 to O6

        minusX[2] = edges.create(HalfEdge(oct7, minusX[1], &particle)); // O3 to O7

        minusX[3] = edges.create(HalfEdge(oct3, minusX[2], &particle)); // O2 to O3

        edges[minusX[0]].next = minusX[3];

        // Now we have to set the flips manually
        edges[plusZ[0]].flip = plusY[3];
        edges[plusZ[1]].flip = minusX[3];
        edges[plusZ[2]].flip = minusY[3];
        edges[plusZ[3]].flip = plusX[0];

        edges[minusZ[0]].flip = plusX[2];
        edges[minusZ[1]].flip = minusY[1];
        edges[minusZ[2]].flip = minusX[1];
        edges[minusZ[3]].flip = plusY[1];

        edges[plusY[0]].flip = plusX[3];
        edges[plusY[1]].flip = minusZ[3];
        edges[plusY[2]].flip = minusX[0];
        edges[plusY[3]].flip = plusZ[0];

        edges[minusY[0]].flip = minusX[2];
        edges[minusY[1]].flip = minusZ[1];
        edges[minusY[2]].flip = plusX[1];
        edges[minusY[3]].flip = plusZ[2];

        edges[plusX[0]].flip = plusZ[3];
        edges[plusX[1]].flip = minusY[2];
        edges[plusX[2]].flip = minusZ[0];
        edges[plusX[3]].flip = plusY[0];

        edges[minusX[0]].flip = plusY[2];
        edges[minusX[1]].flip = minusZ[2];
        edges[minusX[2]].flip = minusY[0];
        edges[minusX[3]].flip = plusZ[1];

        // Set firstEdge to our favorite edge.

        firstEdge = plusX[0];

    }
}

// voronoiCell& voronoiCell::operator=(voronoiCell rhs) {
//     edges = rhs.edges;
//     vertices = rhs.vertices;
//     particle = rhs.particle;
//     firstEdge = rhs.firstEdge;
//     faceVertices = rhs.faceVertices;
//     maxRadius = rhs.maxRadius;

// }

// voronoiCell::~voronoiCell() {

//     // resetEdgesAndVertices(firstEdge);
//     // std::stack<EdgeIndex> edgeStack;
//     // std::stack<VertexIndex> vertexStack;

//     // getEdgeAndVertex(firstEdge, edgeStack, vertexStack);
//     // resetEdgesAndVertices(firstEdge);


//     // while (!edgeStack.empty()) {
//     //     edges.destroy(edgeStack.top());

//     //     edgeStack.pop();
//     // }
//     // while (!vertexStack.empty()) {
//     //     vertices.destroy(vertexStack.top());

//     //     vertexStack.pop();
//     // }

//     while (!faceVertices.empty()) {
//         delete faceVertices.back();

//         faceVertices.pop_back();
//     }
// }

void voronoiCell::getEdgeAndVertex(EdgeIndex testEdge, std::stack<EdgeIndex> &edgeStack,
                                   std::stack<VertexIndex> &vertexStack) {
    if (!edges[testEdge].seen) {

        edges[testEdge].seen = true;

        edgeStack.push(testEdge);

        if (!vertices[edges[testEdge].target].seen) {
            vertices[edges[testEdge].target].seen = true;
            vertexStack.push(edges[testEdge].target);
        }



        getEdgeAndVertex(edges[testEdge].flip, edgeStack, vertexStack);

        getEdgeAndVertex(edges[testEdge].next, edgeStack, vertexStack);
    }
}

side voronoiCell::planeSide(VertexIndex vertex)
{
    float ans = (vertices[vertex].position - (particle.position + neighborParticle.position) / 2).Dot(particle.position - neighborParticle.position);

    if (std::abs(ans) < tolerance) {
        return incident;
    }
    if (ans > 0) {
        return inside;
    }
    if (ans < 0) {
        return outside;
    }

}



//Polyhedron
//    - Add a vertex in the middle of an edge
//    - Connect two vertices, assuming they are both on the same face
//    - Delete a vertex.


void voronoiCell::splitEdge(EdgeIndex edge, VertexIndex newVertex) {

    EdgeIndex newEdge = makeSelfLoopOnVertex(newVertex);
    EdgeIndex otherNewEdge = edges[newEdge].flip;
    EdgeIndex otherOldEdge = edges[edge].flip;
    edges[newEdge].creator = edges[edge].creator;
    edges[otherNewEdge].creator = edges[otherOldEdge].creator;
    std::swap(edges[edge], edges[newEdge]);
    std::swap(edges[otherOldEdge], edges[otherNewEdge]);
}

EdgeIndex voronoiCell::makeSelfLoopOnVertex(VertexIndex newVertex) {
    EdgeIndex edge = makeOneEdgeFace(newVertex);
    EdgeIndex back = makeOneEdgeFace(newVertex);
    setFlip(edge, back);
    return edge;
}

void voronoiCell::setFlip(EdgeIndex forward, EdgeIndex back) {
    edges[forward].flip = back;
    edges[back].flip = forward;
}

EdgeIndex voronoiCell::makeOneEdgeFace(VertexIndex vertex) {
    EdgeIndex edge = edges.create(HalfEdge(vertex));

    edges[edge].next = edge;
    return edge;
}

double voronoiCell::planeDist(VertexIndex vertex) {
    return std::abs((vertices[vertex].position.distanceTo(particle.position))
                     - (vertices[vertex].position.distanceTo(neighborParticle.position)));
}

// This is probably fairly close to being correct, but has a few nagging
// problems.
bool voronoiCell::findSomeIncidentEdge(EdgeIndex &returnEdge) {

    EdgeIndex testEdge = firstEdge;
    VertexIndex testVertex = edges[testEdge].target;

    side testSide = planeSide(testVertex);
    side flipSide = planeSide(edges[edges[testEdge].flip].target);



    double maxDist = planeDist(testVertex);


    // While we have not yet found an edge between two vertices not on
    // the same side of the plane
    while (testSide == flipSide) {

        // If the flip's target is closer, iterate on the flip
        if (planeDist(edges[edges[testEdge].flip].target) < maxDist) {
            testEdge = edges[testEdge].flip;
            testVertex = edges[testEdge].target;

            testSide = planeSide(testVertex);
            flipSide = planeSide(edges[edges[testEdge].flip].target);
            maxDist = planeDist(testVertex);

        } else {
            EdgeIndex nextEdge = edges[testEdge].next;

            // Otherwise, first check to see if the next is closer
            if (planeDist(edges[nextEdge].target) < maxDist) {

                testEdge = nextEdge;
                testVertex = edges[testEdge].target;

                testSide = planeSide(testVertex);

                flipSide = planeSide(edges[edges[testEdge].flip].target);
                maxDist = planeDist(testVertex);

            // and, if not, repeatedly examine the edges leaving the
            // testVertex until we find one that is closer
            } else {
                VertexIndex firstVertex = edges[nextEdge].target;

                // UNLESS we find an incident edge while doing so
                while (planeDist(edges[nextEdge].target) > maxDist &&
                       planeSide(edges[nextEdge].target) == planeSide(edges[edges[nextEdge].flip].target)) {

                    nextEdge = edges[edges[nextEdge].flip].next;

                    if (vertices[firstVertex].position == vertices[edges[nextEdge].target].position){

                        return false;
                    }
                }
                testEdge = nextEdge;
                testVertex = edges[testEdge].target;

                testSide = planeSide(testVertex);

                flipSide = planeSide(edges[edges[testEdge].flip].target);
                maxDist = planeDist(testVertex);
            }
        }
    }

    // We now should have an edge between vertices not on the same side of
    // the plane. We want to return an edge going from either inside->incident
    // or inside->outside (I think?). However, right now, this doesn't work exactly as
    // advertised.

    if (flipSide == inside) {
        returnEdge = testEdge;
        return true;
    } else {
        returnEdge = edges[testEdge].flip;
        return true;
    }
}

EdgeIndex voronoiCell::findNextIncidentEdge(EdgeIndex orig) {
    EdgeIndex cross = edges[orig].next;

    while (planeSide(edges[cross].target) == outside){
        cross = edges[cross].next;
    }

    return edges[cross].flip;
}

// Returns a Vector3 corresponding to the point at which the cutting plane
// intersects a given HalfEdge.
Vector3 voronoiCell::planeEdgeIntersect(EdgeIndex edge) {
    Vector3 p1 = vertices[edges[edge].target].position;
    Vector3 p2 = vertices[edges[edges[edge].flip].target].position;
    Vector3 intersectionPoint(0, 0, 0);
    Vector3 normal = neighborParticle.position - particle.position;
    Vector3 edgeVec = p2 - p1;
    intersectionPoint = p1 + edgeVec * (((particle.position + neighborParticle.position)/2
                        - p1).Dot(normal)) / (edgeVec.Dot(normal));

    return intersectionPoint;
}

// Creates an edge pair between a given pair of vertices
EdgeIndex voronoiCell::addEdgePair(VertexIndex vertex1, VertexIndex vertex2) {


    EdgeIndex forwardEdge = edges.create(HalfEdge(vertex2));

    EdgeIndex backEdge = edges.create(HalfEdge(vertex1));

    edges[forwardEdge].flip = backEdge;
    edges[backEdge].flip = forwardEdge;
    return forwardEdge;
}

EdgeIndex voronoiCell::maintainFirstEdge(EdgeIndex edge) {

    std::vector<EdgeIndex> testEdges;
    testEdges.push_back(edge);
    std::size_t i = 0;
    //     std::cerr << "Edge: " << edge << ", size: " << edges.size();
    //     std::cerr << ", Flip: " << edges[edge].flip  << ", vsize: " << vertices.size() << std::endl;
    while (planeSide(edges[edges[edge].flip].target) != inside) {
        // std::cerr << "Edge: " << edge << ", size: " << edges.size();
        // std::cerr << ", Flip: " << edges[edge].flip  << ", vsize: " << vertices.size() << std::endl;

        if (!edges[edge].seen) {

            testEdges.push_back(edges[edge].next);
            testEdges.push_back(edges[edge].flip);
            edges[edge].seen = true;
        }

        i++;
        edge = testEdges[i];
    }

    return edge;

    // if (planeSide(edge->flip->target) != outside) {

    //     return edge;
    // }
    // if (!edge->seen) {

    //     edge->seen = true;
    //     maintainFirstEdge(edge->next);
    //     maintainFirstEdge(edge->flip);
    // }
}

void voronoiCell::cutCell(const Particle& neighbor) {
    neighborParticle = neighbor;

    // We need to maintain firstEdge so that we can index into the Voronoi cell.
    // Specifically, we do not want its flip's target to be outside the plane,
    // or it will be deleted during the resulting cut. Therefore, we have the following
    // simple recursive routine for doing this:
    firstEdge = maintainFirstEdge(firstEdge);

    resetEdges(firstEdge);

    EdgeIndex startingIncidentEdge;

    // Determine whether there is an edge incident to the plane; if so, find
    // it (startingIncidentEdge is passed by reference)
    bool cuttable = findSomeIncidentEdge(startingIncidentEdge);

    if (!cuttable) {

        return;
    }

    EdgeIndex nextIncidentEdge = startingIncidentEdge;

    std::vector<EdgeIndex> outsideEdges;
    EdgeIndex prevIncidentEdge = nextIncidentEdge;
    EdgeIndex lastInsideEdge;

    do {

        // was inside

        if (planeSide(edges[nextIncidentEdge].target) == outside) {
            // Create a new vertex at the location of the intersection of
            // the crossing edge and the cutting plane.
            Vector3 vertexLoc = planeEdgeIntersect(nextIncidentEdge);
            VertexIndex newVertex = vertices.create(Vertex(vertexLoc));

            // was "splitCrossingEdge", but that was not defined -- should it
            // have been splitEdge?
            splitEdge(nextIncidentEdge, newVertex);

            // The first time, make sure to update the first incident edge
            if (prevIncidentEdge == nextIncidentEdge) {

                // For later creation of the final edge of the new face
                lastInsideEdge = edges[edges[edges[startingIncidentEdge].next].flip].next;
                startingIncidentEdge = edges[startingIncidentEdge].next;
            }
        }
        // Now it's an incident edge

        // we need to write getPrev as well. which is done. or consider storing
        // prev.
        // HalfEdge* prevInsideEdge = prevIncidentEdge->flip->getPrev();
        // HalfEdge* nextInsideEdge = nextIncidentEdge->next;

        // Create a new edge on the cutting plane (unless this is the first
        // incident edge, in which case do nothing)
        if (prevIncidentEdge != nextIncidentEdge) {
            EdgeIndex prevInsideEdge = prevIncidentEdge;
            EdgeIndex nextInsideEdge = edges[nextIncidentEdge].flip;

            // addEdgePair is now implemented, and thus is no longer a FIXME.
            EdgeIndex newInsideEdge = addEdgePair(edges[prevIncidentEdge].target, edges[nextIncidentEdge].target);
            edges[prevInsideEdge].next = newInsideEdge;
            edges[newInsideEdge].next = nextInsideEdge;
            edges[newInsideEdge].creator = edges[edges[newInsideEdge].next].creator;

            outsideEdges.push_back(edges[newInsideEdge].flip);
        }

        prevIncidentEdge = nextIncidentEdge;

        nextIncidentEdge = findNextIncidentEdge(prevIncidentEdge);

    } while(nextIncidentEdge != startingIncidentEdge && nextIncidentEdge != edges[edges[edges[startingIncidentEdge].flip].next].flip);

    // We have to manually add the last edge in the newly-created face.
    EdgeIndex prevInsideEdge = prevIncidentEdge;
    EdgeIndex nextInsideEdge = lastInsideEdge;

    // addEdgePair is now implemented, and thus is no longer a FIXME.
    EdgeIndex newInsideEdge = addEdgePair(edges[prevIncidentEdge].target, edges[edges[lastInsideEdge].flip].target);
    edges[prevInsideEdge].next = newInsideEdge;


    edges[newInsideEdge].next = nextInsideEdge;

    // All edges on a face should be associated with the
    // same neighbor particle.
    edges[newInsideEdge].creator = edges[edges[newInsideEdge].next].creator;

    outsideEdges.push_back(edges[newInsideEdge].flip);

    // Set the next of each edge in the new face, and set
    // the creator of each edge in the new face to the particle
    // we are cutting by.
    for (std::size_t i = 1; i < outsideEdges.size(); ++i) {
        edges[outsideEdges[i]].next = outsideEdges[i-1];
        edges[outsideEdges[i]].creator = &neighbor;



    }
    edges[outsideEdges.front()].next = outsideEdges.back();
    edges[outsideEdges.front()].creator = &neighbor;

    cleanUp(startingIncidentEdge);
}

void voronoiCell::cleanUp(EdgeIndex edge){

    std::stack<EdgeIndex>* deleteStackEdge = new std::stack<EdgeIndex>;

    std::stack<VertexIndex>* deleteStackVertex = new std::stack<VertexIndex>;

    deleteSearch(deleteStackEdge, deleteStackVertex, edge);

    while (!deleteStackEdge->empty()) {
        EdgeIndex topEdge = deleteStackEdge->top();
        deleteStackEdge->pop();

        edges.destroy(topEdge);

    }
    while (!deleteStackVertex->empty()) {
        VertexIndex topVertex = deleteStackVertex->top();

        deleteStackVertex->pop();

        vertices.destroy(topVertex);

    }
    delete deleteStackVertex;

    delete deleteStackEdge;

}

void voronoiCell::deleteSearch(std::stack<EdgeIndex>* deleteStackEdge,
                  std::stack<VertexIndex>* deleteStackVertex,
                  EdgeIndex edge){


        // Currently needs to watch out for problems with edges
        // within the cutting plane
    if ((planeSide(edges[edge].target) == outside || planeSide(edges[edges[edge].flip].target) == outside)
	  && !edges[edge].deleteFlag) {

        deleteStackEdge->push(edge);
        edges[edge].deleteFlag = true;
        deleteSearch(deleteStackEdge, deleteStackVertex, edges[edge].next);
        deleteSearch(deleteStackEdge, deleteStackVertex, edges[edge].flip);
    }

    if (planeSide(edges[edge].target) == outside && !vertices[edges[edge].target].deleteFlag) {

        deleteStackVertex->push(edges[edge].target);
        vertices[edges[edge].target].deleteFlag = true;
    }
}

// Calculates the volume of the current Voronoi cell, although note that
// a better implementation of this is possible if we store face_vertices in a
// persistent manner.
double voronoiCell::volume() {

    double volume = 0;

    // First, we need to get all the face vertices (if they haven't already been
    // calculated)
    if (faceVertices.size() == 0) {
        reset(firstEdge);

        getFaceVertex(firstEdge);
        resetEdgesAndVertices(firstEdge);

        resetEdges(firstEdge);

    }

    // We then calculate the volume through a tetrahedral decomposition.
    // First, we select one vertex, which will serve as one corner of all
    // of the tetrahedra.
    VertexIndex firstVertex = edges[firstEdge].target;

    // Next, loop over the face vertices, decomposing each face into
    // triangles, which can then be used to decompose the cell into tetrahedra.
    for (std::size_t i = 0; i < faceVertices.size(); ++i) {

        EdgeIndex firstBaseEdge = faceVertices[i].edges[0];
        EdgeIndex secondBaseEdge = faceVertices[i].edges[1];
        EdgeIndex thirdBaseEdge = faceVertices[i].edges[2];

        while (thirdBaseEdge != firstBaseEdge) {

            volume += tetVolume(firstVertex, edges[firstBaseEdge].target,
                                edges[secondBaseEdge].target, edges[thirdBaseEdge].target);

            secondBaseEdge = thirdBaseEdge;
            thirdBaseEdge = edges[thirdBaseEdge].next;
        }
    }


    return volume;
}

// Calculates the volume of a tetrahedron composed of four input vertices.
double voronoiCell::tetVolume(VertexIndex vertex1, VertexIndex vertex2,
                              VertexIndex vertex3, VertexIndex vertex4) {
    Vector3 firstVector = vertices[vertex2].position - vertices[vertex1].position;

    Vector3 secondVector = vertices[vertex3].position - vertices[vertex1].position;

    Vector3 thirdVector = vertices[vertex4].position - vertices[vertex1].position;

    return (std::abs(firstVector.Dot(secondVector.Cross(thirdVector)))/ 6);
}

// Takes by reference a vector that will become filled with
// the ids of the neighbors of this cell.
void voronoiCell::neighbors(std::vector<int> &v) {

    // We want to only call face_vertices when it
    // hasn't already been called
    if (faceVertices.size() == 0) {
        std::vector<int> face_v;
        face_vertices(face_v);
    }

    for (std::size_t i = 0; i < faceVertices.size(); ++i) {
        // We defaulted the creator to particle during initialization,
        // so this is a check for if the "neigbor" is the boundary
        // Should default to -1
        if (edges[faceVertices[i].edges[0]].creator->id != particle.id) {
            v.push_back(edges[faceVertices[i].edges[0]].creator->id);
        }
    }

}

void voronoiCell::face_areas(std::vector<double> &v) {
    // First, make the face vertices, if they haven't already been made
    if (faceVertices.size() == 0) {
        std::vector<int> face_v;
        face_vertices(face_v);
    }

    double area = 0.0;

    // Next, loop over the face vertices, decomposing each face into
    // triangles, whose areas will then be calculated
    for (std::size_t i = 0; i < faceVertices.size(); ++i) {
        area = 0.0;

        EdgeIndex firstBaseEdge = faceVertices[i].edges[0];
        EdgeIndex secondBaseEdge = faceVertices[i].edges[1];
        EdgeIndex thirdBaseEdge = faceVertices[i].edges[2];

        while (thirdBaseEdge != firstBaseEdge) {

            area += triArea(edges[firstBaseEdge].target,
                            edges[secondBaseEdge].target, edges[thirdBaseEdge].target);

            secondBaseEdge = thirdBaseEdge;
            thirdBaseEdge = edges[thirdBaseEdge].next;
        }

        v.push_back(area);
    }
}

// Calculates the area of a triangle composed of three input vertices.
double voronoiCell::triArea(VertexIndex vertex1, VertexIndex vertex2,
                              VertexIndex vertex3) {
    Vector3 firstVector = vertices[vertex2].position - vertices[vertex1].position;

    Vector3 secondVector = vertices[vertex3].position - vertices[vertex1].position;

    return (std::abs(firstVector.Cross(secondVector).distanceTo(Vector3(0, 0, 0)))/ 2);
}

void voronoiCell::face_vertices(std::vector<int> &v) {

    resetEdgesAndVertices(firstEdge);
    getFaceVertex(firstEdge);

    for (std::size_t i = 0; i < faceVertices.size(); ++i) {
        // THIS IS NOT THE FORMAT THIS SHOULD BE OUTPUT IN, since unfortunately
        // the voro++ output relies on a well-defined vertex ordering. We need to
        // make a decision over whether or not we will permanently store the
        // vertices, faceVertices after first calculating them.
        v.push_back(i);
    }

    // We eventually need to reset all the seen flags at the end of this. or,
    // more likely, we need to optimize this procedure.
    resetEdges(firstEdge);
};

void voronoiCell::getFaceVertex(EdgeIndex testEdge) {

    if (!edges[testEdge].seen) {

        FaceVertex newFaceVertex;

        newFaceVertex.edges.push_back(testEdge);


        EdgeIndex otherEdge = edges[testEdge].next;

        edges[testEdge].seen = true;

        while (otherEdge != testEdge) {

            newFaceVertex.edges.push_back(otherEdge);
            edges[otherEdge].seen = true;

            otherEdge = edges[otherEdge].next;

        }

        faceVertices.push_back(newFaceVertex);

        getFaceVertex(edges[testEdge].flip);

        otherEdge = edges[testEdge].next;
        while (otherEdge != testEdge) {
            getFaceVertex(edges[otherEdge].flip);
            otherEdge = edges[otherEdge].next;
        }
    }
};

// CustomInterface IO
void voronoiCell::computeVertices(std::vector<double> &v) {
    std::vector<VertexIndex> vs;
    getVertex(firstEdge, vs);
    resetEdgesAndVertices(firstEdge);

    for (std::size_t i = 0; i < vs.size(); ++i) {
        v.push_back(vertices[vs[i]].position.X);
        v.push_back(vertices[vs[i]].position.Y);
        v.push_back(vertices[vs[i]].position.Z);
    }

    // we eventually need to reset all the seen flags at the end of this. or,
    // more likely, we need to optimize this procedure.
    resetEdgesAndVertices(firstEdge);
};

void voronoiCell::getVertex(EdgeIndex testEdge, std::vector<VertexIndex> &vs) {

    if (edges[testEdge].seen) {
        edges[testEdge].seen = true;

        if (!vertices[edges[testEdge].target].seen) {
            vertices[edges[testEdge].target].seen = true;
            vs.push_back(edges[testEdge].target);
        }

        getVertex(edges[testEdge].flip, vs);
        getVertex(edges[testEdge].next, vs);
    }
};

// Resets the seen flags on all edges
void voronoiCell::resetEdges(EdgeIndex edge) {
    if (edges[edge].seen) {

        edges[edge].seen = false;
        resetEdges(edges[edge].flip);
        resetEdges(edges[edge].next);
    }
};

// Resets the seen flags on all edges and vertices
void voronoiCell::resetEdgesAndVertices(EdgeIndex edge) {
    if (edges[edge].seen) {

        edges[edge].seen = false;

        // Reset the associated vertex, if we need to
        if (vertices[edges[edge].target].seen) {

            vertices[edges[edge].target].seen = false;
        }

        resetEdgesAndVertices(edges[edge].flip);
        resetEdgesAndVertices(edges[edge].next);
    }
};

void voronoiCell::reset(EdgeIndex edge){

    std::stack<EdgeIndex>* seenStackEdge = new std::stack<EdgeIndex>;

    std::stack<VertexIndex>* seenStackVertex = new std::stack<VertexIndex>;

    seenSearch(seenStackEdge, seenStackVertex, edge);

    while (!seenStackEdge->empty()) {
        EdgeIndex topEdge = seenStackEdge->top();
        seenStackEdge->pop();

        edges[topEdge].seen = false;

    }
    while (!seenStackVertex->empty()) {
        VertexIndex topVertex = seenStackVertex->top();

        seenStackVertex->pop();

        vertices[topVertex].seen = false;
    }
    delete seenStackVertex;

    delete seenStackEdge;

};

void voronoiCell::seenSearch(std::stack<EdgeIndex>* seenStackEdge,
                  std::stack<VertexIndex>* seenStackVertex,
                  EdgeIndex edge){

    // Currently needs to watch out for problems with edges
    // within the cutting plane
    if (!edges[edge].seen) {
        seenStackEdge->push(edge);
        edges[edge].seen = true;
        seenSearch(seenStackEdge, seenStackVertex, edges[edge].next);
        seenSearch(seenStackEdge, seenStackVertex, edges[edge].flip);
    }

    if (!vertices[edges[edge].target].seen) {
        seenStackVertex->push(edges[edge].target);
        vertices[edges[edge].target].seen = true;
    }
};


/** Outputs the edges of the Voronoi cell in gnuplot format to an output stream.
 * \param[in] (x,y,z) a displacement vector to be added to the cell's position.
 * \param[in] fp a file handle to write to. */
// Shamelessly stolen from voro++
// Remember to resetEdges after calling this
void voronoiCell::drawGnuplot(double dispX, double dispY, double dispZ, FILE* fp, EdgeIndex edge) {
    EdgeIndex currentEdge = edges[edge].next;
    Vertex* target = &vertices[edges[edges[currentEdge].flip].target];
    fprintf(fp, "%g, %g, %g\n", target->position.X + dispX,
                                target->position.Y + dispY,
                                target->position.Z + dispZ);

    // Loop around the current face, drawing all edges.
    while (currentEdge != edge) {
        target = &vertices[edges[currentEdge].target];
        fprintf(fp, "%g, %g, %g\n", target->position.X + dispX,
                                    target->position.Y + dispY,
                                    target->position.Z + dispZ);

        currentEdge = edges[currentEdge].next;
    }

    fprintf(fp, "%g, %g, %g\n", target->position.X + dispX,
                                target->position.Y + dispY,
                                target->position.Z + dispZ);

    fputs("\n\n", fp);

    // Loop around the face again, recursing on faces that have not yet been drawn
    edges[currentEdge].seen = true;
    if (!edges[edges[currentEdge].flip].seen) {
        drawGnuplot(dispX, dispY, dispZ, fp, edges[currentEdge].flip);
    }
    currentEdge = edges[currentEdge].next;

    while (currentEdge != edge) {
        edges[currentEdge].seen = true;
        if (!edges[edges[currentEdge].flip].seen) {
            drawGnuplot(dispX, dispY, dispZ, fp, edges[currentEdge].flip);
        }
        currentEdge = edges[currentEdge].next;
    }
}

std::size_t voronoiCell::get_memory_usage(){
    std::size_t memory = 0;
    memory += sizeof(FaceVertex) * faceVertices.capacity();
    memory += sizeof(voronoiCell);
    memory += edges.get_memory_usage() + vertices.get_memory_usage();
    std::cout << "Edges: " << edges.get_memory_usage() << std::endl;
    std::cout << "Vertices: " << vertices.get_memory_usage() << std::endl;
    std::cout << "voronoi cell: " << memory << std::endl;
    return memory;
}


// -----------------------------CELL CONTAINER------------------------------ //

cellContainer::cellContainer(std::vector<Particle> parts, double defaultLen):
    defaultLength(defaultLen), x_min(-defaultLen/2), x_max(defaultLen/2), y_min(-defaultLen/2), y_max(defaultLen/2),
    z_min(-defaultLen/2), z_max(defaultLen/2) {
    particles = parts;
}

cellContainer::cellContainer(std::vector<Particle> parts, double defaultLen,
                             double xmin, double xmax, double ymin, double ymax,
                             double zmin, double zmax) :
    defaultLength(defaultLen), x_min(xmin), x_max(xmax), y_min(ymin), y_max(ymax),
    z_min(zmin), z_max(zmax) {
    particles = parts;
}

voronoiCell cellContainer::makeCell(Particle particle) {


    // Calculate the initial max radius of a cube
    double maxRadius = sqrt(3) * defaultLength;


    // Initialize the voronoi cell as a cube...
    voronoiCell cell = voronoiCell("cube", defaultLength, particle, maxRadius,
        x_min, x_max, y_min, y_max, z_min, z_max);

    // loop through our array of particles and cut the cell with planes
    // associated with each possible neighbor.

    std::array<double, 3> point = {{particle.position.X, particle.position.Y, particle.position.Z}};

    std::vector<Particle*> NeighborParticles;

    sds.neighborQuery(point, defaultLength, &NeighborParticles);


    for (size_t i = 0; i < NeighborParticles.size(); ++i) {
        // This condition may need to be fixed.

        if ((NeighborParticles[i]->position != particle.position) && (particle.position.distanceTo(NeighborParticles[i]->position) < maxRadius)) {

            cell.cutCell(*NeighborParticles[i]);
        }
    }

    return cell;
}

cellContainer::~cellContainer() {

    // while (!cells.empty()) {
    //     delete cells.back();

    //     cells.pop_back();
    // }
}

double cellContainer::sum_cell_volumes() {
    double sum = 0;

    // Need to put in an initialize function. DO NOT KEEP HERE
    sds.initialize(&particles[0], &particles[particles.size()]);

    // if(!calculated) {
        voronoiCell c;
        for(unsigned int i = 0; i < particles.size(); ++i) {

            c = makeCell(particles[i]);
            double vol = c.volume();
            sum += vol;//c.volume();
            // std::cout << "Cell volume: " << vol << std::endl;
            // if(c.faceVertices.size() > 0) {
            //     std::cout << "FACEVERTEX SIZE GREATER THAN 0: " << c.faceVertices.size() << std::endl;
            // }
        }
    // }
    // else {

    //     for(unsigned int i = 0; i < particles.size(); ++i) {
    //         sum += c.volume();
    //     }
    // }
    return sum;
}


void cellContainer::put(int id, double px, double py, double pz) {
    particles.push_back(Particle(id, px, py, pz));
}

// double cellContainer::findMaxNeighDist() {
//     double maxDist = 0;
//     for(voronoiCell* cellPointer : cells) {
//         std::vector<int> neighborIndexes;
//         cellPointer->neighbors(neighborIndexes);
//         for (int i = 0; i < neighborIndexes.size(); ++i) {
//             double dist = cellPointer->particle.position.distanceTo(particles[neighborIndexes[i]].position);
//             if(dist > maxDist) {
//                 maxDist = dist;
//             }
//         }
//     }
//     return maxDist;

// }

// std::size_t cellContainer::get_memory_usage() {
//     std::size_t memory = 0;
//     memory += sizeof(Particle) * particles.capacity();
//     memory += sizeof(cellContainer);
//     std::cout << "Cell container: " << memory << std::endl;
//     for(auto cell : cells) {
//         memory += cell->get_memory_usage();
//     }
//     return memory;
// }

#endif
