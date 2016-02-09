// added just for one use of abs; let me know if there is a better way
// to do this
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

Particle::Particle( void ){
    id = -1;
}


// halfedge constructors

HalfEdge::HalfEdge( void ) {
    deleteFlag = false;
    seen = false;
}

HalfEdge::HalfEdge(Vertex* vertex) {
    target = vertex;
    deleteFlag = false;
    seen = false;
}

HalfEdge::HalfEdge(Vertex* vertex, Particle* neighbor) {
    target = vertex;
    deleteFlag = false;
    seen = false;
    creator = neighbor;
}

// This constructor builds a HalfEdge with a specified target and next.
HalfEdge::HalfEdge(Vertex* vertex, HalfEdge* edge2, Particle* neighbor) {
    target = vertex;
    next = edge2;
    deleteFlag = false;
    seen = false;
    creator = neighbor;
}

HalfEdge::HalfEdge(Vertex* vertex, HalfEdge* edge1, HalfEdge* edge2, Particle* neighbor) {
    target = vertex;
    flip = edge1;
    next = edge2;
    deleteFlag = false;
    seen = false;
    creator = neighbor;
}

// Intended to get the previous edge, should be used rarely.
HalfEdge* HalfEdge::getPrev() {
    HalfEdge* other = this;

    //TODO: Optimize
    // HalfEdge* other = next->next;
    while (other->next != this) {
        other = other->next;
    }
    return other;
}

// HalfEdge* firstEdge = getFirstEdge();

voronoiCell::voronoiCell(std::string shape, double length, Particle seedParticle, 
                         double maxRadius, double x_min, double x_max,
                         double y_min, double y_max, double z_min, double z_max) {
    particle = seedParticle;
    maxRadius = maxRadius;
    faceVertices = std::vector<FaceVertex*>();

    if (shape == "cube") {
        // Create half-edge arrays corresponding to the faces of the cube.
        HalfEdge *plusZ[4], *minusZ[4], *plusY[4], *minusY[4], *plusX[4], *minusX[4];

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

        Vertex* oct1 = new Vertex(x_max,
                                  y_max,
                                  z_max);
        ////std::cout << "Allocating vertex: initializing cube" << std::endl;
        

        Vertex* oct2 = new Vertex(x_min,
                                  y_max,
                                  z_max);
        ////std::cout << "Allocating vertex: initializing cube" << std::endl;
        

        Vertex* oct3 = new Vertex(x_min,
                                  y_min,
                                  z_max);
        ////std::cout << "Allocating vertex: initializing cube" << std::endl;
        

        Vertex* oct4 = new Vertex(x_max,
                                  y_min,
                                  z_max);
        ////std::cout << "Allocating vertex: initializing cube" << std::endl;
        

        Vertex* oct5 = new Vertex(x_max,
                                  y_max,
                                  z_min);
        ////std::cout << "Allocating vertex: initializing cube" << std::endl;
        

        Vertex* oct6 = new Vertex(x_min,
                                  y_max,
                                  z_min);
        ////std::cout << "Allocating vertex: initializing cube" << std::endl;
        

        Vertex* oct7 = new Vertex(x_min,
                                  y_min,
                                  z_min);
        ////std::cout << "Allocating vertex: initializing cube" << std::endl;
        

        Vertex* oct8 = new Vertex(x_max,
                                  y_min,
                                  z_min);
        ////std::cout << "Allocating vertex: initializing cube" << std::endl;
        


        // Vertex* oct1 = new Vertex(particle.X + length/2,
        //                           particle.Y + length/2,
        //                           particle.Z + length/2);
        // if (particle.X + length/2 > x_max) {
        //     oct1->position.X = x_max;
        // }
        // if (particle.Y + length/2 > y_max) {
        //     oct1->position.Y = y_max;
        // }
        // if (particle.Z + length/2 > z_max) {
        //     oct1->position.Z = z_max;
        // }

        // Vertex* oct2 = new Vertex(particle.X - length/2,
        //                           particle.Y + length/2,
        //                           particle.Z + length/2);
        // if (particle.X - length/2 > x_min) {
        //     oct2->position.X = x_min;
        // }
        // if (particle.Y + length/2 > y_max) {
        //     oct2->position.Y = y_max;
        // }
        // if (particle.Z + length/2 > z_max) {
        //     oct2->position.Z = z_max;
        // }

        // Vertex* oct3 = new Vertex(particle.X - length/2,
        //                           particle.Y - length/2,
        //                           particle.Z + length/2);
        // if (particle.X - length/2 > x_min) {
        //     oct3->position.X = x_min;
        // }
        // if (particle.Y - length/2 > y_min) {
        //     oct3->position.Y = y_min;
        // }
        // if (particle.Z + length/2 > z_max) {
        //     oct3->position.Z = z_max;
        // }

        // Vertex* oct4 = new Vertex(particle.X + length/2,
        //                           particle.Y - length/2,
        //                           particle.Z + length/2);
        // if (particle.X + length/2 > x_max) {
        //     oct4->position.X = x_max;
        // }
        // if (particle.Y - length/2 > y_min) {
        //     oct4->position.Y = y_min;
        // }
        // if (particle.Z + length/2 > z_max) {
        //     oct4->position.Z = z_max;
        // }

        // Vertex* oct5 = new Vertex(particle.X + length/2,
        //                           particle.Y + length/2,
        //                           particle.Z - length/2);
        // if (particle.X + length/2 > x_max) {
        //     oct5->position.X = x_max;
        // }
        // if (particle.Y + length/2 > y_max) {
        //     oct5->position.Y = y_max;
        // }
        // if (particle.Z - length/2 > z_min) {
        //     oct5->position.Z = z_min;
        // }

        // Vertex* oct6 = new Vertex(particle.X - length/2,
        //                           particle.Y + length/2,
        //                           particle.Z - length/2);
        // if (particle.X - length/2 > x_min) {
        //     oct6->position.X = x_min;
        // }
        // if (particle.Y + length/2 > y_max) {
        //     oct6->position.Y = y_max;
        // }
        // if (particle.Z - length/2 > z_min) {
        //     oct6->position.Z = z_min;
        // }

        // Vertex* oct7 = new Vertex(particle.X - length/2,
        //                           particle.Y - length/2,
        //                           particle.Z - length/2);
        // if (particle.X - length/2 > x_min) {
        //     oct7->position.X = x_min;
        // }
        // if (particle.Y - length/2 > y_min) {
        //     oct7->position.Y = y_min;
        // }
        // if (particle.Z - length/2 > z_min) {
        //     oct7->position.Z = z_min;
        // }

        // Vertex* oct8 = new Vertex(particle.X + length/2,
        //                           particle.Y - length/2,
        //                           particle.Z - length/2);
        // if (particle.X + length/2 > x_max) {
        //     oct8->position.X = x_max;
        // }
        // if (particle.Y - length/2 > y_min) {
        //     oct8->position.Y = y_min;
        // }
        // if (particle.Z - length/2 > z_min) {
        //     oct8->position.Z = z_min;
        // }

        // We will construct each face by beginning with the edge whose target
        // is in the lowest-numbered quadrant, then continue around the face
        // by iteratively creating the prev of the most recently created edge.
        //
        // Traversing the nexts around a face leads if you are looking
        // directly at it should lead to clockwise traversal of the face.

        // Note also our choice of default initialization for the particle
        // associated with each HalfEdge; since no planes have cut this cell
        // yet, we associate it with the particle at the center of this cell.
        plusZ[0] = new HalfEdge(oct1, &particle);             // O2 to O1
        ////std::cout << "Allocating halfedge: initializing cube" << std::endl;
        plusZ[1] = new HalfEdge(oct2, plusZ[0], &particle);   // O3 to O2
        ////std::cout << "Allocating halfedge: initializing cube" << std::endl;
        plusZ[2] = new HalfEdge(oct3, plusZ[1], &particle);   // O4 to O3
        ////std::cout << "Allocating halfedge: initializing cube" << std::endl;
        plusZ[3] = new HalfEdge(oct4, plusZ[2], &particle);   // O1 to O4
        ////std::cout << "Allocating halfedge: initializing cube" << std::endl;
        plusZ[0]->next = plusZ[3];

        minusZ[0] = new HalfEdge(oct5, &particle);            // O8 to O5
        ////std::cout << "Allocating halfedge: initializing cube" << std::endl;
        minusZ[1] = new HalfEdge(oct8, minusZ[0], &particle); // O7 to O8
        ////std::cout << "Allocating halfedge: initializing cube" << std::endl;
        minusZ[2] = new HalfEdge(oct7, minusZ[1], &particle); // O6 to O7
        ////std::cout << "Allocating halfedge: initializing cube" << std::endl;
        minusZ[3] = new HalfEdge(oct6, minusZ[2], &particle); // O5 to O6
        ////std::cout << "Allocating halfedge: initializing cube" << std::endl;
        minusZ[0]->next = minusZ[3];

        plusY[0] = new HalfEdge(oct1, &particle);             // O5 to O1
        ////std::cout << "Allocating halfedge: initializing cube" << std::endl;
        plusY[1] = new HalfEdge(oct5, plusY[0], &particle);   // O6 to O5
        ////std::cout << "Allocating halfedge: initializing cube" << std::endl;
        plusY[2] = new HalfEdge(oct6, plusY[1], &particle);   // O2 to O6
        ////std::cout << "Allocating halfedge: initializing cube" << std::endl;
        plusY[3] = new HalfEdge(oct2, plusY[2], &particle);   // O1 to O2
        ////std::cout << "Allocating halfedge: initializing cube" << std::endl;
        plusY[0]->next = plusY[3];

        minusY[0] = new HalfEdge(oct3, &particle);            // O7 to O3
        ////std::cout << "Allocating halfedge: initializing cube" << std::endl;
        minusY[1] = new HalfEdge(oct7, minusY[0], &particle); // O8 to O7
        ////std::cout << "Allocating halfedge: initializing cube" << std::endl;
        minusY[2] = new HalfEdge(oct8, minusY[1], &particle); // O4 to O8
        ////std::cout << "Allocating halfedge: initializing cube" << std::endl;
        minusY[3] = new HalfEdge(oct4, minusY[2], &particle); // O3 to O4
        ////std::cout << "Allocating halfedge: initializing cube" << std::endl;
        minusY[0]->next = minusY[3];

        plusX[0] = new HalfEdge(oct1, &particle);             // O4 to O1 
        ////std::cout << "Allocating halfedge: initializing cube" << std::endl;
        plusX[1] = new HalfEdge(oct4, plusX[0], &particle);   // O8 to O4 
        ////std::cout << "Allocating halfedge: initializing cube" << std::endl;
        plusX[2] = new HalfEdge(oct8, plusX[1], &particle);   // O5 to O8 
        ////std::cout << "Allocating halfedge: initializing cube" << std::endl;
        plusX[3] = new HalfEdge(oct5, plusX[2], &particle);   // O1 to O5 
        ////std::cout << "Allocating halfedge: initializing cube" << std::endl;
        plusX[0]->next = plusX[3];

        minusX[0] = new HalfEdge(oct2, &particle);            // O6 to O2
        ////std::cout << "Allocating halfedge: initializing cube" << std::endl;
        minusX[1] = new HalfEdge(oct6, minusX[0], &particle); // O7 to O6
        ////std::cout << "Allocating halfedge: initializing cube" << std::endl;
        minusX[2] = new HalfEdge(oct7, minusX[1], &particle); // O3 to O7
        ////std::cout << "Allocating halfedge: initializing cube" << std::endl;
        minusX[3] = new HalfEdge(oct3, minusX[2], &particle); // O2 to O3
        ////std::cout << "Allocating halfedge: initializing cube" << std::endl;
        minusX[0]->next = minusX[3];

        // Now we have to set the flips manually
        plusZ[0]->flip = plusY[3];
        plusZ[1]->flip = minusX[3];
        plusZ[2]->flip = minusY[3];
        plusZ[3]->flip = plusX[0];

        minusZ[0]->flip = plusX[2];
        minusZ[1]->flip = minusY[1];
        minusZ[2]->flip = minusX[1];
        minusZ[3]->flip = plusY[1];

        plusY[0]->flip = plusX[3];
        plusY[1]->flip = minusZ[3];
        plusY[2]->flip = minusX[0];
        plusY[3]->flip = plusZ[0];

        minusY[0]->flip = minusX[2];
        minusY[1]->flip = minusZ[1];
        minusY[2]->flip = plusX[1];
        minusY[3]->flip = plusZ[2];

        plusX[0]->flip = plusZ[3];
        plusX[1]->flip = minusY[2];
        plusX[2]->flip = minusZ[0];
        plusX[3]->flip = plusY[0];

        minusX[0]->flip = plusY[2];
        minusX[1]->flip = minusZ[2];
        minusX[2]->flip = minusY[0];
        minusX[3]->flip = plusZ[1];

        //// std::cout << "Hi! My name is O2 to O1 plusZ[0] " << plusZ[0] << " with target "
        //       << plusZ[0]->target->position.X << ", "
        //       << plusZ[0]->target->position.Y << ", "
        //       << plusZ[0]->target->position.Z << " and flip target "
        //       << plusZ[0]->flip->target->position.X << ", "
        //       << plusZ[0]->flip->target->position.Y << ", "
        //       << plusZ[0]->flip->target->position.Z << std::endl;

        //// std::cout << "The next edge of " << plusZ[0] << " has target "
        //       << plusZ[0]->next->target->position.X << ", "
        //       << plusZ[0]->next->target->position.Y << ", "
        //       << plusZ[0]->next->target->position.Z << " and flip target "
        //       << plusZ[0]->next->flip->target->position.X << ", "
        //       << plusZ[0]->next->flip->target->position.Y << ", "
        //       << plusZ[0]->next->flip->target->position.Z << std::endl;

        //// std::cout << "Hi! My name is O5 to O1 plusY[0] " << plusY[0] << " with target "
        //       << plusY[0]->target->position.X << ", "
        //       << plusY[0]->target->position.Y << ", "
        //       << plusY[0]->target->position.Z << " and flip target "
        //       << plusY[0]->flip->target->position.X << ", "
        //       << plusY[0]->flip->target->position.Y << ", "
        //       << plusY[0]->flip->target->position.Z << std::endl;

        //// std::cout << "The next edge of " << plusY[0] << " has target "
        //       << plusY[0]->next->target->position.X << ", "
        //       << plusY[0]->next->target->position.Y << ", "
        //       << plusY[0]->next->target->position.Z << " and flip target "
        //       << plusY[0]->next->flip->target->position.X << ", "
        //       << plusY[0]->next->flip->target->position.Y << ", "
        //       << plusY[0]->next->flip->target->position.Z << std::endl;

        //// std::cout << "Hi! My name is O4 to O1 plusX[0] " << plusX[0] << " with target "
        //       << plusX[0]->target->position.X << ", "
        //       << plusX[0]->target->position.Y << ", "
        //       << plusX[0]->target->position.Z << " and flip target "
        //       << plusX[0]->flip->target->position.X << ", "
        //       << plusX[0]->flip->target->position.Y << ", "
        //       << plusX[0]->flip->target->position.Z << std::endl;

        //// std::cout << "The next edge of " << plusX[0] << " has target "
        //       << plusX[0]->next->target->position.X << ", "
        //       << plusX[0]->next->target->position.Y << ", "
        //       << plusX[0]->next->target->position.Z << " and flip target "
        //       << plusX[0]->next->flip->target->position.X << ", "
        //       << plusX[0]->next->flip->target->position.Y << ", "
        //       << plusX[0]->next->flip->target->position.Z << std::endl;

        // Set firstEdge to our favorite edge.
        //// std::cout << "Setting firstEdge..." << std::endl;
        firstEdge = plusX[0];
        //// std::cout << "First edge creator initialized as: " << firstEdge->creator->id <<std::endl;
    }
}

voronoiCell::~voronoiCell() {
    //// std::cout << "Entering voronoiCell destructor" << std::endl;
    resetEdgesAndVertices(firstEdge);
    std::stack<HalfEdge*> edgeStack;
    std::stack<Vertex*> vertexStack;

    getEdgeAndVertex(firstEdge, edgeStack, vertexStack);
    resetEdgesAndVertices(firstEdge);
    //// std::cout << "Filled delete stacks" << std::endl;

    while (!edgeStack.empty()) {
        delete edgeStack.top();
        ////std::cout << "Deleting halfedge: ~voronoiCell" << std::endl;
        edgeStack.pop();
    }
    while (!vertexStack.empty()) {
        delete vertexStack.top();
        ////std::cout << "Deleting vertex: ~voronoiCell" << std::endl;
        vertexStack.pop();
    }

    while (!faceVertices.empty()) {
        delete faceVertices.back();
        ////std::cout << "Deleting facevertex: ~voronoiCell" << std::endl;
        faceVertices.pop_back();
    }
}

void voronoiCell::getEdgeAndVertex(HalfEdge* testEdge, std::stack<HalfEdge*> &edgeStack,
                                   std::stack<Vertex*> &vertexStack) {
    if (!testEdge->seen) {
        //// std::cout << "Found edge to delete" << std::endl;
        testEdge->seen = true;
        //// std::cout << "seen" << std::endl;
        edgeStack.push(testEdge);
        
        if (!testEdge->target->seen) {
            testEdge->target->seen = true;
            vertexStack.push(testEdge->target);
        }
        
        //// std::cout << "pushback" << std::endl;
        //// std::cout << "flip of " << testEdge << std::endl;
        getEdgeAndVertex(testEdge->flip, edgeStack, vertexStack);
        //// std::cout << "next of " << testEdge << std::endl;
        getEdgeAndVertex(testEdge->next, edgeStack, vertexStack);
    }
}

side voronoiCell::planeSide(Vertex* vertex)
{
    float ans = (vertex->position - (particle.position + neighborParticle.position) / 2).Dot(particle.position - neighborParticle.position);
    
    if (std::abs(ans) < tolerance) {
        //// std::cout << "Found incident: " << vertex->position.X << " " << vertex->position.Y << " " 
        // << vertex->position.Z << " " << std::abs(ans) << " " << ans << std::endl;
        return incident;
    }
    if (ans > 0) {
        return inside;
    }
    if (ans < 0) {
        return outside;
    }
    
    // else {
    //     return incident;
    // }
}



//Polyhedron
//    - Add a vertex in the middle of an edge
//    - Connect two vertices, assuming they are both on the same face
//    - Delete a vertex.


void voronoiCell::splitEdge(HalfEdge* edge, Vertex* newVertex) {
    //// std::cout << "Splitting an edge..." << std::endl;
    HalfEdge* newEdge = makeSelfLoopOnVertex(newVertex);
    HalfEdge* otherNewEdge = newEdge->flip;
    HalfEdge* otherOldEdge = edge->flip;
    newEdge->creator = edge->creator;
    otherNewEdge->creator = otherOldEdge->creator;
    std::swap(*edge, *newEdge);
    std::swap(*otherOldEdge, *otherNewEdge);

    //// std::cout << "Made an edge " << newEdge << " with target "
    //           << newEdge->target->position.X << ", "
    //           << newEdge->target->position.Y << ", "
    //           << newEdge->target->position.Z << " and flip target "
    //           << newEdge->flip->target->position.X << ", "
    //           << newEdge->flip->target->position.Y << ", "
    //           << newEdge->flip->target->position.Z << std::endl;

    //// std::cout << "The next edge of " << newEdge << " has target "
    //           << newEdge->next->target->position.X << ", "
    //           << newEdge->next->target->position.Y << ", "
    //           << newEdge->next->target->position.Z << " and flip target "
    //           << newEdge->next->flip->target->position.X << ", "
    //           << newEdge->next->flip->target->position.Y << ", "
    //           << newEdge->next->flip->target->position.Z << std::endl;

    //// std::cout << "Made an edge " << otherNewEdge << " with target "
    //           << otherNewEdge->target->position.X << ", "
    //           << otherNewEdge->target->position.Y << ", "
    //           << otherNewEdge->target->position.Z << " and flip target "
    //           << otherNewEdge->flip->target->position.X << ", "
    //           << otherNewEdge->flip->target->position.Y << ", "
    //           << otherNewEdge->flip->target->position.Z << std::endl;

    //// std::cout << "The next edge of " << otherNewEdge << " has target "
    //           << otherNewEdge->next->target->position.X << ", "
    //           << otherNewEdge->next->target->position.Y << ", "
    //           << otherNewEdge->next->target->position.Z << " and flip target "
    //           << otherNewEdge->next->flip->target->position.X << ", "
    //           << otherNewEdge->next->flip->target->position.Y << ", "
    //           << otherNewEdge->next->flip->target->position.Z << std::endl;

}

HalfEdge* voronoiCell::makeSelfLoopOnVertex(Vertex* newVertex) {
    HalfEdge* edge = makeOneEdgeFace(newVertex);
    HalfEdge* back = makeOneEdgeFace(newVertex);
    setFlip(edge, back);
    return edge;
}

void voronoiCell::setFlip(HalfEdge* forward, HalfEdge* back) {
    forward->flip = back;
    back->flip = forward;
}

HalfEdge* voronoiCell::makeOneEdgeFace(Vertex* vertex) {
    HalfEdge* edge = new HalfEdge(vertex);
    ////std::cout << "Allocating halfedge: makeOneEdgeFace" << std::endl;
    edge->next = edge;
    return edge;
}

// I think this works?
double voronoiCell::planeDist(Vertex* vertex) {
    return std::abs((vertex->position.distanceTo(particle.position))
                     - (vertex->position.distanceTo(neighborParticle.position)));
}

// This is probably fairly close to being correct, but has a few nagging
// problems.
bool voronoiCell::findSomeIncidentEdge(HalfEdge* &returnEdge) {
    //// std::cout << "Finding some incident edge..." << std::endl;
    HalfEdge* testEdge = firstEdge;
    Vertex* testVertex = testEdge->target;

    side testSide = planeSide(testVertex);
    side flipSide = planeSide(testEdge->flip->target);

    

    double maxDist = planeDist(testVertex);
    //// std::cout << "Found maxDist: " << maxDist << std::endl;

    // While we have not yet found an edge between two vertices not on
    // the same side of the plane
    while (testSide == flipSide) {
        //// std::cout << "Testing edge " << testEdge << " with target "
        //           << testEdge->target->position.X << ", "
        //           << testEdge->target->position.Y << ", "
        //           << testEdge->target->position.Z << " and flip target "
        //           << testEdge->flip->target->position.X << ", "
        //           << testEdge->flip->target->position.Y << ", "
        //           << testEdge->flip->target->position.Z 
        //           << " and target side " << testSide
        //           << " and flip target side " << flipSide << std::endl;
        //// std::cout << "testSide = flipSide" << std::endl;

        // If the flip's target is closer, iterate on the flip
        if (planeDist(testEdge->flip->target) < maxDist) {
            testEdge = testEdge->flip;
            testVertex = testEdge->target;

            testSide = planeSide(testVertex);
            flipSide = planeSide(testEdge->flip->target);
            maxDist = planeDist(testVertex);  
            //// std::cout << "flip is closer" << std::endl; 

        } else {
            HalfEdge* nextEdge = testEdge->next;

            // Otherwise, first check to see if the next is closer
            if (planeDist(nextEdge->target) < maxDist) {
                //// std::cout << "next is closer" << std::endl; 
                testEdge = nextEdge;
                testVertex = testEdge->target;

                testSide = planeSide(testVertex);
                flipSide = planeSide(testEdge->flip->target);
                maxDist = planeDist(testVertex);    

            // and, if not, repeatedly examine the edges leaving the
            // testVertex until we find one that is closer
            } else {
                Vertex* firstVertex = nextEdge->target;
                //// std::cout << "examining edge leaving testVertex " << maxDist << std::endl; 
                // UNLESS we find an incident edge while doing so
                //// std::cout << "First vertex: " << firstVertex->position.X << " "<<firstVertex->position.Y<<" "<<firstVertex->position.Z<<std::endl;

                while (planeDist(nextEdge->target) > maxDist &&
                       planeSide(nextEdge->target) == planeSide(nextEdge->flip->target)) {
                  //  // std::cout << "we're in a loop! " << planeDist(nextEdge->target) << " " << maxDist << " " 
                  //   << planeSide(nextEdge->target) << " " << planeSide(nextEdge->flip->target) << " with target "
                  // << nextEdge->target->position.X << ", "
                  // << nextEdge->target->position.Y << ", "
                  // << nextEdge->target->position.Z << std::endl; 
                  //  // std::cout << planeDist(nextEdge->target) << std::endl;

                    nextEdge = nextEdge->flip->next;

                    if (firstVertex->position == nextEdge->target->position){
                        //// std::cout << "Already seen this edge. Exiting findSomeIncidentEdge" << std::endl;
                        return false;
                    }
                }
                testEdge = nextEdge;
                testVertex = testEdge->target;

                testSide = planeSide(testVertex);
                flipSide = planeSide(testEdge->flip->target);
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
        returnEdge = testEdge->flip;
        return true;
    }
}

HalfEdge* voronoiCell::findNextIncidentEdge(HalfEdge* orig) {
    HalfEdge* cross = orig->next;

    //// std::cout << "\nTesting edge " << cross << " with target "
    //           << cross->target->position.X << ", "
    //           << cross->target->position.Y << ", "
    //           << cross->target->position.Z << " and flip target "
    //           << cross->flip->target->position.X << ", "
    //           << cross->flip->target->position.Y << ", "
    //           << cross->flip->target->position.Z 
    //           << " and target side " << planeSide(cross->target)
    //           << " and flip target side " << planeSide(cross->flip->target) << std::endl;

    //// std::cout << "The next edge of " << cross << " has target "
    //           << cross->next->target->position.X << ", "
    //           << cross->next->target->position.Y << ", "
    //           << cross->next->target->position.Z << " and flip target "
    //           << cross->next->flip->target->position.X << ", "
    //           << cross->next->flip->target->position.Y << ", "
    //           << cross->next->flip->target->position.Z << std::endl;

    // was cross->next->target
    while (planeSide(cross->target) == outside){
        cross = cross->next;

        //// std::cout << "\nTesting edge " << cross << " with target "
        //       << cross->target->position.X << ", "
        //       << cross->target->position.Y << ", "
        //       << cross->target->position.Z << " and flip target "
        //       << cross->flip->target->position.X << ", "
        //       << cross->flip->target->position.Y << ", "
        //       << cross->flip->target->position.Z 
        //       << " and target side " << planeSide(cross->target) << " (dist " << planeDist(cross->target)
        //       << ") and flip target side " << planeSide(cross->flip->target) << " (dist "  << planeDist(cross->flip->target) << ")" << std::endl;
    }
    // was cross
    return cross->flip;
}

// Returns a Vector3 corresponding to the point at which the cutting plane
// intersects a given HalfEdge.
Vector3 voronoiCell::planeEdgeIntersect(HalfEdge* edge) {
    Vector3 p1 = edge->target->position;
    Vector3 p2 = edge->flip->target->position;
    Vector3 intersectionPoint(0, 0, 0);
    Vector3 normal = neighborParticle.position - particle.position;
    Vector3 edgeVec = p2 - p1;
    intersectionPoint = p1 + edgeVec * (((particle.position + neighborParticle.position)/2
                        - p1).Dot(normal)) / (edgeVec.Dot(normal));

    return intersectionPoint;
}

// Creates an edge pair between a given pair of vertices
HalfEdge* voronoiCell::addEdgePair(Vertex* vertex1, Vertex* vertex2) {
    //// std::cout << "Adding edge pair between points " << vertex1->position.X << ", " << vertex1->position.Y << ", " << vertex1->position.Z
    //           << " and " << vertex2->position.X << ", " << vertex2->position.Y << ", " << vertex2->position.Z << std::endl;
    HalfEdge* forwardEdge = new HalfEdge(vertex2);
    ////std::cout << "Allocating halfedge: addEdgePair 1" << std::endl;
    HalfEdge* backEdge = new HalfEdge(vertex1);
    ////std::cout << "Allocating halfedge: addEdgePair 2" << std::endl;
    forwardEdge->flip = backEdge;
    backEdge->flip = forwardEdge;
    return forwardEdge;
}

HalfEdge* voronoiCell::maintainFirstEdge(HalfEdge* edge) {
    //// std::cout << "Maintaining firstEdge" <<std::endl;
    std::vector<HalfEdge*> testEdges;
    testEdges.push_back(edge);
    std::size_t i = 0;
    while (planeSide(edge->flip->target) != inside) {
        if (!edge->seen) {
            //// std::cout << "edge not seen" <<std::endl;
            testEdges.push_back(edge->next);
            testEdges.push_back(edge->flip);
            edge->seen = true;
        }
        //// std::cout << "Looking for firstEdge " << edge << " with target "
        //   << edge->target->position.X << ", "
        //   << edge->target->position.Y << ", "
        //   << edge->target->position.Z << " and flip target "
        //   << edge->flip->target->position.X << ", "
        //   << edge->flip->target->position.Y << ", "
        //   << edge->flip->target->position.Z                   
        //   << " and target side " << planeSide(edge->target)
        //   << " and flip target side " << planeSide(edge->flip->target) << std::endl;
        i++;
        edge = testEdges[i];
    }
    //// std::cout << "Found firstEdge " << edge << " with target "
    //           << edge->target->position.X << ", "
    //           << edge->target->position.Y << ", "
    //           << edge->target->position.Z << " and flip target "
    //           << edge->flip->target->position.X << ", "
    //           << edge->flip->target->position.Y << ", "
    //           << edge->flip->target->position.Z                   
    //           << " and target side " << planeSide(edge->target) << " (dist " << planeDist(edge->target)
    //           << ") and flip target side " << planeSide(edge->flip->target) << " (dist "  << planeDist(edge->flip->target) << ")" << std::endl;
    //// std::cout << "We are maintaining firstEdge to an edge with creator " << edge->creator->id <<std::endl;
    return edge;

    // if (planeSide(edge->flip->target) != outside) {
        //// std::cout << "Returning firstEdge" <<std::endl;
    //     return edge;
    // } 
    // if (!edge->seen) {
        //// std::cout << "Not returning firstEdge" <<std::endl;
    //     edge->seen = true;
    //     maintainFirstEdge(edge->next);
    //     maintainFirstEdge(edge->flip);
    // }
}

void voronoiCell::cutCell(const Particle& neighbor) {

    // DEBUGGING OUTPUT TO CHECK WHAT THE VERTICES ARE AT THE BEGINNING OF EACH
    // PLANE CUT

    // std::vector<double> test_v;
    // vertices(test_v);
    
    //// std::cout << "We obtained some vertices." << std::endl;
    // for (size_t i = 0; i < test_v.size(); ++i) {
    //     if (i % 3 == 0) {
            //// std::cout << "Vertex " << (i / 3) << ":" << std::endl;
    //     }
        //// std::cout << test_v[i] << std::endl;
    // }

    // END DEBUG

    // std::cout << "Cutting cell with point " << neighbor.position.X << ", " << neighbor.position.Y << ", " << neighbor.position.Z << std::endl;
    neighborParticle = neighbor;
    //// std::cout << "neighbor cut: " << neighbor.X << " " << neighbor.Y << " " << neighbor.Z << std::endl;

    // We need to maintain firstEdge so that we can index into the Voronoi cell.
    // Specifically, we do not want its flip's target to be outside the plane,
    // or it will be deleted during the resulting cut. Therefore, we have the following
    // simple recursive routine for doing this:
    firstEdge = maintainFirstEdge(firstEdge);
    //// std::cout << "First edge creator id changed to: " << firstEdge->creator->id <<std::endl;
    resetEdges(firstEdge);

    HalfEdge* startingIncidentEdge;

    // Determine whether there is an edge incident to the plane; if so, find
    // it (startingIncidentEdge is passed by reference)
    bool cuttable = findSomeIncidentEdge(startingIncidentEdge);

    if (!cuttable) {
        //// std::cout << "cutCell returning..." << std::endl;    
        return;
    }

    HalfEdge* nextIncidentEdge = startingIncidentEdge;
    //// std::cout << "Found some incident edge " << startingIncidentEdge << " with target " 
    //           << startingIncidentEdge->target->position.X << ", "
    //           << startingIncidentEdge->target->position.Y << ", "
    //           << startingIncidentEdge->target->position.Z << " and flip target "
    //           << startingIncidentEdge->flip->target->position.X << ", "
    //           << startingIncidentEdge->flip->target->position.Y << ", "
    //           << startingIncidentEdge->flip->target->position.Z << std::endl;

    std::vector<HalfEdge*> outsideEdges;
    HalfEdge* prevIncidentEdge = nextIncidentEdge;
    HalfEdge* lastInsideEdge;

    do {

        // was inside
        if (planeSide(nextIncidentEdge->target) == outside) {      
            // Create a new vertex at the location of the intersection of
            // the crossing edge and the cutting plane.
            Vector3 vertexLoc = planeEdgeIntersect(nextIncidentEdge);
            Vertex* newVertex = new Vertex(vertexLoc);
            ////std::cout << "Allocating vertex: cutCell" << std::endl;

            //// std::cout << "\nMaking a new vertex at "
            //           << newVertex->position.X << ", "
            //           << newVertex->position.Y << ", "
            //           << newVertex->position.Z << std::endl;
            
            // was "splitCrossingEdge", but that was not defined -- should it
            // have been splitEdge?
            splitEdge(nextIncidentEdge, newVertex);

            // The first time, make sure to update the first incident edge
            if (prevIncidentEdge == nextIncidentEdge) {
                // For later creation of the final edge of the new face
                lastInsideEdge = startingIncidentEdge->next->flip->next;
                startingIncidentEdge = startingIncidentEdge->next;
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
            HalfEdge* prevInsideEdge = prevIncidentEdge;
            HalfEdge* nextInsideEdge = nextIncidentEdge->flip;

            // addEdgePair is now implemented, and thus is no longer a FIXME.
            HalfEdge* newInsideEdge = addEdgePair(prevIncidentEdge->target, nextIncidentEdge->target);
            prevInsideEdge->next = newInsideEdge;
            newInsideEdge->next = nextInsideEdge;
            newInsideEdge->creator = newInsideEdge->next->creator;
            //// std::cout << "Setting inside edge creator to: " << newInsideEdge->creator->id << std::endl;

            //// std::cout << "Found prevInsideEdge " << prevInsideEdge << " with target " 
            //           << prevInsideEdge->target->position.X << ", "
            //           << prevInsideEdge->target->position.Y << ", "
            //           << prevInsideEdge->target->position.Z << " and flip target "
            //           << prevInsideEdge->flip->target->position.X << ", "
            //           << prevInsideEdge->flip->target->position.Y << ", "
            //           << prevInsideEdge->flip->target->position.Z << std::endl;

            //// std::cout << "Found nextInsideEdge " << nextInsideEdge << " with target " 
            //           << nextInsideEdge->target->position.X << ", "
            //           << nextInsideEdge->target->position.Y << ", "
            //           << nextInsideEdge->target->position.Z << " and flip target "
            //           << nextInsideEdge->flip->target->position.X << ", "
            //           << nextInsideEdge->flip->target->position.Y << ", "
            //           << nextInsideEdge->flip->target->position.Z << std::endl;

            //// std::cout << "Found newInsideEdge " << newInsideEdge << " with target " 
            //           << newInsideEdge->target->position.X << ", "
            //           << newInsideEdge->target->position.Y << ", "
            //           << newInsideEdge->target->position.Z << " and flip target "
            //           << newInsideEdge->flip->target->position.X << ", "
            //           << newInsideEdge->flip->target->position.Y << ", "
                      // << newInsideEdge->flip->target->position.Z << std::endl;

            outsideEdges.push_back(newInsideEdge->flip);
        }

        prevIncidentEdge = nextIncidentEdge;

        nextIncidentEdge = findNextIncidentEdge(prevIncidentEdge);
        //// std::cout << "Found the next incident edge " << nextIncidentEdge << " with target " 
        //           << nextIncidentEdge->target->position.X << ", "
        //           << nextIncidentEdge->target->position.Y << ", "
        //           << nextIncidentEdge->target->position.Z << " and flip target "
        //           << nextIncidentEdge->flip->target->position.X << ", "
        //           << nextIncidentEdge->flip->target->position.Y << ", "
        //           << nextIncidentEdge->flip->target->position.Z << std::endl;

        //// std::cout << "The next edge of " << nextIncidentEdge << " has target "
        //       << nextIncidentEdge->next->target->position.X << ", "
        //       << nextIncidentEdge->next->target->position.Y << ", "
        //       << nextIncidentEdge->next->target->position.Z << " and flip target "
        //       << nextIncidentEdge->next->flip->target->position.X << ", "
        //       << nextIncidentEdge->next->flip->target->position.Y << ", "
        //       << nextIncidentEdge->next->flip->target->position.Z << std::endl;

        //// std::cout << "Find another edge if next incident edge is not " << startingIncidentEdge << " with target "
        //       << startingIncidentEdge->target->position.X << ", "
        //       << startingIncidentEdge->target->position.Y << ", "
        //       << startingIncidentEdge->target->position.Z << " and flip target "
        //       << startingIncidentEdge->flip->target->position.X << ", "
        //       << startingIncidentEdge->flip->target->position.Y << ", "
        //       << startingIncidentEdge->flip->target->position.Z << std::endl;

        //// std::cout << "AND (bandaid) nextIncidentEdge is also not " << startingIncidentEdge->flip->next->flip << " with target "
        //       << startingIncidentEdge->flip->next->flip->target->position.X << ", "
        //       << startingIncidentEdge->flip->next->flip->target->position.Y << ", "
        //       << startingIncidentEdge->flip->next->flip->target->position.Z << " and flip target "
        //       << startingIncidentEdge->flip->next->flip->flip->target->position.X << ", "
        //       << startingIncidentEdge->flip->next->flip->flip->target->position.Y << ", "
        //       << startingIncidentEdge->flip->next->flip->flip->target->position.Z << " is also not startingIncidentEdge" << std::endl;
    // } while(nextIncidentEdge != startingIncidentEdge->next);

    // Continue to find incident edges until we have traversed the entire
    // new face, created by the intersection of the plane and cell
                                                     // THIS IS A BANDAID to cover up the lack of tolerance
                                                     // tolerance in the code as currently structured
    } while(nextIncidentEdge != startingIncidentEdge && nextIncidentEdge != startingIncidentEdge->flip->next->flip);

    // We have to manually add the last edge in the newly-created face.
    HalfEdge* prevInsideEdge = prevIncidentEdge;
    HalfEdge* nextInsideEdge = lastInsideEdge;

    // addEdgePair is now implemented, and thus is no longer a FIXME.
    HalfEdge* newInsideEdge = addEdgePair(prevIncidentEdge->target, lastInsideEdge->flip->target);
    prevInsideEdge->next = newInsideEdge;
    //// std::cout << "Setting next of insideEdge " << prevInsideEdge
    //     << " to " << newInsideEdge << std::endl;
    newInsideEdge->next = nextInsideEdge;
    //// std::cout << "Setting next of insideEdge " << newInsideEdge
    //     << " to " << nextInsideEdge << std::endl;

    // All edges on a face should be associated with the
    // same neighbor particle.
    newInsideEdge->creator = newInsideEdge->next->creator;
    //// std::cout << "Setting last inside edge creator to: " << newInsideEdge->creator->id << std::endl;

    //// std::cout << "Found prevInsideEdge " << prevInsideEdge << " with target " 
    //           << prevInsideEdge->target->position.X << ", "
    //           << prevInsideEdge->target->position.Y << ", "
    //           << prevInsideEdge->target->position.Z << " and flip target "
    //           << prevInsideEdge->flip->target->position.X << ", "
    //           << prevInsideEdge->flip->target->position.Y << ", "
    //           << prevInsideEdge->flip->target->position.Z << std::endl;

    //// std::cout << "Found nextInsideEdge " << nextInsideEdge << " with target " 
    //           << nextInsideEdge->target->position.X << ", "
    //           << nextInsideEdge->target->position.Y << ", "
    //           << nextInsideEdge->target->position.Z << " and flip target "
    //           << nextInsideEdge->flip->target->position.X << ", "
    //           << nextInsideEdge->flip->target->position.Y << ", "
    //           << nextInsideEdge->flip->target->position.Z << std::endl;

    //// std::cout << "Found newInsideEdge " << newInsideEdge << " with target " 
    //           << newInsideEdge->target->position.X << ", "
    //           << newInsideEdge->target->position.Y << ", "
    //           << newInsideEdge->target->position.Z << " and flip target "
    //           << newInsideEdge->flip->target->position.X << ", "
    //           << newInsideEdge->flip->target->position.Y << ", "
    //           << newInsideEdge->flip->target->position.Z << std::endl;

    outsideEdges.push_back(newInsideEdge->flip);

    // Set the next of each edge in the new face, and set
    // the creator of each edge in the new face to the particle
    // we are cutting by.
    for (std::size_t i = 1; i < outsideEdges.size(); ++i) {
        outsideEdges[i]->next = outsideEdges[i-1];
        outsideEdges[i]->creator = &neighbor;
        //// std::cout << "Setting outisde edge creator to: " << outsideEdges[i]->creator->id << std::endl;
        //// std::cout << "Setting next of outsideEdge " << outsideEdges[i]
        // << " to " << outsideEdges[i-1] << std::endl;
    }
    outsideEdges.front()->next = outsideEdges.back();
    outsideEdges.front()->creator = &neighbor;
    //// std::cout << "Setting first outisde edge creator to: " << outsideEdges.front()->creator->id << std::endl;
    //// std::cout << "Setting next of first outsideEdge " << outsideEdges.front()
    // << " to " << outsideEdges.back() << std::endl;

    // cleanUp(outsideEdges.front());
    // cleanUp(startingIncidentEdge->flip);

    // DISREGARD CLEANING UP (NOTE THAT THIS IS QUITE BAD AND MAKES A MEMORY LEAK)
    cleanUp(startingIncidentEdge);
}

void voronoiCell::cleanUp(HalfEdge* edge){
    //// std::cout << "Cleaning up..." << std::endl;
    std::stack<HalfEdge*>* deleteStackEdge = new std::stack<HalfEdge*>;
    ////std::cout << "Allocating stack of halfedges: cleanup" << std::endl;
    std::stack<Vertex*>* deleteStackVertex = new std::stack<Vertex*>;
    ////std::cout << "Allocating stack of vertices: cleanup" << std::endl;
    deleteSearch(deleteStackEdge, deleteStackVertex, edge);
    //// std::cout << "Found everything to delete" << std::endl;
    while (!deleteStackEdge->empty()) {
        HalfEdge* topEdge = deleteStackEdge->top();
        deleteStackEdge->pop();
        //// std::cout << "Deleting an edge..." << std::endl;
        delete topEdge;
        ////std::cout << "Deleting halfedge: cleanup" << std::endl;
    }
    while (!deleteStackVertex->empty()) {
        Vertex* topVertex = deleteStackVertex->top();
        //// std::cout << "Deleting vertex: "<<topVertex->position.X << " " 
        // << topVertex->position.Y << " " << topVertex->position.Z << " " << planeSide(topVertex) << " "
        // << planeDist(topVertex)<< std::endl; 
        deleteStackVertex->pop();
        //// std::cout << "Deleting a vertex..." << std::endl;
        delete topVertex;
        ////std::cout << "Deleting vertex: cleanup" << std::endl;
    }
    delete deleteStackVertex;
    ////std::cout << "Deleting stack of vertices: cleanup" << std::endl;
    delete deleteStackEdge;
    ////std::cout << "Deleting stack of halfedges: cleanup" << std::endl;
}

void voronoiCell::deleteSearch(std::stack<HalfEdge*>* deleteStackEdge,
                  std::stack<Vertex*>* deleteStackVertex,
                  HalfEdge* edge){
    //// std::cout << "Searching for things to delete..." << std::endl;

        // Currently needs to watch out for problems with edges
        // within the cutting plane
    if ((planeSide(edge->target) == outside || planeSide(edge->flip->target) == outside)
	  && !edge->deleteFlag) {
             //// std::cout << "Found an edge to delete: " << nextEdge << " with target " 
             //      << nextEdge->target->position.X << ", "
             //      << nextEdge->target->position.Y << ", "
             //      << nextEdge->target->position.Z << " and flip target "
             //      << nextEdge->flip->target->position.X << ", "
             //      << nextEdge->flip->target->position.Y << ", "
             //      << nextEdge->flip->target->position.Z << std::endl;
        deleteStackEdge->push(edge);
        edge->deleteFlag = true;
        deleteSearch(deleteStackEdge, deleteStackVertex, edge->next);
        deleteSearch(deleteStackEdge, deleteStackVertex, edge->flip);    
    }

    if (planeSide(edge->target) == outside && !edge->target->deleteFlag) {
        //// std::cout << "Found a vertex to delete: " 
        //           << nextEdge->target->position.X << ", "
        //           << nextEdge->target->position.Y << ", "
        //           << nextEdge->target->position.Z << std::endl;
        deleteStackVertex->push(edge->target);
        edge->target->deleteFlag = true;
    }

    // added at the last minute, might work?
    //if (planeSide(nextEdge->flip->target) == outside && !nextEdge->flip->target->deleteFlag) {
        //// std::cout << "Found a vertex to delete: " 
        //           << nextEdge->flip->target->position.X << ", "
        //           << nextEdge->flip->target->position.Y << ", "
        //           << nextEdge->flip->target->position.Z << std::endl;
    //    deleteStackVertex->push(nextEdge->flip->target);
    //    nextEdge->flip->target->deleteFlag = true;
    //}
}

// Calculates the volume of the current Voronoi cell, although note that
// a better implementation of this is possible if we store face_vertices in a
// persistent manner.
double voronoiCell::volume() {
   // std::cout << "Finding a volume" << std::endl;
    double volume = 0;

    // First, we need to get all the face vertices (if they haven't already been
    // calculated)
    if (faceVertices.size() == 0) {
        reset(firstEdge);
        // std::cout << firstEdge->getPrev()->seen << std::endl;
        getFaceVertex(firstEdge);
        resetEdgesAndVertices(firstEdge);
        // std::cout << "Found face vertex" << std::endl;
        resetEdges(firstEdge);
        // std::cout << "reset edges" << std::endl;
    }

    // We then calculate the volume through a tetrahedral decomposition.
    // First, we select one vertex, which will serve as one corner of all
    // of the tetrahedra.
    Vertex* firstVertex = firstEdge->target;

    // Next, loop over the face vertices, decomposing each face into
    // triangles, which can then be used to decompose the cell into tetrahedra.
    for (std::size_t i = 0; i < faceVertices.size(); ++i) {
       // std::cout << "Starting loop: volume " << faceVertices[i]->edges.size() << std::endl;
        HalfEdge* firstBaseEdge = faceVertices[i]->edges[0];
        HalfEdge* secondBaseEdge = faceVertices[i]->edges[1];
        HalfEdge* thirdBaseEdge = faceVertices[i]->edges[2];

        while (thirdBaseEdge != firstBaseEdge) {
           // std::cout << "About to find tet volume" << std::endl;
            volume += tetVolume(firstVertex, firstBaseEdge->target,
                                secondBaseEdge->target, thirdBaseEdge->target);
           // std::cout << "Found tet volume" << std::endl;
            secondBaseEdge = thirdBaseEdge;
            thirdBaseEdge = thirdBaseEdge->next;
        }
    }
    //// std::cout << volume << std::endl;

    return volume;
}

// Calculates the volume of a tetrahedron composed of four input vertices.
double voronoiCell::tetVolume(Vertex* vertex1, Vertex* vertex2,
                              Vertex* vertex3, Vertex* vertex4) {
    Vector3 firstVector = vertex2->position - vertex1->position;
   // std::cout << "Found a vector with X = " << firstVector.X
              // << " and Y = " << firstVector.Y
              // << " and Z = " << firstVector.Z << std::endl;
    Vector3 secondVector = vertex3->position - vertex1->position;
   // std::cout << "Found a vector with X = " << secondVector.X
              // << " and Y = " << secondVector.Y
              // << " and Z = " << secondVector.Z << std::endl;
    Vector3 thirdVector = vertex4->position - vertex1->position;
   // std::cout << "Found a vector with X = " << thirdVector.X
              // << " and Y = " << thirdVector.Y
              // << " and Z = " << thirdVector.Z << std::endl;

    // Compute the vector triple product.
    //// std::cout << "Found a tetrahedron with volume " 
    //           << std::abs(firstVector.Dot(secondVector.Cross(thirdVector))) / 6 << std::endl;
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

    // std::cout << faceVertices.size() << " " << face_v.size() << std::endl;

    for (std::size_t i = 0; i < faceVertices.size(); ++i) {
        // We defaulted the creator to particle during initialization,
        // so this is a check for if the "neigbor" is the boundary
        //// std::cout << "Pushing a neighbor with id " 
        //           << faceVertices[i]->edges[0]->creator->id << std::endl;
        if (faceVertices[i]->edges[0]->creator->id != particle.id) {
            v.push_back(faceVertices[i]->edges[0]->creator->id);
        }
    }
   // std::cout << "Done pushing neighbors:" << v.size() << std::endl;
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

        // std::cout << "Starting loop: volume " << faceVertices[i]->edges.size() << std::endl;
        HalfEdge* firstBaseEdge = faceVertices[i]->edges[0];
        HalfEdge* secondBaseEdge = faceVertices[i]->edges[1];
        HalfEdge* thirdBaseEdge = faceVertices[i]->edges[2];

        while (thirdBaseEdge != firstBaseEdge) {
            // std::cout << "About to find tet volume" << std::endl;
            area += triArea(firstBaseEdge->target,
                            secondBaseEdge->target, thirdBaseEdge->target);
            // std::cout << "Found tet volume" << std::endl;
            secondBaseEdge = thirdBaseEdge;
            thirdBaseEdge = thirdBaseEdge->next;
        }

        v.push_back(area);
    }
}

// Calculates the area of a triangle composed of three input vertices.
double voronoiCell::triArea(Vertex* vertex1, Vertex* vertex2,
                              Vertex* vertex3) {
    Vector3 firstVector = vertex2->position - vertex1->position;
   // std::cout << "Found a vector with X = " << firstVector.X
              // << " and Y = " << firstVector.Y
              // << " and Z = " << firstVector.Z << std::endl;
    Vector3 secondVector = vertex3->position - vertex1->position;
   // std::cout << "Found a vector with X = " << secondVector.X
              // << " and Y = " << secondVector.Y
              // << " and Z = " << secondVector.Z << std::endl;

    return (std::abs(firstVector.Cross(secondVector).distanceTo(Vector3(0, 0, 0)))/ 2);
}

void voronoiCell::face_vertices(std::vector<int> &v) {
	//// std::cout << "First edge creator id: " << firstEdge->creator->id <<std::endl;
    //// std::cout << "face_vertices called" << std::endl;
    resetEdgesAndVertices(firstEdge);
    getFaceVertex(firstEdge);
    //// std::cout << "faceVertices has a size, which is" << std::endl;
    //// std::cout << faceVertices.size() << std::endl;
    for (std::size_t i = 0; i < faceVertices.size(); ++i) {
        // THIS IS NOT THE FORMAT THIS SHOULD BE OUTPUT IN, since unfortunately
        // the voro++ output relies on a well-defined vertex ordering. We need to
        // make a decision over whether or not we will permanently store the
        // vertices, faceVertices after first calculating them.
        v.push_back(i);
        //// std::cout << "Face vertex neighbor id: " << faceVertices[i]->edges[0]->creator->id << std::endl;

        // \/\/\/ WE USED TO DELETE EVERYTHING BUT NOW WE STORE THEM
        // Right now, since the face vertices are not stored in any persistent
        // manner, we need to deallocate them to avoid memory leaks.
        // delete faceVertices[i];
    }

    // We eventually need to reset all the seen flags at the end of this. or,
    // more likely, we need to optimize this procedure.
    resetEdges(firstEdge);
};

void voronoiCell::getFaceVertex(HalfEdge* testEdge) {
    //// std::cout << "getFaceVertex " << testEdge << ", with vertex " << testEdge->target << ": " 
    //               << testEdge->target->position.X << ", "
    //               << testEdge->target->position.Y << ", "
    //               << testEdge->target->position.Z << std::endl;
    if (!testEdge->seen) {
       // std::cout << "First Edge target: " << testEdge->target->position.X << " " << testEdge->target->position.Y << " " << testEdge->target->position.Z << " " << std::endl;
        FaceVertex* newFaceVertex = new FaceVertex;
        ////std::cout << "Allocating facevertex: getFaceVertex" << std::endl;
        newFaceVertex->edges.push_back(testEdge);
		//// std::cout << "Face vertex neighbor id asdf: " << testEdge->creator->id << std::endl;

        HalfEdge* otherEdge = testEdge->next;

        testEdge->seen = true;
       // std::cout<<testEdge->getPrev()->seen<<std::endl;
        //// std::cout << "seen" << std::endl;

        while (otherEdge != testEdge) {
            
            newFaceVertex->edges.push_back(otherEdge);
            otherEdge->seen = true;
                //// std::cout << "seen" << std::endl;
           // std::cout << "Other Edge target: " << otherEdge->target->position.X << " " << otherEdge->target->position.Y << " " << otherEdge->target->position.Z << " " << std::endl;
            otherEdge = otherEdge->next;
           // std::cout << (otherEdge != testEdge) << otherEdge->seen << std::endl;

        }
       // std::cout << "Face Vertex edges: " << newFaceVertex->edges.size() << std::endl;
        for (size_t i = 0; i < newFaceVertex->edges.size(); ++i){
           // std::cout << "Face Vertex edge target: " << newFaceVertex->edges[i]->target->position.X << " " << newFaceVertex->edges[i]->target->position.Y << " " << newFaceVertex->edges[i]->target->position.Z << " " << std::endl;
        }

        faceVertices.push_back(newFaceVertex);
        //// std::cout << "Face vertex neighbor id asdfasdf: " << newFaceVertex->edges[0]->creator->id << std::endl;

        getFaceVertex(testEdge->flip);

        otherEdge = testEdge->next;
        while (otherEdge != testEdge) {
            getFaceVertex(otherEdge->flip);
            otherEdge = otherEdge->next;
        }
    }
};

// CustomInterface IO
void voronoiCell::vertices(std::vector<double> &v) {
    std::vector<Vertex*> vertices;
    getVertex(firstEdge, vertices);
    resetEdgesAndVertices(firstEdge);

    for (std::size_t i = 0; i < vertices.size(); ++i) {
        v.push_back(vertices[i]->position.X);
        v.push_back(vertices[i]->position.Y);
        v.push_back(vertices[i]->position.Z);
    }

    // we eventually need to reset all the seen flags at the end of this. or,
    // more likely, we need to optimize this procedure.
    resetEdgesAndVertices(firstEdge);
};

void voronoiCell::getVertex(HalfEdge* testEdge, std::vector<Vertex*> &vertices) {
    //// std::cout << "getVertex " << testEdge << ", with vertex " << testEdge->target << ": " 
    //               << testEdge->target->position.X << ", "
    //               << testEdge->target->position.Y << ", "
    //               << testEdge->target->position.Z << std::endl;
    if (!testEdge->seen) {
        testEdge->seen = true;
        //// std::cout << "seen" << std::endl;
        
        if (!testEdge->target->seen) {
            testEdge->target->seen = true;
            vertices.push_back(testEdge->target);
        }
        
        //// std::cout << "pushback" << std::endl;
        //// std::cout << "flip of " << testEdge << std::endl;
        getVertex(testEdge->flip, vertices);
        //// std::cout << "next of " << testEdge << std::endl;
        getVertex(testEdge->next, vertices);
    }
};

// Resets the seen flags on all edges
void voronoiCell::resetEdges(HalfEdge* edge) {
    if (edge->seen) {
        //// std::cout << "Reset an edge." << std::endl;
        edge->seen = false;
        resetEdges(edge->flip);
        resetEdges(edge->next);
    }
};

// Resets the seen flags on all edges and vertices
void voronoiCell::resetEdgesAndVertices(HalfEdge* edge) {
    if (edge->seen) {
        //// std::cout << "Reset an edge (looking for vertices)." << std::endl;
        edge->seen = false;

        // Reset the associated vertex, if we need to
        if (edge->target->seen) {
            //// std::cout << "Reset a vertex." << std::endl;
            edge->target->seen = false;
        }

        resetEdgesAndVertices(edge->flip);
        resetEdgesAndVertices(edge->next);
    }
};

void voronoiCell::reset(HalfEdge* edge){
    //// std::cout << "Cleaning up..." << std::endl;
    std::stack<HalfEdge*>* seenStackEdge = new std::stack<HalfEdge*>;
    ////std::cout << "Allocating stack of halfedges: cleanup" << std::endl;
    std::stack<Vertex*>* seenStackVertex = new std::stack<Vertex*>;
    ////std::cout << "Allocating stack of vertices: cleanup" << std::endl;
    seenSearch(seenStackEdge, seenStackVertex, edge);
    //// std::cout << "Found everything to seen" << std::endl;
    while (!seenStackEdge->empty()) {
        HalfEdge* topEdge = seenStackEdge->top();
        seenStackEdge->pop();
        //// std::cout << "Deleting an edge..." << std::endl;
        topEdge->seen = false;
        ////std::cout << "Deleting halfedge: cleanup" << std::endl;
    }
    while (!seenStackVertex->empty()) {
        Vertex* topVertex = seenStackVertex->top();
        //// std::cout << "Deleting vertex: "<<topVertex->position.X << " " 
        // << topVertex->position.Y << " " << topVertex->position.Z << " " << planeSide(topVertex) << " "
        // << planeDist(topVertex)<< std::endl; 
        seenStackVertex->pop();
        //// std::cout << "Deleting a vertex..." << std::endl;
        topVertex->seen = false;
        ////std::cout << "Deleting vertex: cleanup" << std::endl;
    }
    delete seenStackVertex;
    ////std::cout << "Deleting stack of vertices: cleanup" << std::endl;
    delete seenStackEdge;
    ////std::cout << "Deleting stack of halfedges: cleanup" << std::endl;
};

void voronoiCell::seenSearch(std::stack<HalfEdge*>* seenStackEdge,
                  std::stack<Vertex*>* seenStackVertex,
                  HalfEdge* edge){
    //// std::cout << "Searching for things to seen..." << std::endl;

        // Currently needs to watch out for problems with edges
        // within the cutting plane
    if (!edge->seen) {
             //// std::cout << "Found an edge to seen: " << nextEdge << " with target " 
             //      << nextEdge->target->position.X << ", "
             //      << nextEdge->target->position.Y << ", "
             //      << nextEdge->target->position.Z << " and flip target "
             //      << nextEdge->flip->target->position.X << ", "
             //      << nextEdge->flip->target->position.Y << ", "
             //      << nextEdge->flip->target->position.Z << std::endl;
        seenStackEdge->push(edge);
        edge->seen = true;
        seenSearch(seenStackEdge, seenStackVertex, edge->next);
        seenSearch(seenStackEdge, seenStackVertex, edge->flip);    
    }

    if (!edge->target->seen) {
        //// std::cout << "Found a vertex to seen: " 
        //           << nextEdge->target->position.X << ", "
        //           << nextEdge->target->position.Y << ", "
        //           << nextEdge->target->position.Z << std::endl;
        seenStackVertex->push(edge->target);
        edge->target->seen = true;
    }

    // added at the last minute, might work?
    //if (planeSide(nextEdge->flip->target) == outside && !nextEdge->flip->target->seen) {
        //// std::cout << "Found a vertex to seen: " 
        //           << nextEdge->flip->target->position.X << ", "
        //           << nextEdge->flip->target->position.Y << ", "
        //           << nextEdge->flip->target->position.Z << std::endl;
    //    seenStackVertex->push(nextEdge->flip->target);
    //    nextEdge->flip->target->seen = true;
    //}
};


// WORRY ABOUT THIS LATER

// std::vector<FaceVertex>* voronoiCell::getFaceVertices() {
//     faceVertices = new std::vector<FaceVertex>;
    
//     faceVertices->push(new FaceVertex);
//     firstEdge->face = faceVertices[-1];
//     faceVertices[-1]->edges->push(firstEdge);

//     HalfEdge* nextEdge = firstEdge->next;
//     while (nextEdge != firstEdge) {
//         nextEdge->face = faceVertices[-1];
//         faceVertices[-1]->edges->push(nextEdge);
//     }
// }

/** Outputs the edges of the Voronoi cell in gnuplot format to an output stream.
 * \param[in] (x,y,z) a displacement vector to be added to the cell's position.
 * \param[in] fp a file handle to write to. */
// Shamelessly stolen from voro++
// Remember to resetEdges after calling this
void voronoiCell::drawGnuplot(double dispX, double dispY, double dispZ, FILE* fp, HalfEdge* edge) {
    HalfEdge* currentEdge = edge->next;
    fprintf(fp, "%g, %g, %g\n", currentEdge->flip->target->position.X + dispX,
                                currentEdge->flip->target->position.Y + dispY,
                                currentEdge->flip->target->position.Z + dispZ);

    // Loop around the current face, drawing all edges.
    while (currentEdge != edge) {
        fprintf(fp, "%g, %g, %g\n", currentEdge->target->position.X + dispX,
                                    currentEdge->target->position.Y + dispY,
                                    currentEdge->target->position.Z + dispZ);

        currentEdge = currentEdge->next;
    }

    fprintf(fp, "%g, %g, %g\n", currentEdge->target->position.X + dispX,
                                currentEdge->target->position.Y + dispY,
                                currentEdge->target->position.Z + dispZ);

    fputs("\n\n", fp);

    // Loop around the face again, recursing on faces that have not yet been drawn
    currentEdge->seen = true;
    if (!currentEdge->flip->seen) {
        drawGnuplot(dispX, dispY, dispZ, fp, currentEdge->flip);
    }
    currentEdge = currentEdge->next;

    while (currentEdge != edge) {
        currentEdge->seen = true;
        if (!currentEdge->flip->seen) {
            drawGnuplot(dispX, dispY, dispZ, fp, currentEdge->flip);
        }
        currentEdge = currentEdge->next;
    }

    // int i,j,k,l,m;
    // for(i=1;i<p;i++) for(j=0;j<nu[i];j++) {
    //     k=ed[i][j];
    //     if(k>=0) {
    //         fprintf(fp,"%g %g %g\n",x+0.5*pts[3*i],y+0.5*pts[3*i+1],z+0.5*pts[3*i+2]);
    //         l=i;m=j;
    //         do {
    //             ed[k][ed[l][nu[l]+m]]=-1-l;
    //             ed[l][m]=-1-k;
    //             l=k;
    //             fprintf(fp,"%g %g %g\n",x+0.5*pts[3*k],y+0.5*pts[3*k+1],z+0.5*pts[3*k+2]);
    //         } while (search_edge(l,m,k));
    //         fputs("\n\n",fp);
    //     }
    // }
    // reset_edges();
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

voronoiCell* cellContainer::makeCell(Particle particle) {


    sds.initialize(&particles[0], &particles[particles.size()]);
    // Calculate the initial max radius of a cube
    double maxRadius = sqrt(3) * defaultLength;
//// std::cout<< "Initial particle before cell creation: "<< particle.X << " " << particle.Y << " "<< particle.Z <<std::endl;

    // Initialize the voronoi cell as a cube...
    voronoiCell* cell = new voronoiCell("cube", defaultLength, particle, maxRadius,
        x_min, x_max, y_min, y_max, z_min, z_max);
    ////std::cout << "Allocating voronoiCell: makeCell" << std::endl;
//// std::cout<< "Initial particle: "<< cell->particle.X << " " << cell->particle.Y << " "<< cell->particle.Z <<std::endl;

    // loop through our array of particles and cut the cell with planes
    // associated with each possible neighbor.

    std::array<double, 3> point = {{particle.position.X, particle.position.Y, particle.position.Z}};

    std::vector<Particle*> NeighborParticles;

    sds.neighborQuery(point, defaultLength, &NeighborParticles);

    // std::cout << NeighborParticles.size() << std::endl;
    for (size_t i = 0; i < NeighborParticles.size(); ++i) {
        // This condition may need to be fixed.
        // std::cout << iter->position.X << " " << iter->position.Y << " " << iter->position.Z << std::endl;
        if ((NeighborParticles[i]->position != particle.position) && (particle.position.distanceTo(NeighborParticles[i]->position) < maxRadius)) {
            // std::cout << cell->particle.position.X << " " << cell->particle.position.Y << " "<< cell->particle.position.Z <<std::endl;
            cell->cutCell(*NeighborParticles[i]);
        }
    }
    // std::cout << "Finished making a cell: " << cell->particle.position.X << " " << cell->particle.position.Y << " "<< cell->particle.position.Z  << std::endl;
    return cell;
}

cellContainer::~cellContainer() {
    //// std::cout << "Entering cellContainer destructor" << std::endl;
    while (!cells.empty()) {
        delete cells.back();
        ////std::cout << "Deleting voronoiCell: ~cellContainer" << std::endl;
        cells.pop_back();
    }
}

double cellContainer::sum_cell_volumes() {
    double sum = 0;
    if(!calculated) {
        for(unsigned int i = 0; i < particles.size(); ++i) {
           // std::cout << "Making cell: " << particles[i].position.X << " " << particles[i].position.Y << " " << particles[i].position.Z << std::endl;
            cells.push_back(makeCell(particles[i]));
           // std::cout << "Made cell: " << std::endl;
            sum += cells[i]->volume();
        }
    }
    else {
        //// std::cout << "Don't need to make cells" <<z std::endl;
        for(unsigned int i = 0; i < particles.size(); ++i) {
            sum += cells[i]->volume();
        }
    }
    return sum;
}


void cellContainer::put(int id, double px, double py, double pz) {
    particles.push_back(Particle(id, px, py, pz));
}

double cellContainer::findMaxNeighDist() {
    double maxDist = 0;
    for(voronoiCell* cellPointer : cells) {
        std::vector<int> neighborIndexes;
        cellPointer->neighbors(neighborIndexes);
        for (int i = 0; i < neighborIndexes.size(); ++i) {
            double dist = cellPointer->particle.position.distanceTo(particles[neighborIndexes[i]].position);
            if(dist > maxDist) {
                maxDist = dist;
            }
        }
    }
    return maxDist;

}



// This function returns a random floating point number between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

