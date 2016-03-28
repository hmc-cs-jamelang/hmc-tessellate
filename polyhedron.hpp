#pragma once

#include <vector>
#include <iostream>
#include "structpool.hpp"
#include "verification.hpp"
#include "sds.hpp"

struct Particle;
struct Vector3;
struct ParticleIndex;

class Polyhedron {
protected:
    struct HalfEdge;
    using Vertex = Vector3;
    struct Face;

    StructPool<HalfEdge> edges_;
    StructPool<Vertex> vertices_;
    StructPool<Face> faces_;

    using EdgeIndex = StructPool<HalfEdge>::Index;
    using VertexIndex = StructPool<Vertex>::Index;
    using FaceIndex = StructPool<Face>::Index;

    struct HalfEdge {
        VertexIndex target;
        EdgeIndex flip;
        EdgeIndex next;
        FaceIndex face;

        HalfEdge() = default;
        HalfEdge(VertexIndex vertex) : vertex(vertex) {}
        HalfEdge(VertexIndex vertex, EdgeIndex flip,
                 EdgeIndex next, FaceIndex face)
            : vertex(vertex), flip(flip), next(next), face(face)
        { /* Done */ }

    };

    struct Face {
        EdgeIndex handle;
        ParticleIndex neighbor;
    };

public:


};


class VoronoiCell : Polyhedron {
    Vector position_;
    EdgeIndex handle_;

    enum class PlaneSide {
        OUTSIDE,
        INCIDENT,
        INSIDE
    };

    findAndTrimOutgoingEdge

    EdgeIndex findIngoingEdgeOnSameFace(const Plane& plane, EdgeIndex ei)
    {
        // Declare early so we can use it in the verification
        // expression
        EdgeIndex result VERIFICATION(= INVALID_EDGE);

        VERIFY(outgoing(ei));
        VERIFY_EXIT(ingoing(result));



        return (result = );
    }

    PlaneSide planeSide(const Plane& plane, VertexIndex vi)
    {
        Vector& v = vertices_[vi].position;
        double ans = dot(v, plane.normal) - plane.offset;
        if (std::abs(ans) < TOLERANCE) {
            return Location::INCIDENT;
        }
        else if (ans > 0) {
            return Location::OUTSIDE;
        }
        else {
            return Location::INSIDE;
        }
    }

    Vector intersection(const Plane& plane, EdgeIndex ei)
    {
        Vector& p1 = vector(source(ei));
        Vector& p2 = vector(target(ei));
        Vector edge = p2 - p1;

        // The plane offset as if p1 was the origin.
        auto shiftedOffset = plane.offset - dot(plane.normal, p1);

        // The fraction of the distance from p1 to p2 that we need to
        // move to reach the intersection point.
        auto fractionOfEdge = shiftedOffset / dot(plane.normal, edge);

        // Assert that the edge does indeed intersect the plane
        VERIFY(0 <= fractionOfEdge && fractionOfEdge <= 1);

        return p1 + edge * fractionOfEdge;
    }

    EdgeIndex& next(EdgeIndex ei) {return edges_[ei].next;}
    EdgeIndex& flip(EdgeIndex ei) {return edges_[ei].flip;}
    VertexIndex& target(EdgeIndex ei) {return edges_[ei].target;}

    VERIFICATION(
        bool isOutgoingEdge(const Plane& plane, EdgeIndex ei)
        {

        }
    )

    EdgeIndex findIngoingEdgeOnSameFace(const Plane& plane, EdgeIndex start) {
        VERIFY(isOutgoingEdge(plane, start));
        // Loop around face
        for (EdgeIndex test = next(start); test != start; test = next(test)) {
            if (planeSide(plane, target(test)) == INSIDE) {
                VERIFY(isIngoingEdge(plane, test));
                return test;
            }
        }
        return INVALID_EDGE;
    }

    bool cut(const Particle& neighbor) {
        Plane plane {position, neighbor.position};

        handle_ = findInsideEdge(plane);
        EdgeIndex outgoingEdge = findOutgoingEdge(plane);
        if (outgoingEdge == INVALID_EDGE) {
            return false;
        }

        struct FaceCuttingState {
            EdgeIndex start;
            EdgeIndex current;
            VertexIndex currentVertex;
            VertexIndex prevVertex;
        };

        auto beginCuttingFace(FaceCuttingState& state,
                              EdgeIndex outgoingEdge)
        {
            state.start = state.current = outgoingEdge;
            state.currentVertex = target(outgoingEdge);
        };

        auto nextEdge = [&](FaceCuttingState& state) {
            state.current = next(state.current);
            state.prevVertex = state.currentVertex;
            state.currentVertex = target(state.current);
        };

        auto nextFace = [&](FaceCuttingState& state) {
            beginCuttingFace(state, flip(state.current));
        };

        auto needToCutThisFace = [&](FaceCuttingState& state) -> bool
        {
            do {
                switch (location(state.currentVertex)) {
                    case OUTSIDE: return true;
                    case INSIDE: return false;
                }
                nextEdge(state);
            } while (state.current != state.start);
        };

        auto findIntersection = [&](FaceCuttingState& state) -> VertexIndex
        {
            VERIFY(location(state.currentVertex) == OUTSIDE);

            do {
                nextEdge(state);
            } while (location(state.currentVertex) != INSIDE);

            // If the previous vertex is incident, then we don't want
            // to bother replacing it with a new vertex.
            // If we did, we would start having edges between nearly
            // (or completely) equal vertices if that incident vertex
            // is shared on multiple faces, since each would create their
            // own new version of it.
            // That wouldn't actually break anything in particular AFAIK,
            // but it is undesirable.

            // Since otherwise the loop above would have stopped earlier:
            VERIFY(location(state.prevVertex) != INSIDE);

            if (location(state.prevVertex) == INCIDENT) {
                return state.prevVertex;
            }

            // location is OUTSIDE, since we know it isn't INSIDE.
            else {
                return vertices_.create(
                    intersection(plane,
                                 position(state.prevVertex),
                                 position(state.currentVertex))
                );
            }
        };

        auto cutFace = [&](FaceCuttingState& state, VertexIndex v)
        {
            EdgeIndex newEdge = edges_.create(v);
            next(state.start) = newEdge;
            next(newEdge) = state.current;

            flip(newEdge) = state.outerEdge;
            flip(state.outerEdge) = newEdge;

            EdgeIndex newOuterEdge = edges_.create(v);
            next(newOuterEdge) = state.outerEdge;
            state.outerEdge = newOuterEdge;
        };

        auto finalize = [&](FaceCuttingState& state, EdgeIndex outgoingEdge)
        {
            edges_[state.firstOuterEdge] = edges_[state.outerEdge];
            edges_.destroy(state.outerEdge);
        };

        FaceBuilder newFace = buildFace(neighbor);


        FaceCuttingState state = beginCuttingFace(state, outgoingEdge);
        do {
            if (needToCutThisFace(state)) {
                VertexIndex edgeIntersection = findIntersection(state);
                cutFace(state, edgeIntersection);

                nextFace(state);
                target(state.current) = edgeIntersection;
            }
            else {
                nextFace(state);
            }
        } while (state.start != outgoingEdge);

        // Delete the now-orphaned part of the polyhedron that was cut off
        garbageCollect();
        return true;
    }

    void mark(EdgeIndex ei)
    {
        VERIFY_EXIT(edges_.marked(ei));

        if (edges_.marked(ei)) {return;}
        edges_.setMarked(ei, true);
        mark(flip(ei));
        mark(next(ei));
        vertices_.setMarked(target(ei), true);
    }

    void garbageCollect()
    {
        mark(handle_);
        edges_.destroyUnmarkedAndResetFlags();
        vertices_.destroyUnmarkedAndResetFlags();
    }
};

(v - (c + n)/2) * (n - c)

v * (n - c) - (c + n)/2 * (n - c)

Vector normal = n - c;
Vector halfway = (c + n)/2;
double offset = dot(normal, halfway);

Outgoing means "From INSIDE to OUTSIDE"
It does NOT mean "From INSIDE to INCIDENT"
    no cutting actually has to be done
More crucially it does not mean "From INCIDENT to OUTSIDE"
    then the whole edge should disappear

No wait. What about INSIDE to INCIDENT? That might not be a cut,
yes. So we check the next one, and end up with an INCIDENT to OUTSIDE.
At that point, we say yes indeed it is an outgoing edge. Is that correct?

Almost. So our main problem right now is that we are not totally sure
what to do with incident vertices. Say an incident vertex is shared.
Then each thing should end up pointing to that exact vertex, not some random
vertex nearby. What does that do for us?


OK. So. We have


Ingoing means "From !INSIDE to INSIDE"
It does NOT mean "From OUTSIDE to !OUTSIDE"

keep(vertex) = !OUTSIDE
keep(edge) =
    one vertex INSIDE or both vertices INCIDENT





void cutCell()
{

}



We start on an outgoing edge.
We loop around to find an ingoing edge,
and if we did not find an outside vertex
while doing so, then we are done.

Otherwise, we need to get the intersection as well.


void cut() {

    FaceClipper currentFace;

    do {
        VertexIndex v = findIngoingIntersection();
        if (v != INVALID_VERTEX) {
            clipFace();
        }
        nextFace(v);
    } while (currentFace.outgoingEdge != firstOutgoingEdge);

    bool needToCut = findIngoingEdge();
    if (needToCut) {
        cut(outgoing, ingoing);
    }
    nextFace();
}
