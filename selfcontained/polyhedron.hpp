#pragma once

#include <iostream>
#include <vector>
// #include <queue>
#include <functional>
#include <unordered_set>
#include <iterator>
#include <algorithm>

#include "vectormath.hpp"
#include "structpool.hpp"
#include "verification.hpp"

namespace hmc {
    struct HalfEdge;
    struct Face;

    using EdgePool    = StructPool<HalfEdge>;
    using VertexPool  = StructPool<Vector3>;
    using FacePool    = StructPool<Face>;

    using EdgeIndex   = EdgePool::Index;
    using VertexIndex = VertexPool::Index;
    using FaceIndex   = FacePool::Index;

    static constexpr EdgeIndex INVALID_EDGE     = EdgePool::INVALID_INDEX;
    static constexpr VertexIndex INVALID_VERTEX = VertexPool::INVALID_INDEX;
    static constexpr FaceIndex INVALID_FACE     = FacePool::INVALID_INDEX;

    struct HalfEdge {
        EdgeIndex flip;
        EdgeIndex next;
        VertexIndex target;
        FaceIndex face;

        HalfEdge() = default;
        HalfEdge(FaceIndex f) : flip(), next(), target(), face(f) {}
        HalfEdge(FaceIndex f, VertexIndex v) : flip(), next(), target(v), face(f) {}
        HalfEdge(FaceIndex f, VertexIndex v, EdgeIndex flip, EdgeIndex next)
            : flip(flip), next(next), target(v), face(f)
        { /* Done */ }
    };

    struct Face {
        int id;
        EdgeIndex startingEdge;

        Face(int id)
            : id(id),
              VERIFICATION   ( startingEdge(INVALID_EDGE) )
              NO_VERIFICATION( startingEdge() )
        { /* Done */ }

        Face(int id, EdgeIndex startingEdge) : id(id), startingEdge(startingEdge) {}
    };


    class Polyhedron {
    public:
        struct FaceData {
            FaceIndex face;
            Vector3 weightedNormal;

            FaceData(FaceIndex face, Vector3 weightedNormal)
                : face(face), weightedNormal(weightedNormal)
            { /* Done */ }
        };

        EdgeIndex root_ = INVALID_EDGE;

        EdgePool edges_;
        VertexPool vertices_;
        FacePool faces_;

        std::vector<FaceData> faceData_;

        std::vector<VertexIndex> verticesToDestroy_;
        std::vector<EdgeIndex> edgesToDestroy_;

        VertexIndex maxDistanceVertex_ = INVALID_VERTEX;
        double maximumNeighborDistance_;

        // struct VertexDistance {
        //     VertexIndex vi;
        //     double dist;

        //     VertexDistance(VertexIndex vi, double dist)
        //         : vi(vi), dist(dist)
        //     { /* Done */ }

        //     friend bool operator<(const VertexDistance& a, const VertexDistance& b)
        //     {
        //         return a.dist < b.dist;
        //     }

        //     friend bool operator==(const VertexDistance& a, const VertexDistance& b)
        //     {
        //         return a.vi == b.vi;
        //     }

        //     friend bool operator!=(const VertexDistance& a, const VertexDistance& b)
        //     {
        //         return a.vi != b.vi;
        //     }
        // };

        // struct VertexDistanceHeap {
        //     std::vector<VertexDistance>;



        // } vertexDistances_;

        void clear()
        {
            root_ = INVALID_EDGE;
            edges_.clear();
            vertices_.clear();
            faces_.clear();
            faceData_.clear();
        }

        bool isClear()
        {
            bool isClear = (root_ == INVALID_EDGE);
            if (isClear) {
                VERIFY(edges_.size() == 0);
                VERIFY(vertices_.size() == 0);
                VERIFY(faces_.size() == 0);
                VERIFY(faceData_.size() == 0);
            }
            return isClear;
        }

        void clearComputation()
        {
            faceData_.clear();
        }

        EdgeIndex& flip(EdgeIndex ei) {return edges_[ei].flip;}
        EdgeIndex& next(EdgeIndex ei) {return edges_[ei].next;}
        VertexIndex& target(EdgeIndex ei) {return edges_[ei].target;}
        VertexIndex& source(EdgeIndex ei) {return target(flip(ei));}
        FaceIndex& face(EdgeIndex ei) {return edges_[ei].face;}
        EdgeIndex& startingEdge(FaceIndex fi) {return faces_[fi].startingEdge;}

        const EdgeIndex& flip(EdgeIndex ei) const {return edges_[ei].flip;}
        const EdgeIndex& next(EdgeIndex ei) const {return edges_[ei].next;}
        const VertexIndex& target(EdgeIndex ei) const {return edges_[ei].target;}
        const VertexIndex& source(EdgeIndex ei) const {return target(flip(ei));}
        const FaceIndex& face(EdgeIndex ei) const {return edges_[ei].face;}
        const EdgeIndex& startingEdge(FaceIndex fi) const {return faces_[fi].startingEdge;}

        template <typename... Edge_Constructor_Args>
        EdgeIndex createEdge(Edge_Constructor_Args... args)
        {
            return edges_.create(args...);
        }

        template <typename... Vertex_Constructor_Args>
        VertexIndex createVertex(Vertex_Constructor_Args... args)
        {
            return vertices_.create(args...);
        }

        template <typename... Face_Constructor_Args>
        FaceIndex createFace(Face_Constructor_Args... args)
        {
            return faces_.create(args...);
        }

        void destroy(EdgeIndex ei)
        {
            edges_.destroy(ei);
        }

        void destroy(VertexIndex vi)
        {
            vertices_.destroy(vi);
            if (vi == maxDistanceVertex_) {
                maxDistanceVertex_ = INVALID_VERTEX;
            }
        }

        void destroy(FaceIndex fi)
        {
            faces_.destroy(fi);
        }

        VERIFICATION(
            void verifyIsValidPolyhedron()
            {
                VERIFY(root_ != INVALID_EDGE);
                for (auto eit = edges_.begin(); eit != edges_.end(); ++eit) {
                    VERIFY(eit->target != INVALID_VERTEX);
                    VERIFY(flip(eit->flip) == eit.index());
                    VERIFY(face(eit->next) == eit->face);
                    VERIFY(face(eit->flip) != eit->face);
                    VERIFY(flip(eit->next) != INVALID_EDGE);
                    VERIFY(source(eit->next) == eit->target);
                }
                for (auto fit = faces_.begin(); fit != faces_.end(); ++fit) {
                    VERIFY(face(fit->startingEdge) == fit.index());
                }


                std::unordered_set<EdgeIndex> reachableEdges;
                std::unordered_set<VertexIndex> reachableVertices;
                std::unordered_set<FaceIndex> reachableFaces;

                std::function<void(VertexIndex)> markVertex = [&](VertexIndex vi) -> void {
                    reachableVertices.insert(vi);
                };

                std::function<void(EdgeIndex)> markEdge = [&](EdgeIndex ei) -> void {
                    if (reachableEdges.find(ei) != reachableEdges.end()) { return; }
                    reachableEdges.insert(ei);
                    VERIFY(edges_.active(ei));
                    markEdge(flip(ei));
                    markEdge(next(ei));
                    markVertex(target(ei));
                    reachableFaces.insert(face(ei));
                };


                markEdge(root_);

                for (auto ei : reachableEdges) { VERIFY(edges_.active(ei)); }
                for (auto vi : reachableVertices) { VERIFY(vertices_.active(vi)); }
                for (auto fi : reachableFaces) { VERIFY(faces_.active(fi)); }

                for (auto e = edges_.begin(); e != edges_.end(); ++e) {
                    VERIFY(reachableEdges.find(e.index()) != reachableEdges.end());
                }
                for (auto v = vertices_.begin(); v != vertices_.end(); ++v) {
                    VERIFY(reachableVertices.find(v.index()) != reachableVertices.end());
                }
                for (auto f = faces_.begin(); f != faces_.end(); ++f) {
                    VERIFY(reachableFaces.find(f.index()) != reachableFaces.end());
                }

                // It's easy to get the handedness wrong (-> negative volume).
                VERIFY(temporarilyComputeVolume() > 0);

                if (!isClear()) {
                    // Checking the Euler characteristic, which is probably
                    // unnecessary but might catch SOMETHING odd.
                    auto E = std::distance(edges_.begin(), edges_.end());
                    auto V = std::distance(vertices_.begin(), vertices_.end());
                    auto F = std::distance(faces_.begin(), faces_.end());
                    VERIFY(V - E/2 + F == 2);
                }
            }
        )

        void translate(const Vector3 shift)
        {
            for (auto& v : vertices_) {
                v += shift;
            }
        }

        double computeVolume()
        {
            computeFaceData();
            double volume = 0;
            for (auto& data : faceData_) {
                VERIFY(faces_.active(data.face));
                VERIFY(edges_.active(startingEdge(data.face)));
                VERIFY(vertices_.active(target(startingEdge(data.face))));
                volume += dot(
                    vertices_[target(startingEdge(data.face))],
                    data.weightedNormal
                );
            }
            return volume/6;
        }

        template <typename Collection>
        void computeVertices(Collection& result)
        {
            result.insert(result.end(), vertices_.begin(), vertices_.end());
        }

        template <typename Collection>
        void computeNeighbors(Collection& result)
        {
            std::transform(
                faces_.begin(), faces_.end(),
                std::back_inserter(result),
                [](Face f) {
                    return f.id;
                }
            );
        }

        void computeFaceData()
        {
            if (faceData_.size() > 0) {return;}

            for (auto f = faces_.begin(); f != faces_.end(); ++f) {
                faceData_.emplace_back(f.index(), weightedNormal(f.index()));
            }
        }

        Vector3 weightedNormal(FaceIndex fi)
        {
            Vector3 normal = Vector3(0, 0, 0);

            EdgeIndex start = startingEdge(fi);
            Vector3 anchor = vertices_[target(start)];
            EdgeIndex currentEdge = next(start);
            Vector3 curr = vertices_[target(currentEdge)] - anchor;
            currentEdge = next(currentEdge);
            for(; currentEdge != start; currentEdge = next(currentEdge)) {
                Vector3 prev = curr;
                curr = vertices_[target(currentEdge)] - anchor;
                normal += cross(prev, curr);
            }
            return normal;
        }

        double temporarilyComputeVolume()
        {
            bool alreadyComputed = (faceData_.size() > 0);
            double vol = computeVolume();
            if (!alreadyComputed) {faceData_.clear();}
            return vol;
        }

        double maximumNeighborDistance()
        {
            if (maxDistanceVertex_ == INVALID_VERTEX) {
                maximumNeighborDistance_ = 0;
                for (auto vi = vertices_.begin(); vi != vertices_.end(); ++vi) {
                    double d = mag2(*vi);
                    if (d > maximumNeighborDistance_) {
                        maximumNeighborDistance_ = d;
                        maxDistanceVertex_ = vi.index();
                    }
                }
                maximumNeighborDistance_ *= 4*(1 + 1e-6);
            }
            return maximumNeighborDistance_;
        }

        void buildCube(double xmin, double xMAX,
                       double ymin, double yMAX,
                       double zmin, double zMAX)
        {
            VERIFY(edges_.size() == 0);
            VERIFY(vertices_.size() == 0);
            VERIFY(faces_.size() == 0);

            VERIFY(xmin < xMAX);
            VERIFY(ymin < yMAX);
            VERIFY(zmin < zMAX);

            VERIFY_EXIT(approxEq(temporarilyComputeVolume(), (xMAX-xmin)*(yMAX-ymin)*(zMAX-zmin)));
            VERIFICATION_EXIT(verifyIsValidPolyhedron();)

            // In general, the faces will be referred to using single characters,
            // which are:
            //   Front         : F
            //   Right         : R
            //   Back          : B
            //   Left          : L
            //   Up (Top)      : U
            //   Down (Bottom) : D

            // Each vertex is written as the three faces it touches.
            // The overall layout of the cube is:
            //
            //    BUL-------BUR
            //    /|        /|
            //   / |       / |
            // FUL-------FUR |
            //  |  |      |  |
            //  | BDL-----|-BDR
            //  | /       | /
            //  |/        |/
            // FDL-------FDR

            enum Vertices {
                FDL,
                FDR,
                FUR,
                FUL,
                BDL,
                BDR,
                BUR,
                BUL
            };

            auto V = [&](Vertices VERIFICATION(expected),
                         double x, double y, double z)
            {
                VERIFICATION(VertexIndex v =) createVertex(x, y, z);
                VERIFY(
                    static_cast<VertexPool::SizeType>(v)
                      ==
                    static_cast<VertexPool::SizeType>(expected)
                );
            };

            //    BUL-------BUR
            //    /|        /|
            //   / |       / |
            // FUL-------FUR |
            //  |  |      |  |
            //  | BDL-----|-BDR
            //  | /       | /
            //  |/        |/
            // FDL-------FDR

            V(FDL, xmin, ymin, zmin);
            V(FDR, xMAX, ymin, zmin);
            V(FUR, xMAX, ymin, zMAX);
            V(FUL, xmin, ymin, zMAX);

            V(BDL, xmin, yMAX, zmin);
            V(BDR, xMAX, yMAX, zmin);
            V(BUR, xMAX, yMAX, zMAX);
            V(BUL, xmin, yMAX, zMAX);

            enum Edges {
                FU, FL, FD, FR,
                RU, RF, RD, RB,
                BU, BR, BD, BL,
                LU, LB, LD, LF,
                UF, UR, UB, UL,
                DF, DL, DB, DR
            };

            auto E = [&](FaceIndex face, Edges VERIFICATION(expected),
                         Edges flip, Vertices target, Edges next)
            {
                VERIFICATION(EdgeIndex e =) createEdge(
                    face,
                    VertexIndex {static_cast<VertexPool::SizeType>(target)},
                    EdgeIndex {static_cast<EdgePool::SizeType>(flip)},
                    EdgeIndex {static_cast<EdgePool::SizeType>(next)}
                );
                VERIFY(
                    static_cast<EdgePool::SizeType>(e)
                      ==
                    static_cast<EdgePool::SizeType>(expected)
                );
            };

            auto F = [&](Edges startingEdge) -> FaceIndex
            {
                return createFace(
                    -1,
                    EdgeIndex {static_cast<EdgePool::SizeType>(startingEdge)}
                );
            };

            VERIFICATION(
                auto edgeOrder = [&](Edges e1, Edges e2, Edges e3, Edges e4) -> bool
                {
                    EdgeIndex es[4] = {
                        EdgeIndex {static_cast<EdgePool::SizeType>(e1)},
                        EdgeIndex {static_cast<EdgePool::SizeType>(e2)},
                        EdgeIndex {static_cast<EdgePool::SizeType>(e3)},
                        EdgeIndex {static_cast<EdgePool::SizeType>(e4)}
                    };
                    for (int i = 0; i < 4; ++i) {
                        if (next(es[i]) != es[(i+1)%4]) {return false;}
                    }
                    return true;
                };
            )

            // Front Face
            //            FU
            //    FUL <--------- FUR
            //     |              ^
            //     |              |
            //  FL |              | FR
            //     |              |
            //     v              |
            //    FDL ---------> FDR
            //            FD

            VERIFY_EXIT(edgeOrder(FU, FL, FD, FR));
            FaceIndex f = F(FU);
            E(f, FU, UF, FUL, FL);
            E(f, FL, LF, FDL, FD);
            E(f, FD, DF, FDR, FR);
            E(f, FR, RF, FUR, FU);

            // Right Face
            //            RU
            //    FUR <--------- BUR
            //     |              ^
            //     |              |
            //  RF |              | RB
            //     |              |
            //     v              |
            //    FDR ---------> BDR
            //            RD

            VERIFY_EXIT(edgeOrder(RU, RF, RD, RB));
            f = F(RU);
            E(f, RU, UR, FUR, RF);
            E(f, RF, FR, FDR, RD);
            E(f, RD, DR, BDR, RB);
            E(f, RB, BR, BUR, RU);

            // Back Face
            //            BU
            //    BUR <--------- BUL
            //     |              ^
            //     |              |
            //  BR |              | BL
            //     |              |
            //     v              |
            //    BDR ---------> BDL
            //            BD

            VERIFY_EXIT(edgeOrder(BU, BR, BD, BL));
            f = F(BU);
            E(f, BU, UB, BUR, BR);
            E(f, BR, RB, BDR, BD);
            E(f, BD, DB, BDL, BL);
            E(f, BL, LB, BUL, BU);

            // Left Face
            //            LU
            //    BUL <--------- FUL
            //     |              ^
            //     |              |
            //  LB |              | LF
            //     |              |
            //     v              |
            //    BDL ---------> FDL
            //            LD

            VERIFY_EXIT(edgeOrder(LU, LB, LD, LF));
            f = F(LU);
            E(f, LU, UL, BUL, LB);
            E(f, LB, BL, BDL, LD);
            E(f, LD, DL, FDL, LF);
            E(f, LF, FL, FUL, LU);

            // Up Face
            //            UB
            //    BUL <--------- BUR
            //     |              ^
            //     |              |
            //  UL |              | UR
            //     |              |
            //     v              |
            //    FUL ---------> FUR
            //            UF

            VERIFY_EXIT(edgeOrder(UF, UR, UB, UL));
            f = F(UF);
            E(f, UF, FU, FUR, UR);
            E(f,UR, RU, BUR, UB);
            E(f, UB, BU, BUL, UL);
            E(f, UL, LU, FUL, UF);

            // Down Face
            //            DB
            //    BDL ---------> BDR
            //     ^              |
            //     |              |
            //  DL |              | DR
            //     |              |
            //     |              V
            //    FDL <--------- FDR
            //            DF

            VERIFY_EXIT(edgeOrder(DF, DL, DB, DR));
            f = F(DF);
            E(f, DF, FD, FDL, DL);
            E(f, DL, LD, BDL, DB);
            E(f, DB, BD, BDR, DR);
            E(f, DR, RD, FDR, DF);

            // Finally, we need to pick an arbitrary edge to be
            // our handle to the whole polyhedron.
            root_ = EdgeIndex {FU};
        }

        EdgeIndex findOutgoingEdge(const Plane& plane) {
            auto signedDistance = [&](VertexIndex vi) -> double {
                // return plane.signedDistance(vertices_[vi]);
                return plane.offset(vertices_[vi]);
            };

            auto location = [&](double distance) -> Plane::Location {
                return plane.location(distance - plane.planeOffset);
            };
            // We start by finding some OUTSIDE vertex, i.e. one that will
            // be cut off. If we can't find one, then there will be no outgoing
            // edge (and no need to cut), so we return INVALID_EDGE.
            // Starting from the OUTSIDE vertex, we find an INSIDE one.
            // We will then have an edge pointing from either OUTSIDE or
            // directly on the cutting plane (INCIDENT) to INSIDE.
            // The flip of this edge will be the result.

            // We find an outside vertex by starting at some vertex and iteratively
            // moving to a neighbor that is closer to being outside.
            // If we find a local maximum that is not outside, then we are done.
            // THIS ONLY WORKS BECAUSE THE POLYHEDRON IS CONVEX.

            // We have NOT definitely shown that this is entirely robust to
            // floating point errors.
            // However, we believe such an error should only happen
            // when the cutting plane is close to parallel with some of the
            // faces of the polyhedron.
            // As long as the polyhedron is CLOSE to convex, any missed cuts
            // would thus be very small.

            // We start using an arbitrary edge (root_). Between that edge
            // and its flip, we pick the edge that is moving from less to more
            // outside.
            EdgeIndex edgeToCurrent = root_, edgeToPrevious = flip(edgeToCurrent);
            double currentDistance;
            {
                currentDistance = signedDistance(target(edgeToCurrent));
                double previousDistance = signedDistance(target(edgeToPrevious));

                // Increasing distance is more outside.
                if (previousDistance > currentDistance) {
                    swap(edgeToCurrent, edgeToPrevious);
                    currentDistance = previousDistance;
                }
            }

            // Now we enter the hill-climbing portion to find an outside vertex.
            while (location(currentDistance) != Plane::OUTSIDE) {
                VERIFY_INVARIANT(edgeToCurrent == flip(edgeToPrevious));
                VERIFY_INVARIANT(currentDistance == signedDistance(target(edgeToCurrent)));

                EdgeIndex edgeToNeighbor = next(edgeToCurrent);
                for (; edgeToNeighbor != edgeToPrevious;
                     edgeToNeighbor = next(flip(edgeToNeighbor)))
                {
                    double neighborDistance = signedDistance(target(edgeToNeighbor));

                    // More distance means "more outside"
                    if (neighborDistance > currentDistance) {
                        edgeToCurrent = edgeToNeighbor;
                        currentDistance = neighborDistance;
                        edgeToPrevious = flip(edgeToCurrent);

                        break;
                    }
                }

                // If the loop fell through without breaking, then
                // we are at a local maximum. Hence, there are no
                // outside vertices and we don't have to cut.
                if (edgeToNeighbor == edgeToPrevious) {
                    // Return INVALID_EDGE to signify that we do not
                    // need to cut.
                    return INVALID_EDGE;
                }
            }

            // Now the current vertex is an outside vertex.
            // Next, we start OUTSIDE and then find a vertex that is INSIDE.
            // This means we'll end up with an edge pointing inward across the
            // boundary. (Edge target is INSIDE, edge source is either
            // INCIDENT or OUTSIDE.)

            // It's best if we start with an edge going from less to more INSIDE.
            // edgeToCurrent is an edge going from less to more OUTSIDE, so
            // we just reverse directions and swap edgeToCurrent with edgeToPrevious
            // (updating current distance, of course.)
            VERIFY(edgeToPrevious == flip(edgeToCurrent));
            swap(edgeToCurrent, edgeToPrevious);
            currentDistance = signedDistance(target(edgeToCurrent));
            while (location(currentDistance) != Plane::INSIDE) {
                VERIFY_INVARIANT(edgeToCurrent == flip(edgeToPrevious));
                VERIFY_INVARIANT(currentDistance == signedDistance(target(edgeToCurrent)));

                EdgeIndex edgeToNeighbor = next(edgeToCurrent);
                for (; edgeToNeighbor != edgeToPrevious;
                     edgeToNeighbor = next(flip(edgeToNeighbor)))
                {
                    double neighborDistance = signedDistance(target(edgeToNeighbor));

                    // Less distance means "more inside"
                    if (neighborDistance < currentDistance) {
                        edgeToCurrent = edgeToNeighbor;
                        currentDistance = neighborDistance;
                        edgeToPrevious = flip(edgeToCurrent);

                        break;
                    }
                }

                // There MUST be some INSIDE vertex. So if we reach
                // a local minimum without finding an INSIDE vertex,
                // either we're cutting by something very close -
                // as in, within machine epsilon close -
                // or our polyhedron is not convex.
                VERIFY(edgeToNeighbor != edgeToPrevious);
            }
            return flip(edgeToCurrent);
        }

        bool cutWithPlane(const int faceid, const Plane plane)
        {
            VERIFICATION_INVARIANT(verifyIsValidPolyhedron();)
            VERIFY(faceData_.size() == 0);

            auto location = [&](VertexIndex vi) -> Plane::Location {
                return plane.location(vertices_[vi]);
            };

            auto createIntersection = [&](VertexIndex a, VertexIndex b) -> VertexIndex
            {
                return createVertex(plane.intersection(vertices_[a], vertices_[b]));
            };

            // =================================================================

            EdgeIndex firstOutgoingEdge = findOutgoingEdge(plane);
            if (firstOutgoingEdge == INVALID_EDGE) { return false; }

            // Make sure that root_ is an edge that we know
            // won't end up getting destroyed during cutting.
            root_ = firstOutgoingEdge;

            EdgeIndex outgoingEdge = firstOutgoingEdge;
            VertexIndex previousIntersection = INVALID_VERTEX;

            EdgeIndex firstOutsideFaceEdge = createEdge();
            FaceIndex outsideFace = createFace(faceid, firstOutsideFaceEdge);
            face(firstOutsideFaceEdge) = outsideFace;
            EdgeIndex outsideFaceEdge = firstOutsideFaceEdge;

            VERIFY_INVARIANT(verticesToDestroy_.size() == 0);
            VERIFY_INVARIANT(edgesToDestroy_.size() == 0);

            do {
                VertexIndex previousVertex = target(outgoingEdge);
                VERIFY(location(previousVertex) != Plane::INSIDE);
                verticesToDestroy_.push_back(previousVertex);

                EdgeIndex currentEdge = next(outgoingEdge);
                VertexIndex currentVertex = target(currentEdge);

                bool needToCut = (location(previousVertex) == Plane::OUTSIDE);

                while (location(currentVertex) != Plane::INSIDE) {
                    needToCut = true;
                    verticesToDestroy_.push_back(currentVertex);
                    edgesToDestroy_.push_back(currentEdge);

                    previousVertex = currentVertex;
                    currentEdge = next(currentEdge);
                    currentVertex = target(currentEdge);
                }

                target(outgoingEdge) = previousIntersection;

                if (needToCut) {
                    VertexIndex currentIntersection = createIntersection(previousVertex, currentVertex);
                    FaceIndex currentFace = face(outgoingEdge);
                    // Make sure that the face's starting edge is one that
                    // we know won't end up getting destroyed
                    startingEdge(currentFace) = outgoingEdge;

                    EdgeIndex bridge = createEdge(
                        /* face */ currentFace,
                        /* vertex */ currentIntersection,
                        /* flip edge */ outsideFaceEdge,
                        /* next edge */ currentEdge
                    );

                    flip(outsideFaceEdge) = bridge;
                    next(outgoingEdge) = bridge;

                    outsideFaceEdge = createEdge(
                        /* face */ outsideFace,
                        /* vertex */ currentIntersection,
                        /* flip edge */ INVALID_EDGE,
                        /* next edge */ outsideFaceEdge
                    );

                    previousIntersection = currentIntersection;
                }

                outgoingEdge = flip(currentEdge);
            } while (outgoingEdge != firstOutgoingEdge);

            for (; target(outgoingEdge) == INVALID_VERTEX;
                   outgoingEdge = next(flip(outgoingEdge)))
            {
                target(outgoingEdge) = previousIntersection;
            }

            // firstOutsideFaceEdge and the final outsideFaceEdge
            // are actually supposed to end up becoming the same edge.
            // Since some things point to firstOutsideFaceEdge,
            // but nothing points to the final outsideFaceEdge,
            // we copy over the information needed to complete
            // firstOutsideFaceEdge and then destroy the other.
            VERIFY(firstOutsideFaceEdge != outsideFaceEdge);
            target(firstOutsideFaceEdge) = target(outsideFaceEdge);
            next(firstOutsideFaceEdge) = next(outsideFaceEdge);
            destroy(outsideFaceEdge);

            // Cleanup
            for (VertexIndex vi : verticesToDestroy_) {
                if (vertices_.active(vi)) { destroy(vi); }
            }
            verticesToDestroy_.clear();

            if (edgesToDestroy_.size() > 0) {
                // The whole region to be destroyed is a single
                // connected component.

                // First, we delete a perimeter around the region that will
                // be destroyed, marking out the boundary.
                // At the same time, we store some edges on the inside
                // of the region to be deleted to serve as our starting points.
                // (Any single such edge would suffice, but we can't be sure
                // which edges are part of the perimeter being destroyed right
                // now and which aren't, so we have to store all of the options.)
                std::size_t originalEdgesToDestroy = edgesToDestroy_.size();
                for (std::size_t i = 0; i < originalEdgesToDestroy; ++i) {
                    EdgeIndex ei = edgesToDestroy_[i];
                    // Store an edge inside the deletion region
                    edgesToDestroy_.push_back(flip(ei));
                    destroy(ei);
                }

                // Next, we do an exhaustive search of the region to destroy,
                // and DESTROY IT!
                for (std::size_t i = originalEdgesToDestroy;
                                 i < edgesToDestroy_.size();
                                 ++i)
                {
                    // I had previously thought that:
                    // Since the region to destroy is a connected component,
                    // we only need a single non-perimeter edge to successfully
                    // delete the whole thing.
                    // It would be wasteful, but safe, to not break and just
                    // markSweep the rest as well.

                    //     if (markSweep(edgesToDestroy_[i])) { break; }

                    // However, this turns out to be false.
                    // I do not know why.
                    // It seems that the region to destroy can in general
                    // be several connected components
                    // (well, OK, I understand that describes any
                    //  sensible topological region, but...)
                    // So we have to sweep with EVERY edge on the boundary.
                    // No short-circuiting.
                    markSweep(edgesToDestroy_[i]);
                }

                edgesToDestroy_.clear();
            }

            return true;
        }

        bool markSweep(EdgeIndex ei)
        {
            if (!edges_.active(ei)) { return false; }

            HalfEdge edge = edges_[ei];
            destroy(ei);
            markSweep(edge.flip);
            markSweep(edge.next);
            if (vertices_.active(edge.target)) {
                destroy(edge.target);
            }
            if (faces_.active(edge.face)) {
                destroy(edge.face);
            }
            return true;
        }

        friend std::ostream& operator<< (std::ostream& out, const Polyhedron& poly)
        {
            const auto& edges_ = poly.edges_;
            const auto& vertices_ = poly.vertices_;
            const auto& faces_ = poly.faces_;

            out << "POLYHEDRON" << std::endl;
            out << "EDGES" << std::endl;
            for (EdgePool::SizeType i = 0; i < edges_.size(); ++i) {
                EdgeIndex ei {i};
                out << "Edge chunk " << i << ": ";
                if (edges_.active(ei)) {
                    out << "* Active" << std::endl;
                    out << "  Face: " << poly.face(ei)
                        << (poly.face(ei) != INVALID_FACE && faces_.active(poly.face(ei)) ? " *" : " X") << std::endl;
                    out << "  Next: " << poly.next(ei)
                        << (poly.next(ei) != INVALID_EDGE && edges_.active(poly.next(ei)) ? " *" : " X") << std::endl;
                    out << "  Flip: " << poly.flip(ei)
                        << (poly.flip(ei) != INVALID_EDGE && edges_.active(poly.flip(ei)) ? " *" : " X") << std::endl;
                    out << "  Target: " << poly.target(ei)
                        << (poly.target(ei) != INVALID_VERTEX && vertices_.active(poly.target(ei)) ? " *" : " X") << std::endl;
                }
                else {
                    out << "X Inactive" << std::endl;
                }
            }
            out << "VERTICES" << std::endl;
            for (VertexPool::SizeType i = 0; i < vertices_.size(); ++i) {
                VertexIndex vi {i};
                out << "Vertex chunk " << i << ": ";
                if (vertices_.active(vi)) {
                    out << "* Active" << std::endl;
                }
                else {
                    out << "X Inactive" << std::endl;
                }
            }
            out << "FACES" << std::endl;
            for (FacePool::SizeType i = 0; i < faces_.size(); ++i) {
                FaceIndex fi {i};
                out << "Face " << i << ": ";
                if (faces_.active(fi)) {
                    out << "* Active" << std::endl;
                    out << "  Start: " << poly.startingEdge(fi)
                        << ((poly.startingEdge(fi) != INVALID_EDGE && edges_.active(poly.startingEdge(fi))) ? " *" : " X") << std::endl;
                }
                else {
                    out << "X Inactive" << std::endl;
                }
            }

            return out;
        }
    };

}
