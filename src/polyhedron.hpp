/**
 * \file polyhedron.hpp
 *
 * \author 2015-2016 Sandia Clinic Team
 *
 * \brief Contains the Polyhedron class which does the actual Voronoi calculations.
 *
 */

#pragma once

#include <iostream>
#include <vector>
#include <functional>
#include <unordered_set>
#include <iterator>
#include <algorithm>
#include <cmath>

#include "vectormath.hpp"
#include "structpool.hpp"
#include "verification.hpp"

namespace hmc {
    constexpr double TOLERANCE = 1e-14;

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


    /**
     * \struct HalfEdge
     * \brief
     *  The implementation of a half-edge structure. 2 half-edges make up an edge in the polyhedron
     *  each pointing to one of the two vertices that make up that edge.
     */
    struct HalfEdge {
        EdgeIndex flip;             ///< The other HalfEdge in the edge
        EdgeIndex next;             ///< The next HalfEdge in the polyhedron face
        VertexIndex target;         ///< The vertex this HalfEdge points to
        FaceIndex face;             ///< The polyhedral face that this HalfEdge is in


        /**
         * \brief Constructor for HalfEdge
         *
         * \param[in] f        FaceIndex to access correct face from FacePool.
         * \param[in] v        VertexIndex to access the target vertex from VertexPool.
         * \param[in] flip     EdgeIndex to access the HalfEdge flip from EdgePool.
         * \param[in] next     EdgeIndex to access the next HalfEdge from EdgePool.
         */
        HalfEdge() = default;
        HalfEdge(FaceIndex f) : flip(), next(), target(), face(f) {}
        HalfEdge(FaceIndex f, VertexIndex v) : flip(), next(), target(v), face(f) {}
        HalfEdge(FaceIndex f, VertexIndex v, EdgeIndex flip, EdgeIndex next)
            : flip(flip), next(next), target(v), face(f)
        { /* Done */ }
    };

    /**
     * \struct Face
     *
     * \brief The polyhedral face on the polyhedron
     */
    struct Face {
        int id;                         ///< The index of the neighboring particle that made the face
        EdgeIndex startingEdge;         ///< The first HalfEdge in the face
        Vector3 unitNormal;


        /**
         * \brief Constructor for a Face
         *
         * \param[in] id                The index of the neighboring particle that made the face.
         * \param[in] startingEdge      EdgeIndex to access the first HalfEdge from the EdgePool.
         */
        Face(int id)
            : id(id),
              VERIFICATION   ( startingEdge(INVALID_EDGE) )
              NO_VERIFICATION( startingEdge() )
        { /* Done */ }

        Face(int id, EdgeIndex startingEdge, Vector3 unitNormal)
            : id(id), startingEdge(startingEdge), unitNormal(unitNormal) {}
    };

    /**
     * \class Polyhedron
     *
     * \brief The cell from the Voronoi tessellation. Most of the computation work
     *        necessary to get the desired information is done here.
     * \remark This should not be used directly. Instead, interface through the Cell class.
     */
    class Polyhedron {
    public:

        /**
         * \struct FaceData
         */
        struct FaceData {
            FaceIndex face;           ///< FaceIndex to the desired face in FacePool.
            Vector3 weightedNormal;   ///< The vector normal to the face.

            /**
             * \brief Constructor for FaceData
             *
             * \param[in] face               FaceIndex to the associated face in FacePool
             * \param[in] weightedNormal     The normal vector to the face.
             */
            FaceData(FaceIndex face, Vector3 weightedNormal)
                : face(face), weightedNormal(weightedNormal)
            { /* Done */ }
        };

        EdgeIndex root_ = INVALID_EDGE; ///< The root edge of the polyhedron. Start as the invalid edge

        EdgePool edges_; ///< The memory pool that stores the edges
        VertexPool vertices_; ///< The memory pool that stores the vertices
        FacePool faces_; ///< The memory pool that stores the faces

        std::vector<FaceData> faceData_; ///< A vector storing all of the face data

        std::vector<VertexIndex> verticesToDestroy_; ///< A vector of vertices that need to be deleted
        std::vector<EdgeIndex> edgesToDestroy_; ///< A vector of edges that need to be deleted

        VertexIndex maxDistanceVertex_ = INVALID_VERTEX; ///< The furthest vertex from the particle
        double maximumNeighborDistance_; ///< The furthest another particle can be from the current particle and still cut the polyhedron


        void clear() ///< Removes all edges, faces, and vertices from polyhedron
        {
            root_ = INVALID_EDGE;
            edges_.clear();
            vertices_.clear();
            faces_.clear();
            faceData_.clear();
        }

        bool isClear() ///< Verifies that the polyhedron has been cleared
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

        void clearComputation() ///< Only clear faceData
        {
            faceData_.clear();
        }

        /**
         * \brief Get the flip HalfEdge of the input HalfEdge
         *
         * \param[in] ei        The EdgeIndex of a HalfEdge
         * \result          The EdgeIndex of the flip of the input HalfEdge
         */
        EdgeIndex& flip(EdgeIndex ei) {return edges_[ei].flip;}

        /**
         * \brief Get the flip HalfEdge of the input HalfEdge
         *
         * \param[in] ei        The EdgeIndex of a HalfEdge
         * \result          The const EdgeIndex of the flip of the input HalfEdge
         */
        const EdgeIndex& flip(EdgeIndex ei) const {return edges_[ei].flip;}


        /**
         * \brief Get the next HalfEdge of the input HalfEdge
         *
         * \param[in] ei        The EdgeIndex of a HalfEdge
         * \result          The EdgeIndex of the next of the input HalfEdge
         */
        EdgeIndex& next(EdgeIndex ei) {return edges_[ei].next;}

        /**
         * \brief Get the next HalfEdge of the input HalfEdge
         *
         * \param[in] ei        The EdgeIndex of a HalfEdge
         * \result          The const EdgeIndex of the next of the input HalfEdge
         */
        const EdgeIndex& next(EdgeIndex ei) const {return edges_[ei].next;}

        /**
         * \brief Get the next HalfEdge of the input HalfEdge
         *
         * \param ei        The EdgeIndex of a HalfEdge
         * \result          The const EdgeIndex of the next of the input HalfEdge
         */
        EdgeIndex& nextWithSameTarget(EdgeIndex ei) {return flip(next(ei));}

        /**
         * \brief Get the next HalfEdge of the input HalfEdge
         *
         * \param ei        The EdgeIndex of a HalfEdge
         * \result          The const EdgeIndex of the next of the input HalfEdge
         */
        const EdgeIndex& nextWithSameTarget(EdgeIndex ei) const {return flip(next(ei));}

        /**
         * \brief Get the next HalfEdge of the input HalfEdge
         *
         * \param ei        The EdgeIndex of a HalfEdge
         * \result          The const EdgeIndex of the next of the input HalfEdge
         */
        EdgeIndex& nextWithSameSource(EdgeIndex ei) {return next(flip(ei));}

        /**
         * \brief Get the next HalfEdge of the input HalfEdge
         *
         * \param ei        The EdgeIndex of a HalfEdge
         * \result          The const EdgeIndex of the next of the input HalfEdge
         */
        const EdgeIndex& nextWithSameSource(EdgeIndex ei) const {return next(flip(ei));}


        /**
         * \brief Get the target vertex of the input HalfEdge
         *
         * \param[in] ei        The EdgeIndex of a HalfEdge
         * \result          The VertexIndex of the target of the input HalfEdge
         */
        VertexIndex& target(EdgeIndex ei) {return edges_[ei].target;}

        /**
         * \brief Get the target vertex of the input HalfEdge
         *
         * \param[in] ei        The EdgeIndex of a HalfEdge
         * \result          The const VertexIndex of the target of the input HalfEdge
         */
        const VertexIndex& target(EdgeIndex ei) const {return edges_[ei].target;}

        /**
         * \brief Get the source vertex of the input HalfEdge
         *
         * \param[in] ei        The EdgeIndex of a HalfEdge
         * \result          The VertexIndex of the source of the input HalfEdge
         */
        VertexIndex& source(EdgeIndex ei) {return target(flip(ei));}

        /**
         * \brief Get the source vertex of the input HalfEdge
         *
         * \param[in] ei        The EdgeIndex of a HalfEdge
         * \result          The const VertexIndex of the source of the input HalfEdge
         */
        const VertexIndex& source(EdgeIndex ei) const {return target(flip(ei));}

        /**
         * \brief Get the face of the input HalfEdge
         *
         * \param[in] ei        The EdgeIndex of a HalfEdge
         * \result          The FaceIndex of the face of the input HalfEdge
         */
        FaceIndex& face(EdgeIndex ei) {return edges_[ei].face;}

        /**
         * \brief Get the face of the input HalfEdge
         *
         * \param[in] ei        The EdgeIndex of a HalfEdge
         * \result          The const FaceIndex of the face of the input HalfEdge
         */
        const FaceIndex& face(EdgeIndex ei) const {return edges_[ei].face;}

        /**
         * \brief Get the first HalfEdge of the input face
         *
         * \param[in] fi        The FaceIndex of a face
         * \result          The EdgeIndex of the first HalfEdge in the face
         */
        EdgeIndex& startingEdge(FaceIndex fi) {return faces_[fi].startingEdge;}

        /**
         * \brief Get the first HalfEdge of the input face
         *
         * \param[in] fi        The FaceIndex of a face
         * \result          The const EdgeIndex of the first HalfEdge in the face
         */
        const EdgeIndex& startingEdge(FaceIndex fi) const {return faces_[fi].startingEdge;}


        /**
         * \brief Create edges in edge memory pool
         *
         * \param[in] args  Any number of edge constructor calls.
         * \return        The EdgeIndex of the newly created edge.
         */
        template <typename... Edge_Constructor_Args>
        EdgeIndex createEdge(Edge_Constructor_Args... args)
        {
            return edges_.create(args...);
        }

        /**
         * \brief Create vertices in vertex memory pool
         *
         * \param[in] args  Any number of vertex constructor calls.
         * \return        The VertexIndex of the newly created vertex.
         */
        template <typename... Vertex_Constructor_Args>
        VertexIndex createVertex(Vertex_Constructor_Args... args)
        {
            VERIFICATION(
                Vector3 v {args...};
                VERIFY(std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z));
            )
            return vertices_.create(args...);
        }

        /**
         * \brief Create faces in face memory pool
         *
         * \param[in] args  Any number of face constructor calls.
         * \return        The FaceIndex of the newly created face.
         */
        template <typename... Face_Constructor_Args>
        FaceIndex createFace(Face_Constructor_Args... args)
        {
            return faces_.create(args...);
        }

        /**
         * \brief Delete edge in edge memory pool
         *
         * \param[in] ei  The EdgeIndex to be removed from the edge memory pool.
         */
        void destroy(EdgeIndex ei)
        {
            edges_.destroy(ei);
        }

        /**
         * \brief Delete vertex in vertex memory pool
         *
         * \param[in] vi  The VertezIndex to be removed from the vertex memory pool.
         */
        void destroy(VertexIndex vi)
        {
            vertices_.destroy(vi);
            if (vi == maxDistanceVertex_) {
                maxDistanceVertex_ = INVALID_VERTEX;
            }
        }

        /**
         * \brief Delete face in face memory pool
         *
         * \param[in] fi  The FaceIndex to be removed from the face memory pool.
         */
        void destroy(FaceIndex fi)
        {
            faces_.destroy(fi);
        }

        /**
         * \brief Verification that the polyhedron is a valid, convex polyhedron.
         *        For debugging purposes only.
         *
         * \deprecated Should be wrapped in Verification(). This will allow it to be
         *        used when Verification is turned on.
         *
         */
        bool isValid()
        {
            if (!(root_ != INVALID_EDGE)) {return false;}
            for (auto eit = edges_.begin(); eit != edges_.end(); ++eit) {
                if (!(eit->target != INVALID_VERTEX)) {return false;}
                if (!(flip(eit->flip) == eit.index())) {return false;}
                if (!(face(eit->next) == eit->face)) {return false;}
                if (!(face(eit->flip) != eit->face)) {return false;}
                if (!(flip(eit->next) != INVALID_EDGE)) {return false;}
                if (!(source(eit->next) == eit->target)) {return false;}
            }
            for (auto fit = faces_.begin(); fit != faces_.end(); ++fit) {
                if (!(face(fit->startingEdge) == fit.index())) {return false;}
            }


            std::unordered_set<EdgeIndex> reachableEdges;
            std::unordered_set<VertexIndex> reachableVertices;
            std::unordered_set<FaceIndex> reachableFaces;

            std::function<bool(VertexIndex)> markVertex = [&](VertexIndex vi) -> bool {
                if (!(vertices_.active(vi))) {return false;}
                reachableVertices.insert(vi);
                return true;
            };

            std::function<bool(EdgeIndex)> markEdge = [&](EdgeIndex ei) -> bool {
                if (reachableEdges.find(ei) != reachableEdges.end()) { return true; }
                reachableEdges.insert(ei);
                if (!edges_.active(ei)) {return false;}
                if (!markEdge(flip(ei))) {return false;}
                if (!markEdge(next(ei))) {return false;}
                if (!markVertex(target(ei))) {return false;}
                reachableFaces.insert(face(ei));
                return true;
            };


            if (!markEdge(root_)) {return false;}

            for (auto ei : reachableEdges) { if (!(edges_.active(ei))) {return false;} }
            for (auto vi : reachableVertices) { if (!(vertices_.active(vi))) {return false;} }
            for (auto fi : reachableFaces) { if (!(faces_.active(fi))) {return false;} }

            for (auto e = edges_.begin(); e != edges_.end(); ++e) {
                if (!(reachableEdges.find(e.index()) != reachableEdges.end())) {return false;}
            }
            for (auto v = vertices_.begin(); v != vertices_.end(); ++v) {
                if (!(reachableVertices.find(v.index()) != reachableVertices.end())) {return false;}
            }
            for (auto f = faces_.begin(); f != faces_.end(); ++f) {
                if (!(reachableFaces.find(f.index()) != reachableFaces.end())) {return false;}
            }

            // It's easy to get the handedness wrong (-> negative volume).
            if (!(temporarilyComputeVolume() > 0)) {return false;}

            if (!isClear()) {
                // Checking the Euler characteristic, which is probably
                // unnecessary but might catch SOMETHING odd.
                auto E = std::distance(edges_.begin(), edges_.end());
                auto V = std::distance(vertices_.begin(), vertices_.end());
                auto F = std::distance(faces_.begin(), faces_.end());
                if (!(V - E/2 + F == 2)) {return false;}
            }

            return true;
        }



        /**
         * \brief Shift all vertices in the polyhedron by a vector.
         *
         * \param[in] shift  The vector that the polyhedron is shifted by.
         */
        void translate(const Vector3 shift)
        {
            for (auto& v : vertices_) {
                v += shift;
            }
        }

        /**
         * \brief Computes the volume of the polyheron.
         *
         * \return The volume of the polyhedron.
         */
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

        /**
         * \brief Adds the vertices of the polyhedron to the input collection
         *
         * \param[out] result The collection that the vertices are appended to.
         */
        template <typename Collection>
        void computeVertices(Collection& result)
        {
            auto output = std::back_inserter(result);
			for (auto v : vertices_) {
				*output++ = v.x;
				*output++ = v.y;
				*output++ = v.z;
			}
        }

        /**
         * \brief Adds the indices of the neighbors of the polyhedron to the input collection
         *
         * \param[out] result The collection that the neighbor particle indices are appended to.
         */
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

        /**
         * \brief Computes the data for each face. Nothing is returned.
         *
         */
        void computeFaceData()
        {
            if (faceData_.size() > 0) {return;}

            for (auto f = faces_.begin(); f != faces_.end(); ++f) {
                faceData_.emplace_back(f.index(), weightedNormal(f.index()));
            }
        }

        /**
         * \brief Calculates the normal face vector
         *
         * \param[in] fi The FaceIndex of the face the normal is being calculated for.
         * \result       The normal vector for the face.
         */
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

        template <typename F>
        void temporarilyCompute(F computation)
        {
            bool alreadyComputed = (faceData_.size() > 0);
            computation(*this);
            if (!alreadyComputed) {faceData_.clear();}
        }

        /**
         * \brief Computes the volume of the polyhedron. If the face data
         *        was not already computed, erase the face data.
         *
         * \result The volume of the cell.
         */
        double temporarilyComputeVolume()
        {
            double volume;
            temporarilyCompute([&volume](Polyhedron& c){
                volume = c.computeVolume();
            });
            return volume;
        }



        /**
         * \brief Finds the distance of the furthest neighbor.
         *
         * \result The distance of the furthest neighbor.
         */
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

        void buildCube(BoundingBox bounds)
        {
            buildCube(bounds.low.x, bounds.high.x,
                      bounds.low.y, bounds.high.y,
                      bounds.low.z, bounds.high.z);
        }

        /**
         * \brief Constructs the initial polyhedron as a cube
         *
         * \param[in] xmin          The minimum x-coordinate of the uncut cell
         * \param[in] xMAX          The maximum x-coordinate of the uncut cell
         * \param[in] ymin          The minimum y-coordinate of the uncut cell
         * \param[in] yMAX          The maximum y-coordinate of the uncut cell
         * \param[in] zmin          The minimum z-coordinate of the uncut cell
         * \param[in] zMAX          The maximum z-coordinate of the uncut cell
         */
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

            VERIFY_EXIT(approxRelEq(temporarilyComputeVolume(), (xMAX-xmin)*(yMAX-ymin)*(zMAX-zmin), TOLERANCE));
            VERIFY_EXIT(isValid());

            clear();

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

            //    BUL-------BUR               +z
            //    /|        /|                |   +y
            //   / |       / |                |  /
            // FUL-------FUR |                | /
            //  |  |      |  |       -x_______|/_______+x
            //  | BDL-----|-BDR              /|
            //  | /       | /               / |
            //  |/        |/               /  |
            // FDL-------FDR             -y   |
            //                                -z

            V(FDL, xmin, ymin, zmin);
            V(FDR, xMAX, ymin, zmin);
            V(FUR, xMAX, ymin, zMAX);
            V(FUL, xmin, ymin, zMAX);

            V(BDL, xmin, yMAX, zmin);
            V(BDR, xMAX, yMAX, zmin);
            V(BUR, xMAX, yMAX, zMAX);
            V(BUL, xmin, yMAX, zMAX);

            //    BUL-------BUR               +z
            //    /|        /|                |   +y
            //   / |       / |                |  /
            // FUL-------FUR |                | /
            //  |  |      |  |       -x_______|/_______+x
            //  | BDL-----|-BDR              /|
            //  | /       | /               / |
            //  |/        |/               /  |
            // FDL-------FDR             -y   |
            //                                -z

            // Unit normal vectors for each face
            Vector3
                backNormal  { 0, -1,  0},
                frontNormal { 0,  1,  0},
                upNormal    { 0,  0,  1},
                downNormal  { 0,  0, -1},
                rightNormal { 1,  0,  0},
                leftNormal  {-1,  0,  0};

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

            auto F = [&](Edges startingEdge, Vector3 unitNormal) -> FaceIndex
            {
                return createFace(
                    -1,
                    EdgeIndex {static_cast<EdgePool::SizeType>(startingEdge)},
                    unitNormal
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
            FaceIndex f = F(FU, frontNormal);
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
            f = F(RU, rightNormal);
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
            f = F(BU, backNormal);
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
            f = F(LU, leftNormal);
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
            f = F(UF, upNormal);
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
            f = F(DF, downNormal);
            E(f, DF, FD, FDL, DL);
            E(f, DL, LD, BDL, DB);
            E(f, DB, BD, BDR, DR);
            E(f, DR, RD, FDR, DF);

            // Finally, we need to pick an arbitrary edge to be
            // our handle to the whole polyhedron.
            root_ = EdgeIndex {FU};
        }

        /**
         * \brief Find an edge going through the cutting plane
         *
         * \remark Not robust for floating point errors.
         */
        EdgeIndex findOutgoingEdge(const Plane& plane) {
            auto signedDistance = [&](VertexIndex vi) -> double {
                // return plane.signedDistance(vertices_[vi]);
                return plane.offset(vertices_[vi]);
            };

            auto location = [&](double distance) -> Plane::Location {
                return plane.location(distance - plane.planeOffset, TOLERANCE);
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
                // we are at a local maximum. Hence, ASSUMING CONVEXITY,
                // there are no outside vertices and we don't have to cut.
                if (edgeToNeighbor == edgeToPrevious) {
                    // This is only an optimization. We expect to cut next
                    // using other points that are close to the one we
                    // just cut with, and we want to start as close to
                    // lating cutting planes as possible.
                    // So we start out at the previous maximum.
                    root_ = edgeToCurrent;

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

                // So, really we have encountered an error state here.
                // We can still attempt to at least do SOMETHING, however.
                if (edgeToNeighbor == edgeToPrevious) {
                    auto locationV = [&](VertexIndex vi) -> Plane::Location {
                        return plane.location(vertices_[vi], TOLERANCE);
                    };

                    auto isOutgoingEdge = [&](EdgeIndex ei) -> bool {
                        return locationV(target(ei)) != Plane::INSIDE
                            && locationV(source(ei)) == Plane::INSIDE;
                    };

                    for (auto eit = edges_.begin(); eit != edges_.end(); ++eit) {
                        if (isOutgoingEdge(eit.index())) {
                            return eit.index();
                        }
                    }
                }

                // There is no outgoing edge. However, there ARE
                // outside vertices if it got to this point,
                // so this means there aren't any inside vertices.

                VERIFICATION(
                    auto someVertexIsInside = [&]() -> bool {
                        for (auto v : vertices_) {
                            if (plane.location(v, TOLERANCE)) {
                                return true;
                            }
                        }
                        return false;
                    };

                    VERIFY(someVertexIsInside());
                )

                // TODO: What should we do if we enter this error
                // state, where we're trying to cut off the entire
                // polyhedron?

                // At the moment, by simply falling through,
                // we end up infinite looping because cutWithPlane
                // expects to be able to find outside vertices.
            }
            return flip(edgeToCurrent);
        }

        /**
         * \brief Cuts the polyhedron's face with the given plane.
         *
         * \param[in] faceid    The index of the face to be cut.
         * \param[in] plane     The plane object of the cutting plane.
         * \return True if the plane cut the face.
         */
        bool cutWithPlane(const int faceid, const Plane plane, bool verbose = false)
        {
            VERIFY_INVARIANT(isValid());
            VERIFY(faceData_.size() == 0);

            auto location = [&](VertexIndex vi) -> Plane::Location {
                return plane.location(vertices_[vi], TOLERANCE);
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
            FaceIndex outsideFace = createFace(faceid, firstOutsideFaceEdge, plane.unitNormal);
            face(firstOutsideFaceEdge) = outsideFace;
            EdgeIndex outsideFaceEdge = firstOutsideFaceEdge;

            VERIFY_INVARIANT(verticesToDestroy_.size() == 0);
            VERIFY_INVARIANT(edgesToDestroy_.size() == 0);

            auto locationV = [&](VertexIndex vi) -> Plane::Location {
                return plane.location(vertices_[vi], TOLERANCE);
            };

            auto outV = [&](VertexIndex vi) -> void {
                if (vi == INVALID_VERTEX) {
                    std::cerr << "INVALID VERTEX";
                    return;
                }
                std::cerr << "[" << locationV(vi) << "] "
                    << "Vertex #" << vi << " " << vertices_[vi];
            };

            auto outE = [&](EdgeIndex ei) -> void {
                if (ei == INVALID_EDGE) {
                    std::cerr << "INVALID EDGE";
                    return;
                }
                std::cerr << "Edge #" << ei << " (" << face(ei) << "), "
                    << "flip " << flip(ei) << " (" << face(flip(ei)) << "): ";
                    outV(source(ei)); std::cerr << " --> "; outV(target(ei));
            };

            if (verbose) {
                std::cerr << "Beginning to cut. Creating face " << outsideFace << std::endl;
                std::cerr << "First outgoing is: "; outE(firstOutgoingEdge); std::cerr << std::endl;
            }

            do {
                VertexIndex previousVertex = target(outgoingEdge);
                VERIFY(location(previousVertex) != Plane::INSIDE);
                verticesToDestroy_.push_back(previousVertex);

                EdgeIndex currentEdge = next(outgoingEdge);
                VertexIndex currentVertex = target(currentEdge);

                if (verbose) {
                    std::cerr << std::endl;
                    std::cerr << "Cutting face " << face(outgoingEdge) << std::endl;
                    std::cerr << "Outgoing: "; outE(outgoingEdge); std::cerr << std::endl;
                }

                bool needToCut = (location(previousVertex) == Plane::OUTSIDE);

                Plane::Location currentLocation, previousLocation = location(previousVertex);
                while ((currentLocation = location(currentVertex)) != Plane::INSIDE) {
                    needToCut = true;
                    verticesToDestroy_.push_back(currentVertex);
                    edgesToDestroy_.push_back(currentEdge);

                    previousVertex = currentVertex;
                    previousLocation = currentLocation;
                    currentEdge = next(currentEdge);
                    currentVertex = target(currentEdge);
                }

                target(outgoingEdge) = previousIntersection;

                if (verbose) {
                    std::cerr << "Need to cut: " << (needToCut ? "Yes" : "No") << std::endl;
                }

                if (needToCut) {
                    VertexIndex currentIntersection;
                    if (previousLocation == Plane::INCIDENT) {
                        currentIntersection = createVertex(vertices_[previousVertex]);
                    }
                    else {
                        currentIntersection = createIntersection(previousVertex, currentVertex);
                    }
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

                    if (verbose) {
                        std::cerr << "Ingoing: "; outE(currentEdge); std::cerr << std::endl;
                        std::cerr << "Intersection: "; outV(currentIntersection); std::cerr << std::endl;
                        std::cerr << "Bridge: "; outE(bridge); std::cerr << std::endl;
                    }
                }

                outgoingEdge = flip(currentEdge);
            } while (outgoingEdge != firstOutgoingEdge);

            for (; target(outgoingEdge) == INVALID_VERTEX;
                   outgoingEdge = nextWithSameTarget(outgoingEdge))
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

        /**
         * \brief Delete portion of the polyhedron connected to input edge.
         *
         * \param[in] ei    The index of the input edge.
         *
         * \remark Should be used in the plane cutting algorithm to clean up portions cut off by a
         *         polyhedron. If used on an edge in the main part of the polyhedron, the entire
         *         polyhedron will be destroyed.
         *
         */
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

        /**
         * \brief Print function for the polyhedron
         *
         * \param[in] out      The ostream object
         * \param[in] poly     The polyhedron to be printed
         *
         * \remark Printed in the following format. Items in quotes are variable output
         *         (*i = an index, "a/b" = a or b)
         *
         * \code{.unparsed}
         * POLYHEDRON
         * EDGES
         * Edge chunk "ei": "* Active/X Inactive"
         *   Face: "fi" " * / X"
         *   Next: "ni" " * / X"
         *   Flip: "fi" " * / X"
         *   Target: "ti" " * / X"
         * .
         * .
         * .
         * VERTICES
         * Vertex chunk "vi": "* Active/X Inactive"
         * .
         * .
         * .
         * FACES
         * Face "fi": "* Active/X Inactive"
         *   Start: "fei" "* /X"
         * .
         * .
         * .
         * \endcode
         */
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

        void outputGnuplot(std::ostream& out, Vector3 shift = VECTOR_ZERO) const
        {
            for (auto eit = edges_.begin(); eit != edges_.end(); ++eit) {
                if (eit.index().value > eit->flip.value) {
                    Vector3 v = vertices_[eit->target] + shift;
                    out << v.x << " " << v.y << " " << v.z << std::endl;
                    v = vertices_[source(eit.index())] + shift;
                    out << v.x << " " << v.y << " " << v.z << std::endl;
                    out << std::endl << std::endl;
                }
            }
        }

        void outputCuttingGraph(std::ostream& out, const Plane& plane) const
        {
            auto location = [&](VertexIndex vi) -> Plane::Location {
                return plane.location(vertices_[vi], TOLERANCE);
            };

            auto outV = [&](VertexIndex vi) -> void {
                if (vi == INVALID_VERTEX) {
                    out << "INVALID VERTEX";
                    return;
                }
                out << "[" << location(vi) << "] "
                    << "Vertex #" << vi << " " << vertices_[vi];
            };

            auto outE = [&](EdgeIndex ei) -> void {
                if (ei == INVALID_EDGE) {
                    out << "INVALID EDGE";
                    return;
                }
                out << "Edge #" << ei << " (" << face(ei) << "), "
                    << "flip " << flip(ei) << " (" << face(flip(ei)) << "): ";
                    outV(source(ei)); out << " --> "; outV(target(ei));
            };

            out << "Cutting with " << plane << std::endl;
            out << "Root: "; outE(root_);

            for (auto fit = faces_.begin(); fit != faces_.end(); ++fit) {
                out << std::endl << std::endl;
                out << "Face #" << fit.index() << ":" << std::endl;
                auto ei = fit->startingEdge;
                do {
                    out << "        ";
                    outE(ei);
                    out << std::endl;
                    out << "    ";
                    outV(target(ei));
                    out << std::endl;
                    ei = next(ei);
                } while(ei != fit->startingEdge);
            }
        }
    };

}
