#pragma once

#include <iterator>
#include <utility>
#include <vector>
#include <limits>
#include <type_traits>

namespace hmc {
    // Temporary types that will be supplied by implementation code
    struct Point {
        Point(double x, double y, double z) {}
    };

    using VoronoiPolyhedron = int;



    // ParticleIndex is just an internal typedef for readability.
    // ParticleIndices are used to refer to and access particular
    // particles within a diagram.
    using ParticleIndex = std::size_t;

    // INVALID_INDEX is essentially the equivalent of nullptr
    // for ParticleIndices.
    constexpr ParticleIndex INVALID_INDEX =
        std::numeric_limits<ParticleIndex>::max();

    // TODO: Describe Group
    // TODO: int or unsigned int?
    using Group = int;

    constexpr Group NO_GROUP =
        std::numeric_limits<Group>::max();


    // Ids are arbitrary metadata put onto particles by
    // the user, for the user, for whatever reasons they have to
    // do so. The package itself is totally agnostic of its use.
    // TODO: int or unsigned int?
    using Id = int;



    // =============
    // |  Options  |
    // =============

    // These options are bitflags which are passed to various functions
    // to modify behavior. For example,

    //     cell.computeNeighbors(hmc::Options::INCLUDE_SELF |
    //                           hmc::Options::INCLUDE_BOUNDARY_NEIGHBORS);

    // Note that since they are used as bitflags, you must not assume that
    // the enum actually provides an exhaustive list. For example,
    // INCLUDE_SELF | INCLUDE_BOUNDARY_NEIGHBORS produces an object of type
    // Options that is nonetheless not equal to NONE, INCLUDE_SELF or any
    // of the other enumerated flags.

    // Therefore, you should compare options by using the `selected` function,
    // like so:

    //     if(selected(given_options, Options::INCLUDE_BOUNDARY_NEIGHBORS)) {
    //         // blah
    //     }

    // or equivalently

    //     if(selected(Options::INCLUDE_BOUNDARY_NEIGHBORS, given_options)) {
    //         // blah
    //     }

    // `selected` merely masks one bitflag with the other,
    // and produces a boolean from the result.
    // You CANNOT use `selected` to check against Options::NONE, which will
    // simply always return `false`. Instead, in that instance you are
    // trying to get an exact equality, so you would just use `operator==`.

    enum class Options : int {
        NONE = 0,
        INCLUDE_SELF = 1 << 0,
        INCLUDE_BOUNDARY_NEIGHBORS = 1 << 1,
    };

// Temporarily-defined macros to make the code cleaner.

// OPTIONS_NUMERIC will evaluate to the underlying numeric type of
// the enum, allowing us to use bitwise operations.
#define OPTIONS_NUMERIC std::underlying_type<Options>::type

// Generate code for an ordinary binary operator
#define HMC_OPTIONS_GEN_OPERATOR_BIN(op) \
inline Options operator op (Options lhs, Options rhs) \
{ \
    return static_cast<Options>( \
        static_cast<OPTIONS_NUMERIC>(lhs) \
        op \
        static_cast<OPTIONS_NUMERIC>(rhs) \
    ); \
}

// Generate code for an update binary operator (e.g., x &= y)
#define HMC_OPTIONS_GEN_OPERATOR_EQ(op) \
inline Options& operator op ## = (Options& lhs, Options rhs) \
{ \
    lhs op ## = static_cast<Options>( \
        static_cast<OPTIONS_NUMERIC>(lhs) \
        op \
        static_cast<OPTIONS_NUMERIC>(rhs) \
    ); \
    return lhs; \
}

// Generate both the normal and update operator for a binary operator
#define HMC_OPTIONS_GEN_OPERATORS(op) \
HMC_OPTIONS_GEN_OPERATOR_BIN(op) \
HMC_OPTIONS_GEN_OPERATOR_EQ(op)

// Generate code for a unary operator, namely ~ and !.
// Since they have different return types, provide that as well.
#define HMC_OPTIONS_GEN_OPERATOR_UNARY(op, result_type) \
inline result_type operator op (Options opt) \
{ \
    return static_cast<result_type>(op static_cast<OPTIONS_NUMERIC>(opt)); \
}

    // Generate functions
    HMC_OPTIONS_GEN_OPERATORS(|)
    HMC_OPTIONS_GEN_OPERATORS(&)
    HMC_OPTIONS_GEN_OPERATORS(^)

    HMC_OPTIONS_GEN_OPERATOR_UNARY(~, Options)
    HMC_OPTIONS_GEN_OPERATOR_UNARY(!, bool)

    bool selected(Options a, Options b)
    {
        return !!(a & b);
    }

// Clean up macro definitions to avoid unnecessary pollution.
#undef OPTIONS_NUMERIC
#undef HMC_OPTIONS_GEN_OPERATOR_BIN
#undef HMC_OPTIONS_GEN_OPERATOR_EQ
#undef HMC_OPTIONS_GEN_OPERATORS
#undef HMC_OPTIONS_GEN_OPERATOR_UNARY



    // ==========
    // |  Cell  |
    // ==========

    // The Cell class represents and computes a single Voronoi cell within
    // the Voronoi diagram.

    // TODO: More comments here.

    // Forward declaration necessary for Cell
    class Diagram;

    class Cell {
    private:
        // TODO: Member variables are certainly part of the
        // implementation, so revisit this.
        Diagram* diagram_;
        ParticleIndex index_;
        Point position_;
        VoronoiPolyhedron poly_;
        bool polyComputed_;
        // TODO: Multiple target groups
        Group target_group_;

    public:
        // TODO: Are these all necessary?
        Cell() = default;
        Cell(const Cell& other) = default;
        Cell(Cell&& other) noexcept = default;
        Cell& operator=(const Cell& other) = default;
        Cell& operator=(Cell&& other) = default;
        ~Cell() = default;

        Cell(Diagram* diagram, double x, double y, double z,
             ParticleIndex index, Group target_group)
        : diagram_(diagram), index_(index), position_(x, y, z),
          poly_(), polyComputed_(false), target_group_(target_group)
        { /* Done */ }

        double computeVolume()
        {
            // TODO: Implement
            return 0;
        }

// Temporary macro.
// For some arbitrary computation function that takes in output
// iterators, creates some convenience overloads for that function
// which allow you to use a std::vector directly, or get one created
// and returned if you don't supply an output collection at all.
#define HMC_CELL_GEN_CONVENIENCE_OVERLOADS(compute, T) \
template <typename BackInsertableCollection, typename... Arguments> \
void compute(Arguments... args, BackInsertableCollection coll) \
{ \
    compute(std::back_inserter(args..., coll)); \
} \
\
template <typename... Arguments> \
std::vector<T> compute(Arguments... args) \
{ \
    std::vector<T> output; \
    compute(args..., output); \
    return output; \
}


        template <typename DoubleOutputIterator>
        void computeVertices(DoubleOutputIterator out)
        {
            // TODO: Implement
        }
        HMC_CELL_GEN_CONVENIENCE_OVERLOADS(computeVertices, double)


        template <typename ParticleIndexOutputIterator>
        void computeNeighbors(ParticleIndexOutputIterator out,
                              Options options = Options::NONE)
        {
            // TODO: Implement
        }
        HMC_CELL_GEN_CONVENIENCE_OVERLOADS(computeNeighbors, ParticleIndex)


        template <typename ParticleIndexOutputIterator>
        void computeNeighborCloud(double distance,
                                  ParticleIndexOutputIterator out,
                                  Options options = Options::NONE)
        {
            // TODO: Implement
        }
        HMC_CELL_GEN_CONVENIENCE_OVERLOADS(computeNeighborCloud, ParticleIndex)

#undef HMC_CELL_GEN_CONVENIENCE_OVERLOADS
    };



    // =============
    // |  Diagram  |
    // =============

    // TODO: Comments
    class Diagram {
    private:

    public:
        // TODO: Rename
        void initializeDiagram()
        {
            // TODO: Implement
        }

        void addParticle(double x, double y, double z, Id id, Group group = NO_GROUP)
        {
            // TODO: Implement
        }

        Cell getCell(double x, double y, double z, Group target_group = NO_GROUP)
        {
            // TODO: Implement correctly
            return Cell {this, x, y, z, INVALID_INDEX, target_group};
        }

        Cell getCell(ParticleIndex index, Group target_group = NO_GROUP)
        {
            // TODO: Implement correctly
            return Cell {this, 0, 0, 0, index, target_group};
        }

        // TODO: Rename
        // TODO: Multiple source groups
        // TODO: Return type (type of source groups in general)
        std::vector<ParticleIndex> sourceGroups(Group source_group)
        {
            // TODO: Implement
            return {};
        }

        // TODO: Rename
        // TODO: Multiple target groups
        // TODO: Return type (type of target groups in general)
        Group targetGroups(Group target_group)
        {
            // TODO: Implement
            return target_group;
        }
    };
}
