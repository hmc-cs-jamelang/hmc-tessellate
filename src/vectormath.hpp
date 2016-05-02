#pragma once

#include <cmath>
#include <iostream>

namespace hmc {
    static constexpr double VECTORMATH_TOLERANCE = 1e-12;

    inline constexpr bool approxEq(double a, double b, double tolerance = VECTORMATH_TOLERANCE)
    {
        return std::abs(a - b) <= tolerance;
    }

    inline constexpr bool approxRelEq(double a, double b, double tolerance = VECTORMATH_TOLERANCE)
    {
        return (a == 0 || b == 0)
                    ? (a == b)
                    : approxEq(a / b, 1, tolerance);
    }

    template <typename Real>
    struct Vector3_of {
        Real x;
        Real y;
        Real z;

        Vector3_of() = default;

        inline constexpr Vector3_of(double x, double y, double z)
            : x(x), y(y), z(z)
        { /* Done */ }

        inline friend std::ostream& operator<<(std::ostream& out, const Vector3_of& v)
        {
            return out << "(" << v.x << ", " << v.y << ", " << v.z << ")";
        }

        Vector3_of& operator+=(const Vector3_of& v)
        {
            x += v.x;
            y += v.y;
            z += v.z;
            return *this;
        }

        Vector3_of& operator-=(const Vector3_of& v)
        {
            x -= v.x;
            y -= v.y;
            z -= v.z;
            return *this;
        }

        Vector3_of& operator*=(const double a)
        {
            x *= a;
            y *= a;
            z *= a;
            return *this;
        }

        Vector3_of& operator/=(const double a)
        {
            x /= a;
            y /= a;
            z /= a;
            return *this;
        }

        Vector3_of operator-() const
        {
            return Vector3_of(-x, -y, -z);
        }


        inline friend constexpr Vector3_of operator+(const Vector3_of& a, const Vector3_of& b)
        {
            return Vector3_of(a.x + b.x, a.y + b.y, a.z + b.z);
        }

        inline friend constexpr Vector3_of operator-(const Vector3_of& a, const Vector3_of& b)
        {
            return Vector3_of(a.x - b.x, a.y - b.y, a.z - b.z);
        }


        inline friend constexpr Vector3_of operator*(const Vector3_of& v, const double a)
        {
            return Vector3_of(v.x * a, v.y * a, v.z * a);
        }

        inline friend constexpr Vector3_of operator*(const double a, const Vector3_of& v)
        {
            return Vector3_of(v.x * a, v.y * a, v.z * a);
        }

        inline friend constexpr Vector3_of operator/(const Vector3_of& v, const double a)
        {
            // TODO: Consider computing 1/a, and then multiplying.
            // I believe this would gain speed at the cost of accuracy.
            // Remember to also update /= if any changes are made.
            return Vector3_of(v.x / a, v.y / a, v.z / a);
        }

        // Note: `normalize` produces a unit vector IN PLACE.
        // That is, it MODIFIES the vector it is called on.
        // Thus to avoid confusion it is void.
        void normalize()
        {
            *this /= mag(*this);
        }

        // Computes the dot product of two vectors
        inline friend constexpr double dot(const Vector3_of& a, const Vector3_of& b)
        {
            return a.x * b.x   +   a.y * b.y   +   a.z * b.z;
        }

        // Computes the squared magnitude (squared length) of a vector
        inline friend constexpr double mag2(const Vector3_of& v)
        {
            return dot(v, v);
        }

        // Computes the magnitude (length) of a vector
        // Can't be constexpr because sqrt isn't! >:(
        // TODO: A constexpr sqrt would allow this and the
        //       functions below to be constexpr.
        //       See also some functions in Plane.
        inline friend double mag(const Vector3_of& v)
        {
            return std::sqrt(mag2(v));
        }

        // Computes a unit vector in the direction of v
        inline friend Vector3_of unit(const Vector3_of& v)
        {
            // TODO: Consider using a fast inverse sqrt and then
            // multiplying. This would gain speed at the cost of
            // accuracy. Potentially have a 'approximateUnit'
            // function to accomodate this.
            return v / mag(v);
        }

        // Computes the cross product
        inline friend constexpr Vector3_of cross(const Vector3_of& a, const Vector3_of& b)
        {
            return Vector3_of(
                a.y * b.z   -   a.z * b.y,
                a.z * b.x   -   a.x * b.z,
                a.x * b.y   -   a.y * b.x
            );
        }

        // Computes a unit vector in the direction of the cross product
        inline friend Vector3_of unitCross(const Vector3_of& a, const Vector3_of& b)
        {
            return unit(cross(a, b));
        }

        // Computes the squared distance between two vectors
        inline friend constexpr double squaredDistance(const Vector3_of& a, const Vector3_of& b)
        {
            return mag2(a - b);
        }

        // Computes the distance between two vectors
        inline friend double distance(const Vector3_of& a, const Vector3_of& b)
        {
            return mag(a - b);
        }

        // Computes the midpoint (halfway point) between two vectors
        inline friend constexpr Vector3_of midpoint(const Vector3_of& a, const Vector3_of& b)
        {
            return (a + b) * 0.5;
        }
    };

    using Vector3 = Vector3_of<double>;

    constexpr Vector3 VECTOR_ZERO {0, 0, 0};

    constexpr Vector3 VECTOR_MAX {
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max()
    };

    constexpr Vector3 VECTOR_MIN {
        std::numeric_limits<double>::lowest(),
        std::numeric_limits<double>::lowest(),
        std::numeric_limits<double>::lowest()
    };



    struct Plane {
        Vector3 unitNormal;
        double planeOffset;

        enum Location {
            OUTSIDE,
            INCIDENT,
            INSIDE
        };

        inline friend std::ostream& operator<<(std::ostream& out, Location loc)
        {
            switch (loc) {
                case OUTSIDE: {
                    return out << "OUTSIDE";
                } break;
                case INCIDENT: {
                    return out << "INCIDENT";
                } break;
                case INSIDE: {
                    return out << "INSIDE";
                } break;
                default: {
                    return out << "<Unknown Plane::Location>";
                } break;
            }
        }

        constexpr Plane(const Vector3& unitNormal,  const double offset)
            : unitNormal(unitNormal), planeOffset(offset)
        { /* Done */ }

        static constexpr Plane normalAndPoint(const Vector3& unitNormal,
                                              const Vector3& point)
        {
            return Plane(unitNormal, dot(unitNormal, point));
        }

        static Plane halfwayFromOriginTo(const Vector3& v)
        {
            return Plane::normalAndPoint(unit(v), v/2);
        }

        // Not constexpr because unit is not constexpr.
        static Plane between(const Vector3& a, const Vector3& b)
        {
            return Plane::normalAndPoint(unit(b - a), midpoint(a, b));
        }

        constexpr double offset(const Vector3& p) const
        {
            return dot(unitNormal, p);
        }

        constexpr double signedDistance(const Vector3& p) const
        {
            return offset(p) - planeOffset;
        }

        Vector3 intersection(const Vector3& a, const Vector3& b) const
        {
            double aOffset = offset(a);
            double bOffset = offset(b);
            return a + (b - a) * (planeOffset - aOffset) / (bOffset - aOffset);
        }

        constexpr Location location(const Vector3& p, double tolerance = VECTORMATH_TOLERANCE) const
        {
            return location(signedDistance(p), tolerance);
        }

        constexpr Location location(const double signedDistance, double tolerance = VECTORMATH_TOLERANCE) const
        {
            return signedDistance >  tolerance ? OUTSIDE
                :  signedDistance < -tolerance ? INSIDE
                :                                INCIDENT;
        }

        inline friend std::ostream& operator<<(std::ostream& out, const Plane& plane)
        {
            return out << "Plane {normal = " << plane.unitNormal
                        << ", offset = " << plane.planeOffset
                        << "}";
        }
    };


    struct BoundingBox {
        Vector3 low = VECTOR_MAX;
        Vector3 high = VECTOR_MIN;

        void reset()
        {
            low = VECTOR_MAX;
            high = VECTOR_MIN;
        }

        void adjustToContain(Vector3 newPoint)
        {
            if (newPoint.x < low.x) { low.x = newPoint.x; }
            if (newPoint.y < low.y) { low.y = newPoint.y; }
            if (newPoint.z < low.z) { low.z = newPoint.z; }

            if (newPoint.x > high.x) { high.x = newPoint.x; }
            if (newPoint.y > high.y) { high.y = newPoint.y; }
            if (newPoint.z > high.z) { high.z = newPoint.z; }
        }

        void pad(double padding = 1e-6)
        {
            low.x -= padding;
            low.y -= padding;
            low.z -= padding;

            high.x += padding;
            high.y += padding;
            high.z += padding;
        }

        inline friend std::ostream&
        operator<<(std::ostream& out, const BoundingBox& bb)
        {
            return out << "Bounding box "
                        << bb.low << " by " << bb.high;
        }
    };
}
