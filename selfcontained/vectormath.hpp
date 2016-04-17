#pragma once

#include <cmath>
#include <iostream>

namespace hmc {
    static constexpr double VECTORMATH_TOLERANCE = 1e-12;

    constexpr bool approxEq(const double a, const double b)
    {
        return std::abs(a - b) <= VECTORMATH_TOLERANCE;
    }

    constexpr bool approxEq(const double a, const double b, const double tolerance)
    {
        return std::abs(a - b) <= tolerance;
    }

    struct Vector3 {
        double x;
        double y;
        double z;

        Vector3() = default;

        constexpr Vector3(double x, double y, double z)
            : x(x), y(y), z(z)
        { /* Done */ }

        friend std::ostream& operator<<(std::ostream& out, const Vector3& v)
        {
            return out << "(" << v.x << ", " << v.y << ", " << v.z << ")";
        }

        Vector3& operator+=(const Vector3& v)
        {
            x += v.x;
            y += v.y;
            z += v.z;
            return *this;
        }

        Vector3& operator-=(const Vector3& v)
        {
            x -= v.x;
            y -= v.y;
            z -= v.z;
            return *this;
        }

        Vector3& operator*=(const double a)
        {
            x *= a;
            y *= a;
            z *= a;
            return *this;
        }

        Vector3& operator/=(const double a)
        {
            x /= a;
            y /= a;
            z /= a;
            return *this;
        }

        Vector3 operator-() const
        {
            return Vector3(-x, -y, -z);
        }


        friend constexpr Vector3 operator+(const Vector3& a, const Vector3& b)
        {
            return Vector3(a.x + b.x, a.y + b.y, a.z + b.z);
        }

        friend constexpr Vector3 operator-(const Vector3& a, const Vector3& b)
        {
            return Vector3(a.x - b.x, a.y - b.y, a.z - b.z);
        }


        friend constexpr Vector3 operator*(const Vector3& v, const double a)
        {
            return Vector3(v.x * a, v.y * a, v.z * a);
        }
        friend constexpr Vector3 operator*(const double a, const Vector3& v)
        {
            return Vector3(v.x * a, v.y * a, v.z * a);
        }

        friend constexpr Vector3 operator/(const Vector3& v, const double a)
        {
            // TODO: Consider computing 1/a, and then multiplying.
            // I believe this would gain speed at the cost of accuracy.
            // Remember to also update /= if any changes are made.
            return Vector3(v.x / a, v.y / a, v.z / a);
        }

        // Note: `normalize` produces a unit vector IN PLACE.
        // That is, it MODIFIES the vector it is called on.
        // Thus to avoid confusion it is void.
        friend double mag(const Vector3&);
        void normalize()
        {
            *this /= mag(*this);
        }
    };

    constexpr Vector3 VECTOR_ZERO {0, 0, 0};

    // Computes the dot product of two vectors
    constexpr double dot(const Vector3& a, const Vector3& b)
    {
        return a.x * b.x   +   a.y * b.y   +   a.z * b.z;
    }

    // Computes the squared magnitude (squared length) of a vector
    constexpr double mag2(const Vector3& v)
    {
        return dot(v, v);
    }

    // Computes the magnitude (length) of a vector
    // Can't be constexpr because sqrt isn't! >:(
    // TODO: A constexpr sqrt would allow this and the
    //       functions below to be constexpr.
    //       See also some functions in Plane.
    double mag(const Vector3& v)
    {
        return std::sqrt(mag2(v));
    }

    // Computes a unit vector in the direction of v
    Vector3 unit(const Vector3& v)
    {
        // TODO: Consider using a fast inverse sqrt and then
        // multiplying. This would gain speed at the cost of
        // accuracy. Potentially have a 'approximateUnit'
        // function to accomodate this.
        return v / mag(v);
    }

    // Computes the cross product
    constexpr Vector3 cross(const Vector3& a, const Vector3& b)
    {
        return Vector3(
            a.y * b.z   -   a.z * b.y,
            a.z * b.x   -   a.x * b.z,
            a.x * b.y   -   a.y * b.x
        );
    }

    // Computes a unit vector in the direction of the cross product
    Vector3 unitCross(const Vector3& a, const Vector3& b)
    {
        return unit(cross(a, b));
    }

    // Computes the squared distance between two vectors
    constexpr double squaredDistance(const Vector3& a, const Vector3& b)
    {
        return mag2(a - b);
    }

    // Computes the distance between two vectors
    double distance(const Vector3& a, const Vector3& b)
    {
        return mag(a - b);
    }

    // Computes the midpoint (halfway point) between two vectors
    constexpr Vector3 midpoint(const Vector3& a, const Vector3& b)
    {
        return (a + b) * 0.5;
    }


    struct Plane {
        Vector3 unitNormal;
        double planeOffset;

        enum Location {
            OUTSIDE,
            INCIDENT,
            INSIDE
        };

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

        constexpr Vector3 intersection(const Vector3& a, const Vector3& b) const
        {
            return a + (b - a) * (planeOffset - dot(a, unitNormal)) / dot(b - a, unitNormal);
        }

        constexpr Location location(const Vector3& p) const
        {
            return location(signedDistance(p));
        }

        constexpr Location location(const double signedDistance) const
        {
            return signedDistance >  VECTORMATH_TOLERANCE ? OUTSIDE
                :  signedDistance < -VECTORMATH_TOLERANCE ? INSIDE
                :                                           INCIDENT;
        }
    };
}
