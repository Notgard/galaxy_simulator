#pragma once

#include <cmath> // For std::sqrt
#include <ostream>

namespace quadtree
{
    template <typename T>
    class Vector2
    {
    public:
        T x;
        T y;

        constexpr Vector2(T X = 0, T Y = 0) noexcept : x(X), y(Y) {}

        constexpr Vector2<T> &operator=(const Vector2<T> &other) noexcept
        {
            if (this != &other)
            {
                x = other.x;
                y = other.y;
            }
            return *this;
        }

        constexpr Vector2<T> &operator+=(const Vector2<T> &other) noexcept
        {
            x += other.x;
            y += other.y;
            return *this;
        }

        constexpr Vector2<T> &operator-=(const Vector2<T> &other) noexcept
        {
            x -= other.x;
            y -= other.y;
            return *this;
        }

        constexpr Vector2<T> &operator/=(T t) noexcept
        {
            x /= t;
            y /= t;
            return *this;
        }

        constexpr Vector2<T> &operator*=(T t) noexcept
        {
            x *= t;
            y *= t;
            return *this;
        }

        constexpr Vector2<T> operator-() const noexcept
        {
            return Vector2<T>(-x, -y);
        }

        constexpr bool operator==(const Vector2<T> &other) const noexcept
        {
            return x == other.x && y == other.y;
        }

        constexpr bool operator!=(const Vector2<T> &other) const noexcept
        {
            return !(*this == other);
        }

        static T length(const Vector2<T> &c) noexcept
        {
            return std::sqrt(c.x * c.x + c.y * c.y);
        }

        T length() const noexcept
        {
            return std::sqrt(x * x + y * y);
        }

        Vector2<T> normalized() const noexcept
        {
            T len = length();
            return (len == 0) ? Vector2<T>(0, 0) : Vector2<T>(x / len, y / len);
        }

        static T dot(const Vector2<T> &a, const Vector2<T> &b) noexcept
        {
            return a.x * b.x + a.y * b.y;
        }

        static T cross(const Vector2<T> &a, const Vector2<T> &b) noexcept
        {
            return a.x * b.y - a.y * b.x;
        }

        static T distance(const Vector2<T> &c1, const Vector2<T> &c2) noexcept
        {
            return (c1 - c2).length();
        }

        friend std::ostream &operator<<(std::ostream &os, const Vector2<T> &vec)
        {
            return os << "(" << vec.x << ", " << vec.y << ")";
        }
    };

    template <typename T>
    constexpr Vector2<T> operator+(Vector2<T> lhs, const Vector2<T> &rhs) noexcept
    {
        lhs += rhs;
        return lhs;
    }

    template <typename T>
    constexpr Vector2<T> operator-(Vector2<T> lhs, const Vector2<T> &rhs) noexcept
    {
        lhs -= rhs;
        return lhs;
    }

    template <typename T>
    constexpr Vector2<T> operator/(Vector2<T> vec, T t) noexcept
    {
        vec /= t;
        return vec;
    }

    template <typename T>
    constexpr Vector2<T> operator*(Vector2<T> vec, T t) noexcept
    {
        vec *= t;
        return vec;
    }

    template <typename T>
    constexpr Vector2<T> operator*(T t, Vector2<T> vec) noexcept
    {
        return vec * t;
    }
}
