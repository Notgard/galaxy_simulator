#pragma once

#include <cmath> // For std::sqrt

namespace quadtree
{

template<typename T>
class Vector2
{
public:
    T x;
    T y;

    constexpr Vector2<T>(T X = 0, T Y = 0) noexcept : x(X), y(Y)
    {

    }

    constexpr Vector2<T>& operator+=(const Vector2<T>& other) noexcept
    {
        x += other.x;
        y += other.y;
        return *this;
    }

    constexpr Vector2<T>& operator/=(T t) noexcept
    {
        x /= t;
        y /= t;
        return *this;
    }

    constexpr Vector2<T>& operator*=(T t) noexcept
    {
        x *= t;
        y *= t;
        return *this;
    }

    T length() const noexcept
    {
        return std::sqrt(x * x + y * y);
    }

    Vector2<T> direction() const noexcept
    {
        T len = length();
        if (len == 0) return Vector2<T>(0, 0); // Avoid division by zero
        return Vector2<T>(x / len, y / len);
    }
};

template<typename T>
constexpr Vector2<T> operator+(Vector2<T> lhs, const Vector2<T>& rhs) noexcept
{
    lhs += rhs;
    return lhs;
}

template<typename T>
constexpr Vector2<T> operator/(Vector2<T> vec, T t) noexcept
{
    vec /= t;
    return vec;
}

template<typename T>
constexpr Vector2<T> operator*(Vector2<T> vec, T t) noexcept
{
    vec *= t;
    return vec;
}

}