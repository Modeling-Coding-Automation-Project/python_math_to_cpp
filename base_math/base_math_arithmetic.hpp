#ifndef __BASE_MATH_ARITHMETIC_HPP__
#define __BASE_MATH_ARITHMETIC_HPP__

#include "base_math_macros.hpp"

#include "base_math_mathematical_constants.hpp"

#include <cstddef>

#ifdef __BASE_MATH_USE_STD_MATH__
#include <cmath>
#else  // __BASE_MATH_USE_STD_MATH__
#endif // __BASE_MATH_USE_STD_MATH__

namespace Base {
namespace Math {

constexpr double ARITHMETIC_DIVISION_MIN = 1.0e-10;

constexpr std::size_t MOD_REPEAT_MAX_NUMBER = 100;

/* abs */
template <typename T> inline T abs_base_math(const T &x) {
  return x < static_cast<T>(0) ? -x : x;
}

template <typename T> inline T abs(const T &x) {

#ifdef __BASE_MATH_USE_STD_MATH__
  return std::abs(x);
#else  // __BASE_MATH_USE_STD_MATH__
  return Base::Math::abs_base_math(x);
#endif // __BASE_MATH_USE_STD_MATH__
}

/* mod */
template <typename T>
inline T mod_with_casting_integer(const T &x, const T &y,
                                  const T &division_min) {
  if (y < division_min) {
    return static_cast<T>(0);
  }

  T abs_y = Base::Math::abs_base_math(y);
  T quotient = static_cast<T>(static_cast<int>(x / abs_y));
  T result = x - quotient * abs_y;

  return result;
}

template <typename T> inline T mod(const T &x, const T &y) {

#ifdef __BASE_MATH_USE_STD_MATH__
  return std::fmod(x, y);
#else  // __BASE_MATH_USE_STD_MATH__
  return Base::Math::mod_with_casting_integer(
      x, y, static_cast<T>(Base::Math::ARITHMETIC_DIVISION_MIN));
#endif // __BASE_MATH_USE_STD_MATH__
}

} // namespace Math
} // namespace Base

#endif // __BASE_MATH_ARITHMETIC_HPP__
