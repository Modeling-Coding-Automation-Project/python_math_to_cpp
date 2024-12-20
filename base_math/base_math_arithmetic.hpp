#ifndef BASE_MATH_ARITHMETIC_HPP
#define BASE_MATH_ARITHMETIC_HPP

#include "base_math_macros.hpp"

#include "base_math_mathematical_constants.hpp"

#include <cstddef>

#ifdef BASE_MATH_USE_STD_MATH
#include <cmath>
#else  // BASE_MATH_USE_STD_MATH
#endif // BASE_MATH_USE_STD_MATH

namespace Base {
namespace Math {

constexpr double ARITHMETIC_DIVISION_MIN = 1.0e-10;

constexpr std::size_t MOD_REPEAT_MAX_NUMBER = 100;

/* abs */
template <typename T> inline T abs_base_math(const T &x) {
  return x < static_cast<T>(0) ? -x : x;
}

template <typename T> inline T abs(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::abs(x);
#else  // BASE_MATH_USE_STD_MATH
  return Base::Math::abs_base_math(x);
#endif // BASE_MATH_USE_STD_MATH
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

#ifdef BASE_MATH_USE_STD_MATH
  return std::fmod(x, y);
#else  // BASE_MATH_USE_STD_MATH
  return Base::Math::mod_with_casting_integer(
      x, y, static_cast<T>(Base::Math::ARITHMETIC_DIVISION_MIN));
#endif // BASE_MATH_USE_STD_MATH
}

} // namespace Math
} // namespace Base

#endif // BASE_MATH_ARITHMETIC_HPP
