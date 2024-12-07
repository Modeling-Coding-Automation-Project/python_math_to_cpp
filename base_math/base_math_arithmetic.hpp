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
template <typename T> inline T mod_base_math(const T &x, const T &y) {
  T result = static_cast<T>(0);
  T abs_y = Base::Math::abs_base_math(y);

  if (x > static_cast<T>(0)) {
    T dif = x - abs_y;

    if (dif <= static_cast<T>(0)) {
      result = x;
    } else {
      for (size_t i = 0; i < Base::Math::MOD_REPEAT_MAX_NUMBER; i++) {
        if (dif < abs_y) {
          result = dif;
          break;
        } else {
          dif -= abs_y;
        }
      }
    }
  } else {
    T dif = -x - abs_y;

    if (dif <= static_cast<T>(0)) {
      result = x;
    } else {
      for (size_t i = 0; i < Base::Math::MOD_REPEAT_MAX_NUMBER; i++) {
        if (dif < abs_y) {
          result = -dif;
          break;
        } else {
          dif -= abs_y;
        }
      }
    }
  }

  return result;
}

template <typename T> inline T mod(const T &x, const T &y) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::fmod(x, y);
#else  // BASE_MATH_USE_STD_MATH
  return Base::Math::mod_base_math(x, y);
#endif // BASE_MATH_USE_STD_MATH
}

} // namespace Math
} // namespace Base

#endif // BASE_MATH_ARITHMETIC_HPP
