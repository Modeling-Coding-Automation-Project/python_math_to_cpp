#ifndef BASE_MATH_TRIGONOMETRIC_HPP
#define BASE_MATH_TRIGONOMETRIC_HPP

#include "base_math_arithmetic.hpp"
#include "base_math_exponential_logarithmic.hpp"
#include "base_math_macros.hpp"
#include "base_math_mathematical_constants.hpp"
#include "base_math_utility.hpp"

#include <cstddef>

#ifdef BASE_MATH_USE_STD_MATH
#include <cmath>
#else  // BASE_MATH_USE_STD_MATH
#endif // BASE_MATH_USE_STD_MATH

namespace Base {
namespace Math {

#ifdef BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS
constexpr std::size_t SIN_REPEAT_NUMBER = 5;
constexpr std::size_t COS_REPEAT_NUMBER = 6;
constexpr std::size_t ATAN_REPEAT_NUMBER = 3;
#else  // BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS
constexpr std::size_t SIN_REPEAT_NUMBER = 8;
constexpr std::size_t COS_REPEAT_NUMBER = 9;
constexpr std::size_t ATAN_REPEAT_NUMBER = 8;
#endif // BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS

static constexpr double CHEBYSHEV_COEFFICIENT_FOR_ATAN[11] = {
    static_cast<double>(1.0),
    static_cast<double>(-0.3333333333333333),
    static_cast<double>(0.2),
    static_cast<double>(-0.14285714285714285),
    static_cast<double>(0.1111111111111111),
    static_cast<double>(-0.09090909090909091),
    static_cast<double>(0.07692307692307693),
    static_cast<double>(-0.06666666666666667),
    static_cast<double>(0.058823529411764705),
    static_cast<double>(-0.05263157894736842),
    static_cast<double>(0.047619047619047616)};

constexpr double TRIGONOMETRIC_DIVISION_MIN = 1e-10;

/* wrap */
template <typename T> inline T wrap_value_in_minus_pi_and_pi(const T &x) {
  T result = static_cast<T>(0);

  if (x >= static_cast<T>(0)) {
    result = Base::Math::mod_base_math(x + static_cast<T>(Base::Math::PI),
                                       static_cast<T>(2) *
                                           static_cast<T>(Base::Math::PI)) -
             static_cast<T>(Base::Math::PI);
  } else {
    result = Base::Math::mod_base_math(x - static_cast<T>(Base::Math::PI),
                                       static_cast<T>(2) *
                                           static_cast<T>(Base::Math::PI)) +
             static_cast<T>(Base::Math::PI);
  }

  return result;
}

/* sin */
template <typename T, std::size_t LOOP_NUMBER>
inline T sin_base_math(const T &x) {

  T x_wrapped = Base::Math::wrap_value_in_minus_pi_and_pi(x);

  T term = x_wrapped;
  T result = x_wrapped;
  T x_squared = x_wrapped * x_wrapped;

  for (std::size_t n = 1; n < LOOP_NUMBER; n++) {
    term *= -x_squared / static_cast<T>((2 * n) * (2 * n + 1));
    result += term;
  }

  return result;
}

template <typename T> inline T sin(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::sin(x);
#else
  return Base::Math::sin_base_math<T, Base::Math::SIN_REPEAT_NUMBER>(x);
#endif
}

/* cos */
template <typename T, std::size_t LOOP_NUMBER>
inline T cos_base_math(const T &x) {

  T x_wrapped = wrap_value_in_minus_pi_and_pi(x);

  T term = static_cast<T>(1);
  T result = static_cast<T>(1);
  T x_squared = x_wrapped * x_wrapped;

  for (std::size_t n = 1; n < LOOP_NUMBER; n++) {
    term *= -x_squared / static_cast<T>((2 * n - 1) * (2 * n));
    result += term;
  }

  return result;
}

template <typename T> inline T cos(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::cos(x);
#else
  return Base::Math::cos_base_math<T, Base::Math::COS_REPEAT_NUMBER>(x);
#endif
}

/* tan */
template <typename T, std::size_t SIN_LOOP_NUMBER, std::size_t COS_LOOP_NUMBER>
inline T tan_base_math(const T &x) {

  return Base::Math::sin_base_math<T, SIN_LOOP_NUMBER>(x) /
         Base::Math::avoid_zero_divide(
             Base::Math::cos_base_math<T, COS_LOOP_NUMBER>(x),
             static_cast<T>(Base::Math::TRIGONOMETRIC_DIVISION_MIN));
}

template <typename T> inline T tan(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::tan(x);
#else
  return Base::Math::tan_base_math<T, Base::Math::SIN_REPEAT_NUMBER,
                                   Base::Math::COS_REPEAT_NUMBER>(x);
#endif
}

/* atan */
template <typename T, std::size_t LOOP_NUMBER>
inline T atan_base_math(const T &x) {

  T result = static_cast<T>(0);

  if (x > static_cast<T>(1)) {
    result = static_cast<T>(Base::Math::HALF_PI) -
             Base::Math::atan_base_math(static_cast<T>(1) / x);
  } else if (x < static_cast<T>(-1)) {
    result = -static_cast<T>(Base::Math::HALF_PI) -
             Base::Math::atan_base_math(static_cast<T>(1) / x);
  } else {
    if ((x > static_cast<T>(0.5)) || (x < static_cast<T>(-0.5))) {
      T half_x = x * static_cast<T>(0.5);

      result = Base::Math::atan_base_math(half_x) +
               Base::Math::atan_base_math(
                   half_x /
                   (static_cast<T>(1) + static_cast<T>(2) * half_x * half_x));

    } else {
      T term = x;
      T x_squared = x * x;

      result = x;

      for (std::size_t n = 1; n < LOOP_NUMBER; n++) {
        term *=
            -x_squared * static_cast<T>(2 * n - 1) / static_cast<T>(2 * n + 1);
        result += term;
      }
    }
  }

  return result;
}

template <typename T, int N> struct AtanChebyshevLoop {
  static T compute(const T &x_squared, const T &y) {
    return AtanChebyshevLoop<T, N - 1>::compute(
        x_squared,
        y * x_squared +
            static_cast<T>(Base::Math::CHEBYSHEV_COEFFICIENT_FOR_ATAN[N]));
  }
};

template <typename T> struct AtanChebyshevLoop<T, 0> {
  static T compute(const T &x_squared, const T &y) {
    return y * x_squared +
           static_cast<T>(Base::Math::CHEBYSHEV_COEFFICIENT_FOR_ATAN[0]);
  }
};

template <typename T, std::size_t LOOP_NUMBER>
inline T atan_chebyshev_core(const T &x) {

  T x_squared = x * x;
  T y = static_cast<T>(0);

  y = AtanChebyshevLoop<T, LOOP_NUMBER - 1>::compute(x_squared, y);

  return y * x;
}

template <typename T, std::size_t LOOP_NUMBER>
inline T atan_chebyshev_core_wide(const T &x) {

  if ((x > static_cast<T>(0.5)) || (x < static_cast<T>(-0.5))) {
    T half_x = x * static_cast<T>(0.5);

    return Base::Math::atan_chebyshev_core<T, LOOP_NUMBER>(half_x) +
           Base::Math::atan_chebyshev_core<T, LOOP_NUMBER>(
               half_x /
               (static_cast<T>(1) + static_cast<T>(2) * half_x * half_x));

  } else {

    return Base::Math::atan_chebyshev_core<T, LOOP_NUMBER>(x);
  }
}

template <typename T, std::size_t LOOP_NUMBER>
inline T atan_chebyshev(const T &x) {

  T result = static_cast<T>(0);

  if (x > static_cast<T>(1)) {
    result = static_cast<T>(Base::Math::HALF_PI) -
             Base::Math::atan_chebyshev_core_wide<T, LOOP_NUMBER>(
                 static_cast<T>(1) / x);
  } else if (x < static_cast<T>(-1)) {
    result = -static_cast<T>(Base::Math::HALF_PI) -
             Base::Math::atan_chebyshev_core_wide<T, LOOP_NUMBER>(
                 static_cast<T>(1) / x);
  } else {

    result = Base::Math::atan_chebyshev_core_wide<T, LOOP_NUMBER>(x);
  }

  return result;
}

template <typename T> inline T atan(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::atan(x);
#else
#ifdef BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS
  return Base::Math::atan_chebyshev<T, Base::Math::ATAN_REPEAT_NUMBER>(x);
#else  // BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS
  return Base::Math::atan_chebyshev<T, Base::Math::ATAN_REPEAT_NUMBER>(x);
#endif // BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS
#endif
}

/* atan2 */
template <typename T, std::size_t LOOP_NUMBER>
inline T atan2_base_math(const T &y, const T &x) {
  T result = static_cast<T>(0);

  if (Base::Math::near_zero(
          x, static_cast<T>(Base::Math::TRIGONOMETRIC_DIVISION_MIN)) &&
      Base::Math::near_zero(
          y, static_cast<T>(Base::Math::TRIGONOMETRIC_DIVISION_MIN))) {

    return static_cast<T>(0);

  } else if (x > static_cast<T>(0)) {

    result = Base::Math::atan_chebyshev<T, LOOP_NUMBER>(y / x);

  } else if (x < static_cast<T>(0) && y >= static_cast<T>(0)) {

    result = Base::Math::atan_chebyshev<T, LOOP_NUMBER>(y / x) +
             static_cast<T>(Base::Math::PI);

  } else if (x < static_cast<T>(0) && y < static_cast<T>(0)) {

    result = Base::Math::atan_chebyshev<T, LOOP_NUMBER>(y / x) -
             static_cast<T>(Base::Math::PI);

  } else if (Base::Math::near_zero(
                 x, static_cast<T>(Base::Math::TRIGONOMETRIC_DIVISION_MIN)) &&
             y > static_cast<T>(0)) {

    result = static_cast<T>(Base::Math::HALF_PI);

  } else if (Base::Math::near_zero(
                 x, static_cast<T>(Base::Math::TRIGONOMETRIC_DIVISION_MIN)) &&
             y < static_cast<T>(0)) {

    result = -static_cast<T>(Base::Math::HALF_PI);
  }

  return result;
}

template <typename T> inline T atan2(const T &y, const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::atan2(y, x);
#else
  return Base::Math::atan2_base_math<T, Base::Math::ATAN_REPEAT_NUMBER>(y, x);
#endif
}

/* asin */
template <typename T, std::size_t ATAN_LOOP_NUMBER,
          std::size_t SQRT_LOOP_NUMBER>
inline T asin_base_math(const T &x) {

  T result = static_cast<T>(0);

  if (x >= static_cast<T>(1)) {

    result = static_cast<T>(Base::Math::HALF_PI);

  } else if (x <= static_cast<T>(-1)) {

    result = -static_cast<T>(Base::Math::HALF_PI);

  } else {

    result = static_cast<T>(2) *
             Base::Math::atan2_base_math<T, ATAN_LOOP_NUMBER>(
                 x, static_cast<T>(1) +
                        Base::Math::sqrt_base_math<T, SQRT_LOOP_NUMBER>(
                            static_cast<T>(1) - x * x));
  }

  return result;
}

template <typename T> inline T asin(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::asin(x);
#else
  return Base::Math::asin_base_math<T, Base::Math::ATAN_REPEAT_NUMBER,
                                    Base::Math::SQRT_REPEAT_NUMBER>(x);
#endif
}

/* acos */
template <typename T, std::size_t ATAN_LOOP_NUMBER,
          std::size_t SQRT_LOOP_NUMBER>
inline T acos_base_math(const T &x) {

  return static_cast<T>(Base::Math::HALF_PI) -
         Base::Math::asin_base_math<T, ATAN_LOOP_NUMBER, SQRT_LOOP_NUMBER>(x);
}

template <typename T> inline T acos(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::acos(x);
#else
  return Base::Math::acos_base_math<T, Base::Math::ATAN_REPEAT_NUMBER,
                                    Base::Math::SQRT_REPEAT_NUMBER>(x);
#endif
}

/* sinh */
template <typename T, std::size_t LOOP_NUMBER>
inline T sinh_base_math(const T &x) {
  return (Base::Math::exp_base_math<T, LOOP_NUMBER>(x) -
          Base::Math::exp_base_math<T, LOOP_NUMBER>(-x)) *
         static_cast<T>(0.5);
}

template <typename T> inline T sinh(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::sinh(x);
#else
  return Base::Math::sinh_base_math<T, Base::Math::EXP_REPEAT_NUMBER>(x);
#endif
}

/* cosh */
template <typename T, std::size_t LOOP_NUMBER>
inline T cosh_base_math(const T &x) {
  return (Base::Math::exp_base_math<T, LOOP_NUMBER>(x) +
          Base::Math::exp_base_math<T, LOOP_NUMBER>(-x)) *
         static_cast<T>(0.5);
}

template <typename T> inline T cosh(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::cosh(x);
#else
  return Base::Math::cosh_base_math<T, Base::Math::EXP_REPEAT_NUMBER>(x);
#endif
}

/* tanh */
template <typename T, std::size_t LOOP_NUMBER>
inline T tanh_base_math(const T &x) {

  T a = Base::Math::exp_base_math<T, LOOP_NUMBER>(x);
  T b = Base::Math::exp_base_math<T, LOOP_NUMBER>(-x);

  return (a - b) /
         avoid_zero_divide(a + b, static_cast<T>(TRIGONOMETRIC_DIVISION_MIN));
}

template <typename T> inline T tanh(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::tanh(x);
#else
  return Base::Math::tanh_base_math<T, Base::Math::EXP_REPEAT_NUMBER>(x);
#endif
}

} // namespace Math
} // namespace Base

#endif // BASE_MATH_TRIGONOMETRIC_HPP
