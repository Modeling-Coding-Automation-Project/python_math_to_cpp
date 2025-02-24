#ifndef __BASE_MATH_TRIGONOMETRIC_HPP__
#define __BASE_MATH_TRIGONOMETRIC_HPP__

#include "base_math_macros.hpp"

#include "base_math_arithmetic.hpp"
#include "base_math_exponential_logarithmic.hpp"
#include "base_math_mathematical_constants.hpp"
#include "base_utility.hpp"

#include <cstddef>

#ifdef __BASE_MATH_USE_STD_MATH__
#include <cmath>
#else  // __BASE_MATH_USE_STD_MATH__
#endif // __BASE_MATH_USE_STD_MATH__

namespace Base {
namespace Math {

#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
constexpr std::size_t SIN_REPEAT_NUMBER = 5;
constexpr std::size_t COS_REPEAT_NUMBER = 6;
constexpr std::size_t ATAN_REPEAT_NUMBER = 3;

constexpr std::size_t COS_MACLAURIN_DOUBLEANGLE_REPEAT_NUMBER = 1;
constexpr std::size_t SIN_MACLAURIN_DOUBLEANGLE_REPEAT_NUMBER =
    COS_MACLAURIN_DOUBLEANGLE_REPEAT_NUMBER;

constexpr std::size_t SINCOS_MACLAURIN_DOUBLEANGLE_REPEAT_NUMBER = 2;

#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
constexpr std::size_t SIN_REPEAT_NUMBER = 8;
constexpr std::size_t COS_REPEAT_NUMBER = 9;
constexpr std::size_t ATAN_REPEAT_NUMBER = 8;

constexpr std::size_t COS_MACLAURIN_DOUBLEANGLE_REPEAT_NUMBER = 3;
constexpr std::size_t SIN_MACLAURIN_DOUBLEANGLE_REPEAT_NUMBER =
    COS_MACLAURIN_DOUBLEANGLE_REPEAT_NUMBER;

constexpr std::size_t SINCOS_MACLAURIN_DOUBLEANGLE_REPEAT_NUMBER = 3;

#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__

constexpr std::size_t COS_MACLAURIN_DOUBLEANGLE_FACTOR_MAX_SIZE = 6;

using COS_MACLAURIN_FACTOR_LIST = typename CosMaclaurinFactor::MakeList<
    COS_MACLAURIN_DOUBLEANGLE_FACTOR_MAX_SIZE>::type;

constexpr auto COS_MACLAURIN_DOUBLEANGLE_FACTOR =
    CosMaclaurinFactor::to_array(COS_MACLAURIN_FACTOR_LIST{});

constexpr std::size_t SIN_MACLAURIN_DOUBLEANGLE_FACTOR_MAX_SIZE = 6;

using SIN_MACLAURIN_FACTOR_LIST = typename SinMaclaurinFactor::MakeList<
    SIN_MACLAURIN_DOUBLEANGLE_FACTOR_MAX_SIZE>::type;

constexpr auto SIN_MACLAURIN_DOUBLEANGLE_FACTOR =
    SinMaclaurinFactor::to_array(SIN_MACLAURIN_FACTOR_LIST{});

constexpr std::size_t CHEBYSHEV_COEFFICIENT_FOR_ATAN_SIZE = 11;
static constexpr double
    CHEBYSHEV_COEFFICIENT_FOR_ATAN[CHEBYSHEV_COEFFICIENT_FOR_ATAN_SIZE] = {
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
    result = Base::Math::mod(x + static_cast<T>(Base::Math::PI),
                             static_cast<T>(Base::Math::TWO_PI)) -
             static_cast<T>(Base::Math::PI);
  } else {
    result = Base::Math::mod(x - static_cast<T>(Base::Math::PI),
                             static_cast<T>(Base::Math::TWO_PI)) +
             static_cast<T>(Base::Math::PI);
  }

  return result;
}

/* cos maclaurin expansion with DoubleAngleFormula */
template <typename T, std::size_t N,
          std::size_t MACLAURIN_EXPANSION_REPEAT_NUMBER>
struct CosMaclaurinExpansionFirstLoop {
  static void compute(T &y, T &x_wrapped) {
    CosMaclaurinExpansionFirstLoop<
        T, N - 1, MACLAURIN_EXPANSION_REPEAT_NUMBER>::compute(y, x_wrapped);
    y = y * x_wrapped +
        static_cast<T>(
            COS_MACLAURIN_DOUBLEANGLE_FACTOR[MACLAURIN_EXPANSION_REPEAT_NUMBER -
                                             1 - N]);
  }
};

template <typename T, std::size_t MACLAURIN_EXPANSION_REPEAT_NUMBER>
struct CosMaclaurinExpansionFirstLoop<T, 0, MACLAURIN_EXPANSION_REPEAT_NUMBER> {
  static void compute(T &y, T &x_wrapped) {
    y = y * x_wrapped +
        static_cast<T>(
            COS_MACLAURIN_DOUBLEANGLE_FACTOR[MACLAURIN_EXPANSION_REPEAT_NUMBER -
                                             1]);
  }
};

template <typename T, std::size_t N> struct CosMaclaurinExpansionSecondLoop {
  static void compute(T &y) {
    y = y * (static_cast<T>(4) - y);
    CosMaclaurinExpansionSecondLoop<T, N - 1>::compute(y);
  }
};

template <typename T> struct CosMaclaurinExpansionSecondLoop<T, 0> {
  static void compute(T &y) {
    static_cast<void>(y);
    /* Do Nothing. */
  }
};

template <typename T, std::size_t MACLAURIN_EXPANSION_REPEAT_NUMBER>
inline T cos_maclaurin_expansion_with_DoubleAngleFormula(const T &x) {
  static_assert(MACLAURIN_EXPANSION_REPEAT_NUMBER <
                    COS_MACLAURIN_DOUBLEANGLE_FACTOR_MAX_SIZE,
                "MACLAURIN_EXPANSION_REPEAT_NUMBER is too large.");

  T x_wrapped = Base::Math::wrap_value_in_minus_pi_and_pi(x);
  T y = static_cast<T>(0);

  x_wrapped =
      x_wrapped /
      static_cast<T>(Factorial_2<MACLAURIN_EXPANSION_REPEAT_NUMBER>::value);
  x_wrapped = x_wrapped * x_wrapped;

  y = static_cast<T>(
      COS_MACLAURIN_DOUBLEANGLE_FACTOR[MACLAURIN_EXPANSION_REPEAT_NUMBER]);

  CosMaclaurinExpansionFirstLoop<
      T, (MACLAURIN_EXPANSION_REPEAT_NUMBER - 1),
      MACLAURIN_EXPANSION_REPEAT_NUMBER>::compute(y, x_wrapped);

  y = y * x_wrapped;

  CosMaclaurinExpansionSecondLoop<
      T, MACLAURIN_EXPANSION_REPEAT_NUMBER>::compute(y);

  return static_cast<T>(1) - y * static_cast<T>(0.5);
}

/* sin cos Maclaurin expansion with DoubleAngleFormula */
template <typename T, std::size_t N, std::size_t I>
struct SinCosMcLoughlinExpansionFirstLoop {
  static void compute(T &c, T &s, const T &z) {
    SinCosMcLoughlinExpansionFirstLoop<T, N, I - 1>::compute(c, s, z);
    c = c * z + static_cast<T>(COS_MACLAURIN_DOUBLEANGLE_FACTOR[N - 1 - I]);
    s = s * z + static_cast<T>(SIN_MACLAURIN_DOUBLEANGLE_FACTOR[N - 1 - I]);
  }
};

template <typename T, std::size_t N>
struct SinCosMcLoughlinExpansionFirstLoop<T, N, 0> {
  static void compute(T &c, T &s, const T &z) {
    c = c * z + static_cast<T>(COS_MACLAURIN_DOUBLEANGLE_FACTOR[N - 1]);
    s = s * z + static_cast<T>(SIN_MACLAURIN_DOUBLEANGLE_FACTOR[N - 1]);
  }
};

template <typename T, std::size_t N>
struct SinCosMcLoughlinExpansionSecondLoop {
  static void compute(T &s, T &c) {
    s = s * (static_cast<T>(2) - c);
    c = c * (static_cast<T>(4) - c);
    SinCosMcLoughlinExpansionSecondLoop<T, N - 1>::compute(s, c);
  }
};

template <typename T> struct SinCosMcLoughlinExpansionSecondLoop<T, 0> {
  static void compute(T &s, T &c) {
    static_cast<void>(s);
    static_cast<void>(c);
    /* Do Nothing. */
  }
};

template <typename T, std::size_t MACLAURIN_EXPANSION_REPEAT_NUMBER>
inline void sincos_maclaurin_expansion_with_DoubleAngleFormula(const T &theta,
                                                               T &cos_value,
                                                               T &sin_value) {
  static_assert(MACLAURIN_EXPANSION_REPEAT_NUMBER <
                    COS_MACLAURIN_DOUBLEANGLE_FACTOR_MAX_SIZE,
                "MACLAURIN_EXPANSION_REPEAT_NUMBER is too large.");

  T theta_wrapped = Base::Math::wrap_value_in_minus_pi_and_pi(theta);

  T c = static_cast<T>(
      COS_MACLAURIN_DOUBLEANGLE_FACTOR[MACLAURIN_EXPANSION_REPEAT_NUMBER]);
  T s = static_cast<T>(
      SIN_MACLAURIN_DOUBLEANGLE_FACTOR[MACLAURIN_EXPANSION_REPEAT_NUMBER]);
  T z = static_cast<T>(0);

  theta_wrapped =
      theta_wrapped /
      static_cast<T>(Factorial_2<MACLAURIN_EXPANSION_REPEAT_NUMBER>::value);
  z = theta_wrapped * theta_wrapped;

  SinCosMcLoughlinExpansionFirstLoop<T, MACLAURIN_EXPANSION_REPEAT_NUMBER,
                                     (MACLAURIN_EXPANSION_REPEAT_NUMBER -
                                      1)>::compute(c, s, z);

  c = c * z;
  s = s * theta_wrapped;

  SinCosMcLoughlinExpansionSecondLoop<
      T, MACLAURIN_EXPANSION_REPEAT_NUMBER>::compute(s, c);

  cos_value = static_cast<T>(1) - c * static_cast<T>(0.5);
  sin_value = s;
}

/* sin maclaurin expansion with DoubleAngleFormula */
template <typename T, std::size_t MACLAURIN_EXPANSION_REPEAT_NUMBER>
inline T sin_maclaurin_expansion_with_DoubleAngleFormula(const T &x) {

  return Base::Math::cos_maclaurin_expansion_with_DoubleAngleFormula<
      T, MACLAURIN_EXPANSION_REPEAT_NUMBER>(
      x - static_cast<T>(Base::Math::HALF_PI));
}

/* tan maclaurin expansion with DoubleAngleFormula */
template <typename T, std::size_t MACLAURIN_EXPANSION_REPEAT_NUMBER>
inline T tan_maclaurin_expansion_with_DoubleAngleFormula(const T &x) {
  T cos_value = static_cast<T>(0);
  T sin_value = static_cast<T>(0);

  Base::Math::sincos_maclaurin_expansion_with_DoubleAngleFormula<
      T, MACLAURIN_EXPANSION_REPEAT_NUMBER>(x, cos_value, sin_value);

  return sin_value /
         Base::Utility::avoid_zero_divide(
             cos_value, static_cast<T>(Base::Math::TRIGONOMETRIC_DIVISION_MIN));
}

/* sin maclaurin expansion */
template <typename T, std::size_t LOOP_MAX, std::size_t N>
struct SinMaclaurinLoop {
  static void compute(const T &x_squared, T &term, T &result) {
    term *=
        -x_squared * static_cast<T>(static_cast<T>(1) /
                                    static_cast<T>((2 * (LOOP_MAX - N)) *
                                                   (2 * (LOOP_MAX - N) + 1)));
    result += term;

    SinMaclaurinLoop<T, LOOP_MAX, N - 1>::compute(x_squared, term, result);
  }
};

template <typename T, std::size_t LOOP_MAX>
struct SinMaclaurinLoop<T, LOOP_MAX, 0> {
  static void compute(const T &x_squared, T &term, T &result) {
    /* Do Nothing. */
    static_cast<void>(x_squared);
    static_cast<void>(term);
    static_cast<void>(result);
  }
};

template <typename T, std::size_t LOOP_NUMBER>
inline T sin_maclaurin_expansion(const T &x) {

  T x_wrapped = Base::Math::wrap_value_in_minus_pi_and_pi(x);

  T term = x_wrapped;
  T result = x_wrapped;
  T x_squared = x_wrapped * x_wrapped;

  SinMaclaurinLoop<T, LOOP_NUMBER, LOOP_NUMBER - 1>::compute(x_squared, term,
                                                             result);

  return result;
}

/* sin */
template <typename T> inline T sin(const T &x) {

#ifdef __BASE_MATH_USE_STD_MATH__
  return std::sin(x);
#else // __BASE_MATH_USE_STD_MATH__

  return Base::Math::sin_maclaurin_expansion_with_DoubleAngleFormula<
      T, Base::Math::SIN_MACLAURIN_DOUBLEANGLE_REPEAT_NUMBER>(x);

#endif // __BASE_MATH_USE_STD_MATH__
}

/* cos maclaurin expansion */
template <typename T, std::size_t LOOP_MAX, std::size_t N>
struct CosMaclaurinLoop {
  static void compute(const T &x_squared, T &term, T &result) {
    term *=
        -x_squared * static_cast<T>(static_cast<T>(1) /
                                    static_cast<T>((2 * (LOOP_MAX - N) - 1) *
                                                   (2 * (LOOP_MAX - N))));
    result += term;

    CosMaclaurinLoop<T, LOOP_MAX, N - 1>::compute(x_squared, term, result);
  }
};

template <typename T, std::size_t LOOP_MAX>
struct CosMaclaurinLoop<T, LOOP_MAX, 0> {
  static void compute(const T &x_squared, T &term, T &result) {
    /* Do Nothing. */
    static_cast<void>(x_squared);
    static_cast<void>(term);
    static_cast<void>(result);
  }
};

template <typename T, std::size_t LOOP_NUMBER>
inline T cos_maclaurin_expansion(const T &x) {

  T x_wrapped = wrap_value_in_minus_pi_and_pi(x);

  T term = static_cast<T>(1);
  T result = static_cast<T>(1);
  T x_squared = x_wrapped * x_wrapped;

  CosMaclaurinLoop<T, LOOP_NUMBER, LOOP_NUMBER - 1>::compute(x_squared, term,
                                                             result);

  return result;
}

/* cos */
template <typename T> inline T cos(const T &x) {

#ifdef __BASE_MATH_USE_STD_MATH__
  return std::cos(x);
#else // __BASE_MATH_USE_STD_MATH__

  return Base::Math::cos_maclaurin_expansion_with_DoubleAngleFormula<
      T, Base::Math::COS_MACLAURIN_DOUBLEANGLE_REPEAT_NUMBER>(x);

#endif // __BASE_MATH_USE_STD_MATH__
}

/* tan maclaurin expansion */
template <typename T, std::size_t SIN_LOOP_NUMBER, std::size_t COS_LOOP_NUMBER>
inline T tan_maclaurin_expansion(const T &x) {

  return Base::Math::sin_maclaurin_expansion<T, SIN_LOOP_NUMBER>(x) /
         Base::Utility::avoid_zero_divide(
             Base::Math::cos_maclaurin_expansion<T, COS_LOOP_NUMBER>(x),
             static_cast<T>(Base::Math::TRIGONOMETRIC_DIVISION_MIN));
}

/* tan */
template <typename T> inline T tan(const T &x) {

#ifdef __BASE_MATH_USE_STD_MATH__
  return std::tan(x);
#else // __BASE_MATH_USE_STD_MATH__

  return Base::Math::tan_maclaurin_expansion_with_DoubleAngleFormula<
      T, Base::Math::SINCOS_MACLAURIN_DOUBLEANGLE_REPEAT_NUMBER>(x);

#endif // __BASE_MATH_USE_STD_MATH__
}

/* atan maclaurin expansion */
template <typename T, std::size_t LOOP_NUMBER>
inline T atan_maclaurin_expansion(const T &x) {

  T result = static_cast<T>(0);

  if (x > static_cast<T>(1)) {
    result = static_cast<T>(Base::Math::HALF_PI) -
             Base::Math::atan_maclaurin_expansion(static_cast<T>(1) / x);
  } else if (x < static_cast<T>(-1)) {
    result = -static_cast<T>(Base::Math::HALF_PI) -
             Base::Math::atan_maclaurin_expansion(static_cast<T>(1) / x);
  } else {
    if ((x > static_cast<T>(0.5)) || (x < static_cast<T>(-0.5))) {
      T half_x = x * static_cast<T>(0.5);

      result = Base::Math::atan_maclaurin_expansion(half_x) +
               Base::Math::atan_maclaurin_expansion(
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

/* atan Chebyshev */
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

  static_assert(LOOP_NUMBER <= Base::Math::CHEBYSHEV_COEFFICIENT_FOR_ATAN_SIZE,
                "LOOP_NUMBER is too large.");

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

/* atan */
template <typename T> inline T atan(const T &x) {

#ifdef __BASE_MATH_USE_STD_MATH__
  return std::atan(x);
#else // __BASE_MATH_USE_STD_MATH__

  return Base::Math::atan_chebyshev<T, Base::Math::ATAN_REPEAT_NUMBER>(x);

#endif // __BASE_MATH_USE_STD_MATH__
}

/* atan2 Chebyshev */
template <typename T, std::size_t LOOP_NUMBER>
inline T atan2_chebyshev(const T &y, const T &x) {
  T result = static_cast<T>(0);

  if (Base::Utility::near_zero(
          x, static_cast<T>(Base::Math::TRIGONOMETRIC_DIVISION_MIN)) &&
      Base::Utility::near_zero(
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

  } else if (Base::Utility::near_zero(
                 x, static_cast<T>(Base::Math::TRIGONOMETRIC_DIVISION_MIN)) &&
             y > static_cast<T>(0)) {

    result = static_cast<T>(Base::Math::HALF_PI);

  } else if (Base::Utility::near_zero(
                 x, static_cast<T>(Base::Math::TRIGONOMETRIC_DIVISION_MIN)) &&
             y < static_cast<T>(0)) {

    result = -static_cast<T>(Base::Math::HALF_PI);
  }

  return result;
}

/* atan2 */
template <typename T> inline T atan2(const T &y, const T &x) {

#ifdef __BASE_MATH_USE_STD_MATH__
  return std::atan2(y, x);
#else // __BASE_MATH_USE_STD_MATH__

  return Base::Math::atan2_chebyshev<T, Base::Math::ATAN_REPEAT_NUMBER>(y, x);

#endif // __BASE_MATH_USE_STD_MATH__
}

/* asin Chebyshev */
template <typename T, std::size_t ATAN_LOOP_NUMBER,
          std::size_t SQRT_LOOP_NUMBER>
inline T asin_chebyshev(const T &x) {

  T result = static_cast<T>(0);

  if (x >= static_cast<T>(1)) {

    result = static_cast<T>(Base::Math::HALF_PI);

  } else if (x <= static_cast<T>(-1)) {

    result = -static_cast<T>(Base::Math::HALF_PI);

  } else {

    result = static_cast<T>(2) *
             Base::Math::atan_chebyshev<T, ATAN_LOOP_NUMBER>(
                 x / (static_cast<T>(1) +
                      Base::Math::sqrt_newton_method<T, SQRT_LOOP_NUMBER>(
                          static_cast<T>(1) - x * x)));
  }

  return result;
}

/* asin */
template <typename T> inline T asin(const T &x) {

#ifdef __BASE_MATH_USE_STD_MATH__

  if (x >= static_cast<T>(1)) {

    return static_cast<T>(Base::Math::HALF_PI);

  } else if (x <= static_cast<T>(-1)) {

    return -static_cast<T>(Base::Math::HALF_PI);

  } else {

    return std::asin(x);
  }
#else // __BASE_MATH_USE_STD_MATH__

  return Base::Math::asin_chebyshev<T, Base::Math::ATAN_REPEAT_NUMBER,
                                    Base::Math::SQRT_REPEAT_NUMBER>(x);

#endif // __BASE_MATH_USE_STD_MATH__
}

/* acos Chebyshev */
template <typename T, std::size_t ATAN_LOOP_NUMBER,
          std::size_t SQRT_LOOP_NUMBER>
inline T acos_chebyshev(const T &x) {

  return static_cast<T>(Base::Math::HALF_PI) -
         Base::Math::asin_chebyshev<T, ATAN_LOOP_NUMBER, SQRT_LOOP_NUMBER>(x);
}

/* acos */
template <typename T> inline T acos(const T &x) {

#ifdef __BASE_MATH_USE_STD_MATH__

  if (x > static_cast<T>(1)) {

    return static_cast<T>(0);

  } else if (x < static_cast<T>(-1)) {

    return static_cast<T>(Base::Math::PI);

  } else {

    return std::acos(x);
  }
#else // __BASE_MATH_USE_STD_MATH__

  return Base::Math::acos_chebyshev<T, Base::Math::ATAN_REPEAT_NUMBER,
                                    Base::Math::SQRT_REPEAT_NUMBER>(x);

#endif // __BASE_MATH_USE_STD_MATH__
}

/* sinh Maclaurin Expansion with table */
template <typename T, std::size_t LOOP_NUMBER>
inline T sinh_maclaurin_expansion_with_table(const T &x) {
  return (Base::Math::exp_maclaurin_expansion_with_table<T, LOOP_NUMBER>(x) -
          Base::Math::exp_maclaurin_expansion_with_table<T, LOOP_NUMBER>(-x)) *
         static_cast<T>(0.5);
}

/* sinh maclaurin expansion */
template <typename T, std::size_t LOOP_NUMBER>
inline T sinh_maclaurin_expansion(const T &x) {
  return (Base::Math::exp_maclaurin_expansion<T, LOOP_NUMBER>(x) -
          Base::Math::exp_maclaurin_expansion<T, LOOP_NUMBER>(-x)) *
         static_cast<T>(0.5);
}

/* sinh */
template <typename T> inline T sinh(const T &x) {

#ifdef __BASE_MATH_USE_STD_MATH__
  return std::sinh(x);
#else // __BASE_MATH_USE_STD_MATH__

  return Base::Math::sinh_maclaurin_expansion<T, Base::Math::EXP_REPEAT_NUMBER>(
      x);

#endif // __BASE_MATH_USE_STD_MATH__
}

/* cosh Maclaurin Expansion with table */
template <typename T, std::size_t LOOP_NUMBER>
inline T cosh_maclaurin_expansion_with_table(const T &x) {
  return (Base::Math::exp_maclaurin_expansion_with_table<T, LOOP_NUMBER>(x) +
          Base::Math::exp_maclaurin_expansion_with_table<T, LOOP_NUMBER>(-x)) *
         static_cast<T>(0.5);
}

/* cosh maclaurin expansion */
template <typename T, std::size_t LOOP_NUMBER>
inline T cosh_maclaurin_expansion(const T &x) {
  return (Base::Math::exp_maclaurin_expansion<T, LOOP_NUMBER>(x) +
          Base::Math::exp_maclaurin_expansion<T, LOOP_NUMBER>(-x)) *
         static_cast<T>(0.5);
}

/* cosh */
template <typename T> inline T cosh(const T &x) {

#ifdef __BASE_MATH_USE_STD_MATH__
  return std::cosh(x);
#else // __BASE_MATH_USE_STD_MATH__

  return Base::Math::cosh_maclaurin_expansion<T, Base::Math::EXP_REPEAT_NUMBER>(
      x);

#endif // __BASE_MATH_USE_STD_MATH__
}

/* tanh maclaurin expansion with table */
template <typename T, std::size_t LOOP_NUMBER>
inline T tanh_maclaurin_expansion_with_table(const T &x) {

  T a = Base::Math::exp_maclaurin_expansion_with_table<T, LOOP_NUMBER>(x);
  T b = Base::Math::exp_maclaurin_expansion_with_table<T, LOOP_NUMBER>(-x);

  return (a - b) / Base::Utility::avoid_zero_divide(
                       a + b, static_cast<T>(TRIGONOMETRIC_DIVISION_MIN));
}

/* tanh maclaurin expansion */
template <typename T, std::size_t LOOP_NUMBER>
inline T tanh_maclaurin_expansion(const T &x) {

  T a = Base::Math::exp_maclaurin_expansion<T, LOOP_NUMBER>(x);
  T b = Base::Math::exp_maclaurin_expansion<T, LOOP_NUMBER>(-x);

  return (a - b) / Base::Utility::avoid_zero_divide(
                       a + b, static_cast<T>(TRIGONOMETRIC_DIVISION_MIN));
}

/* tanh */
template <typename T> inline T tanh(const T &x) {

#ifdef __BASE_MATH_USE_STD_MATH__
  return std::tanh(x);
#else // __BASE_MATH_USE_STD_MATH__

  return Base::Math::tanh_maclaurin_expansion<T, Base::Math::EXP_REPEAT_NUMBER>(
      x);

#endif // __BASE_MATH_USE_STD_MATH__
}

} // namespace Math
} // namespace Base

#endif // __BASE_MATH_TRIGONOMETRIC_HPP__
