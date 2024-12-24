#ifndef BASE_MATH_TRIGONOMETRIC_HPP
#define BASE_MATH_TRIGONOMETRIC_HPP

#include "base_math_macros.hpp"

#include "base_math_arithmetic.hpp"
#include "base_math_exponential_logarithmic.hpp"
#include "base_math_mathematical_constants.hpp"
#include "base_utility.hpp"

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

constexpr std::size_t COS_MCLOUGHLIN_DOUBLEANGLE_REPEAT_NUMBER = 1;
constexpr std::size_t SIN_MCLOUGHLIN_DOUBLEANGLE_REPEAT_NUMBER =
    COS_MCLOUGHLIN_DOUBLEANGLE_REPEAT_NUMBER;

constexpr std::size_t SINCOS_MCLOUGHLIN_DOUBLEANGLE_REPEAT_NUMBER = 2;

#else // BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS
constexpr std::size_t SIN_REPEAT_NUMBER = 8;
constexpr std::size_t COS_REPEAT_NUMBER = 9;
constexpr std::size_t ATAN_REPEAT_NUMBER = 8;

constexpr std::size_t COS_MCLOUGHLIN_DOUBLEANGLE_REPEAT_NUMBER = 3;
constexpr std::size_t SIN_MCLOUGHLIN_DOUBLEANGLE_REPEAT_NUMBER =
    COS_MCLOUGHLIN_DOUBLEANGLE_REPEAT_NUMBER;

constexpr std::size_t SINCOS_MCLOUGHLIN_DOUBLEANGLE_REPEAT_NUMBER = 3;

#endif // BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS

constexpr std::size_t COS_MCLOUGHLIN_DOUBLEANGLE_FACTOR_MAX_SIZE = 6;

using COS_MCLOUGHLIN_FACTOR_LIST = typename MakeCosMcloughlinFactorList<
    COS_MCLOUGHLIN_DOUBLEANGLE_FACTOR_MAX_SIZE>::type;

constexpr auto COS_MCLOUGHLIN_DOUBLEANGLE_FACTOR =
    Base::Math::to_cos_mcloughlin_factor_array(COS_MCLOUGHLIN_FACTOR_LIST{});

constexpr std::size_t SIN_MCLOUGHLIN_DOUBLEANGLE_FACTOR_MAX_SIZE = 6;

using SIN_MCLOUGHLIN_FACTOR_LIST = typename MakeSinMcloughlinFactorList<
    SIN_MCLOUGHLIN_DOUBLEANGLE_FACTOR_MAX_SIZE>::type;

constexpr auto SIN_MCLOUGHLIN_DOUBLEANGLE_FACTOR =
    Base::Math::to_sin_mcloughlin_factor_array(SIN_MCLOUGHLIN_FACTOR_LIST{});

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

/* cos mcloughlin expansion with DoubleAngleFormula */
template <typename T, std::size_t N,
          std::size_t MCLOUGHLIN_EXPANSION_REPEAT_NUMBER>
struct CosMcloughlinExpansionFirstLoop {
  static void compute(T &y, T &x_wrapped) {
    CosMcloughlinExpansionFirstLoop<
        T, N - 1, MCLOUGHLIN_EXPANSION_REPEAT_NUMBER>::compute(y, x_wrapped);
    y = y * x_wrapped +
        static_cast<T>(COS_MCLOUGHLIN_DOUBLEANGLE_FACTOR
                           [MCLOUGHLIN_EXPANSION_REPEAT_NUMBER - 1 - N]);
  }
};

template <typename T, std::size_t MCLOUGHLIN_EXPANSION_REPEAT_NUMBER>
struct CosMcloughlinExpansionFirstLoop<T, 0,
                                       MCLOUGHLIN_EXPANSION_REPEAT_NUMBER> {
  static void compute(T &y, T &x_wrapped) {
    y = y * x_wrapped +
        static_cast<T>(COS_MCLOUGHLIN_DOUBLEANGLE_FACTOR
                           [MCLOUGHLIN_EXPANSION_REPEAT_NUMBER - 1]);
  }
};

template <typename T, std::size_t N> struct CosMcloughlinExpansionSecondLoop {
  static void compute(T &y) {
    y = y * (static_cast<T>(4) - y);
    CosMcloughlinExpansionSecondLoop<T, N - 1>::compute(y);
  }
};

template <typename T> struct CosMcloughlinExpansionSecondLoop<T, 0> {
  static void compute(T &y) {
    static_cast<void>(y);
    /* Do Nothing. */
  }
};

template <typename T, std::size_t MCLOUGHLIN_EXPANSION_REPEAT_NUMBER>
inline T cos_mcloughlin_expansion_with_DoubleAngleFormula(const T &x) {
  static_assert(MCLOUGHLIN_EXPANSION_REPEAT_NUMBER <
                    COS_MCLOUGHLIN_DOUBLEANGLE_FACTOR_MAX_SIZE,
                "MCLOUGHLIN_EXPANSION_REPEAT_NUMBER is too large.");

  T x_wrapped = Base::Math::wrap_value_in_minus_pi_and_pi(x);
  T y = static_cast<T>(0);

  x_wrapped =
      x_wrapped /
      static_cast<T>(Factorial_2<MCLOUGHLIN_EXPANSION_REPEAT_NUMBER>::value);
  x_wrapped = x_wrapped * x_wrapped;

  y = static_cast<T>(
      COS_MCLOUGHLIN_DOUBLEANGLE_FACTOR[MCLOUGHLIN_EXPANSION_REPEAT_NUMBER]);

  CosMcloughlinExpansionFirstLoop<
      T, (MCLOUGHLIN_EXPANSION_REPEAT_NUMBER - 1),
      MCLOUGHLIN_EXPANSION_REPEAT_NUMBER>::compute(y, x_wrapped);

  y = y * x_wrapped;

  CosMcloughlinExpansionSecondLoop<
      T, MCLOUGHLIN_EXPANSION_REPEAT_NUMBER>::compute(y);

  return static_cast<T>(1) - y * static_cast<T>(0.5);
}

/* sin cos Mcloughlin expansion with DoubleAngleFormula */
template <typename T, std::size_t N, std::size_t I>
struct SinCosMcLoughlinExpansionFirstLoop {
  static void compute(T &c, T &s, const T &z) {
    SinCosMcLoughlinExpansionFirstLoop<T, N, I - 1>::compute(c, s, z);
    c = c * z + static_cast<T>(COS_MCLOUGHLIN_DOUBLEANGLE_FACTOR[N - 1 - I]);
    s = s * z + static_cast<T>(SIN_MCLOUGHLIN_DOUBLEANGLE_FACTOR[N - 1 - I]);
  }
};

template <typename T, std::size_t N>
struct SinCosMcLoughlinExpansionFirstLoop<T, N, 0> {
  static void compute(T &c, T &s, const T &z) {
    c = c * z + static_cast<T>(COS_MCLOUGHLIN_DOUBLEANGLE_FACTOR[N - 1]);
    s = s * z + static_cast<T>(SIN_MCLOUGHLIN_DOUBLEANGLE_FACTOR[N - 1]);
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

template <typename T, std::size_t MCLOUGHLIN_EXPANSION_REPEAT_NUMBER>
inline void sincos_mcloughlin_expansion_with_DoubleAngleFormula(const T &theta,
                                                                T &cos_value,
                                                                T &sin_value) {
  static_assert(MCLOUGHLIN_EXPANSION_REPEAT_NUMBER <
                    COS_MCLOUGHLIN_DOUBLEANGLE_FACTOR_MAX_SIZE,
                "MCLOUGHLIN_EXPANSION_REPEAT_NUMBER is too large.");

  T theta_wrapped = Base::Math::wrap_value_in_minus_pi_and_pi(theta);

  T c = static_cast<T>(
      COS_MCLOUGHLIN_DOUBLEANGLE_FACTOR[MCLOUGHLIN_EXPANSION_REPEAT_NUMBER]);
  T s = static_cast<T>(
      SIN_MCLOUGHLIN_DOUBLEANGLE_FACTOR[MCLOUGHLIN_EXPANSION_REPEAT_NUMBER]);
  T z = static_cast<T>(0);

  theta_wrapped =
      theta_wrapped /
      static_cast<T>(Factorial_2<MCLOUGHLIN_EXPANSION_REPEAT_NUMBER>::value);
  z = theta_wrapped * theta_wrapped;

  SinCosMcLoughlinExpansionFirstLoop<T, MCLOUGHLIN_EXPANSION_REPEAT_NUMBER,
                                     (MCLOUGHLIN_EXPANSION_REPEAT_NUMBER -
                                      1)>::compute(c, s, z);

  c = c * z;
  s = s * theta_wrapped;

  SinCosMcLoughlinExpansionSecondLoop<
      T, MCLOUGHLIN_EXPANSION_REPEAT_NUMBER>::compute(s, c);

  cos_value = static_cast<T>(1) - c * static_cast<T>(0.5);
  sin_value = s;
}

/* sin mcloughlin expansion with DoubleAngleFormula */
template <typename T, std::size_t MCLOUGHLIN_EXPANSION_REPEAT_NUMBER>
inline T sin_mcloughlin_expansion_with_DoubleAngleFormula(const T &x) {

  return Base::Math::cos_mcloughlin_expansion_with_DoubleAngleFormula<
      T, MCLOUGHLIN_EXPANSION_REPEAT_NUMBER>(
      x - static_cast<T>(Base::Math::HALF_PI));
}

/* tan mcloughlin expansion with DoubleAngleFormula */
template <typename T, std::size_t MCLOUGHLIN_EXPANSION_REPEAT_NUMBER>
inline T tan_mcloughlin_expansion_with_DoubleAngleFormula(const T &x) {
  T cos_value = static_cast<T>(0);
  T sin_value = static_cast<T>(0);

  Base::Math::sincos_mcloughlin_expansion_with_DoubleAngleFormula<
      T, MCLOUGHLIN_EXPANSION_REPEAT_NUMBER>(x, cos_value, sin_value);

  return sin_value /
         Base::Utility::avoid_zero_divide(
             cos_value, static_cast<T>(Base::Math::TRIGONOMETRIC_DIVISION_MIN));
}

/* sin mcloughlin expansion */
template <typename T, std::size_t LOOP_MAX, std::size_t N>
struct SinMcloughlinLoop {
  static void compute(const T &x_squared, T &term, T &result) {
    term *=
        -x_squared * static_cast<T>(static_cast<T>(1) /
                                    static_cast<T>((2 * (LOOP_MAX - N)) *
                                                   (2 * (LOOP_MAX - N) + 1)));
    result += term;

    SinMcloughlinLoop<T, LOOP_MAX, N - 1>::compute(x_squared, term, result);
  }
};

template <typename T, std::size_t LOOP_MAX>
struct SinMcloughlinLoop<T, LOOP_MAX, 0> {
  static void compute(const T &x_squared, T &term, T &result) {
    /* Do Nothing. */
    static_cast<void>(x_squared);
    static_cast<void>(term);
    static_cast<void>(result);
  }
};

template <typename T, std::size_t LOOP_NUMBER>
inline T sin_mcloughlin_expansion(const T &x) {

  T x_wrapped = Base::Math::wrap_value_in_minus_pi_and_pi(x);

  T term = x_wrapped;
  T result = x_wrapped;
  T x_squared = x_wrapped * x_wrapped;

  SinMcloughlinLoop<T, LOOP_NUMBER, LOOP_NUMBER - 1>::compute(x_squared, term,
                                                              result);

  return result;
}

template <typename T> inline T sin(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::sin(x);
#else
  return Base::Math::sin_mcloughlin_expansion<T, Base::Math::SIN_REPEAT_NUMBER>(
      x);
#endif
}

/* cos mcloughlin expansion */
template <typename T, std::size_t LOOP_MAX, std::size_t N>
struct CosMcloughlinLoop {
  static void compute(const T &x_squared, T &term, T &result) {
    term *=
        -x_squared * static_cast<T>(static_cast<T>(1) /
                                    static_cast<T>((2 * (LOOP_MAX - N) - 1) *
                                                   (2 * (LOOP_MAX - N))));
    result += term;

    CosMcloughlinLoop<T, LOOP_MAX, N - 1>::compute(x_squared, term, result);
  }
};

template <typename T, std::size_t LOOP_MAX>
struct CosMcloughlinLoop<T, LOOP_MAX, 0> {
  static void compute(const T &x_squared, T &term, T &result) {
    /* Do Nothing. */
    static_cast<void>(x_squared);
    static_cast<void>(term);
    static_cast<void>(result);
  }
};

template <typename T, std::size_t LOOP_NUMBER>
inline T cos_mcloughlin_expansion(const T &x) {

  T x_wrapped = wrap_value_in_minus_pi_and_pi(x);

  T term = static_cast<T>(1);
  T result = static_cast<T>(1);
  T x_squared = x_wrapped * x_wrapped;

  CosMcloughlinLoop<T, LOOP_NUMBER, LOOP_NUMBER - 1>::compute(x_squared, term,
                                                              result);

  return result;
}

template <typename T> inline T cos(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::cos(x);
#else
  return Base::Math::cos_mcloughlin_expansion<T, Base::Math::COS_REPEAT_NUMBER>(
      x);
#endif
}

/* tan mcloughlin expansion */
template <typename T, std::size_t SIN_LOOP_NUMBER, std::size_t COS_LOOP_NUMBER>
inline T tan_mcloughlin_expansion(const T &x) {

  return Base::Math::sin_mcloughlin_expansion<T, SIN_LOOP_NUMBER>(x) /
         Base::Utility::avoid_zero_divide(
             Base::Math::cos_mcloughlin_expansion<T, COS_LOOP_NUMBER>(x),
             static_cast<T>(Base::Math::TRIGONOMETRIC_DIVISION_MIN));
}

template <typename T> inline T tan(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::tan(x);
#else
  return Base::Math::tan_mcloughlin_expansion<T, Base::Math::SIN_REPEAT_NUMBER,
                                              Base::Math::COS_REPEAT_NUMBER>(x);
#endif
}

/* atan mcloughlin expansion */
template <typename T, std::size_t LOOP_NUMBER>
inline T atan_mcloughlin_expansion(const T &x) {

  T result = static_cast<T>(0);

  if (x > static_cast<T>(1)) {
    result = static_cast<T>(Base::Math::HALF_PI) -
             Base::Math::atan_mcloughlin_expansion(static_cast<T>(1) / x);
  } else if (x < static_cast<T>(-1)) {
    result = -static_cast<T>(Base::Math::HALF_PI) -
             Base::Math::atan_mcloughlin_expansion(static_cast<T>(1) / x);
  } else {
    if ((x > static_cast<T>(0.5)) || (x < static_cast<T>(-0.5))) {
      T half_x = x * static_cast<T>(0.5);

      result = Base::Math::atan_mcloughlin_expansion(half_x) +
               Base::Math::atan_mcloughlin_expansion(
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

template <typename T> inline T atan2(const T &y, const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::atan2(y, x);
#else
  return Base::Math::atan2_chebyshev<T, Base::Math::ATAN_REPEAT_NUMBER>(y, x);
#endif
}

/* asin */
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
             Base::Math::atan2_chebyshev<T, ATAN_LOOP_NUMBER>(
                 x, static_cast<T>(1) +
                        Base::Math::sqrt_newton_method<T, SQRT_LOOP_NUMBER>(
                            static_cast<T>(1) - x * x));
  }

  return result;
}

template <typename T> inline T asin(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::asin(x);
#else
  return Base::Math::asin_chebyshev<T, Base::Math::ATAN_REPEAT_NUMBER,
                                    Base::Math::SQRT_REPEAT_NUMBER>(x);
#endif
}

/* acos */
template <typename T, std::size_t ATAN_LOOP_NUMBER,
          std::size_t SQRT_LOOP_NUMBER>
inline T acos_chebyshev(const T &x) {

  return static_cast<T>(Base::Math::HALF_PI) -
         Base::Math::asin_chebyshev<T, ATAN_LOOP_NUMBER, SQRT_LOOP_NUMBER>(x);
}

template <typename T> inline T acos(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::acos(x);
#else
  return Base::Math::acos_chebyshev<T, Base::Math::ATAN_REPEAT_NUMBER,
                                    Base::Math::SQRT_REPEAT_NUMBER>(x);
#endif
}

/* sinh */
template <typename T, std::size_t LOOP_NUMBER>
inline T sinh_mcloughlin_expansion(const T &x) {
  return (Base::Math::exp_mcloughlin_expansion<T, LOOP_NUMBER>(x) -
          Base::Math::exp_mcloughlin_expansion<T, LOOP_NUMBER>(-x)) *
         static_cast<T>(0.5);
}

template <typename T> inline T sinh(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::sinh(x);
#else
  return Base::Math::sinh_mcloughlin_expansion<T,
                                               Base::Math::EXP_REPEAT_NUMBER>(
      x);
#endif
}

/* cosh */
template <typename T, std::size_t LOOP_NUMBER>
inline T cosh_mcloughlin_expansion(const T &x) {
  return (Base::Math::exp_mcloughlin_expansion<T, LOOP_NUMBER>(x) +
          Base::Math::exp_mcloughlin_expansion<T, LOOP_NUMBER>(-x)) *
         static_cast<T>(0.5);
}

template <typename T> inline T cosh(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::cosh(x);
#else
  return Base::Math::cosh_mcloughlin_expansion<T,
                                               Base::Math::EXP_REPEAT_NUMBER>(
      x);
#endif
}

/* tanh */
template <typename T, std::size_t LOOP_NUMBER>
inline T tanh_mcloughlin_expansion(const T &x) {

  T a = Base::Math::exp_mcloughlin_expansion<T, LOOP_NUMBER>(x);
  T b = Base::Math::exp_mcloughlin_expansion<T, LOOP_NUMBER>(-x);

  return (a - b) / Base::Utility::avoid_zero_divide(
                       a + b, static_cast<T>(TRIGONOMETRIC_DIVISION_MIN));
}

template <typename T> inline T tanh(const T &x) {

#ifdef BASE_MATH_USE_STD_MATH
  return std::tanh(x);
#else
  return Base::Math::tanh_mcloughlin_expansion<T,
                                               Base::Math::EXP_REPEAT_NUMBER>(
      x);
#endif
}

} // namespace Math
} // namespace Base

#endif // BASE_MATH_TRIGONOMETRIC_HPP
