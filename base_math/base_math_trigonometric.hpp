/**
 * @file base_math_trigonometric.hpp
 * @brief Trigonometric and Hyperbolic Math Functions and Approximations
 *
 * This header provides a comprehensive set of mathematical functions for
 * trigonometric (sin, cos, tan, atan, asin, acos, atan2) and hyperbolic (sinh,
 * cosh, tanh) operations, along with their various polynomial and Chebyshev
 * approximations. The implementations are designed to be highly configurable,
 * supporting both standard library and custom approximation methods for
 * performance and precision tuning.
 *
 * The code is organized under the `Base::Math` namespace and includes:
 *   - Maclaurin series expansions for sine, cosine, tangent, and their
 * hyperbolic counterparts.
 *   - Double angle formula optimizations for improved convergence.
 *   - Chebyshev polynomial approximations for arctangent and related functions.
 *   - Utility functions for angle wrapping and safe division.
 *   - Compile-time configuration for approximation accuracy and speed.
 *
 * Main Classes and Namespaces:
 *   - `Base::Math`: Main namespace containing all trigonometric and hyperbolic
 * functions.
 *   - `CosMaclaurinExpansion`: Internal namespace for cosine Maclaurin
 * expansion with double angle formula.
 *   - `SinCosMcLoughlinExpansion`: Internal namespace for simultaneous sine and
 * cosine Maclaurin expansion.
 *   - `SinMaclaurin`, `CosMaclaurin`: Internal namespaces for Maclaurin series
 * loop implementations.
 *   - `AtanChebyshev`: Internal namespace for Chebyshev polynomial evaluation
 * of arctangent.
 *
 * The code is highly templated to support various numeric types and
 * compile-time loop unrolling. It is intended for use in performance-critical
 * or numerically sensitive applications where standard library math functions
 * may not be suitable.
 */
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
/**
 * @brief Wraps an angle value to the range [-π, π).
 *
 * This function takes an input angle `x` (in radians) and returns an equivalent
 * angle within the interval [-π, π). It uses modular arithmetic to ensure the
 * result is always within this range, which is useful for trigonometric
 * calculations and angle normalization.
 *
 * @tparam T Numeric type of the input and output (e.g., float, double).
 * @param x The angle value to be wrapped (in radians).
 * @return The wrapped angle value in the range [-π, π).
 */
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

namespace CosMaclaurinExpansion {

/**
 * @brief Recursive template struct to compute a value using a loop unrolling
 * technique.
 *
 * This struct performs a recursive computation, typically used for evaluating
 * Maclaurin series expansions (such as for trigonometric functions) at compile
 * time. The recursion is controlled by the template parameter N, which
 * decrements with each instantiation until the base case is reached (not shown
 * in this snippet).
 *
 * @tparam T The numeric type used for computation (e.g., float, double).
 * @tparam N The current recursion depth or loop index.
 * @tparam MACLAURIN_EXPANSION_REPEAT_NUMBER The total number of terms in the
 * Maclaurin expansion.
 *
 * @param y Reference to the result variable, updated at each recursion step.
 * @param x_wrapped Reference to the input value, typically the variable for the
 * Maclaurin expansion.
 */
template <typename T, std::size_t N,
          std::size_t MACLAURIN_EXPANSION_REPEAT_NUMBER>
struct FirstLoop {
  static void compute(T &y, T &x_wrapped) {
    FirstLoop<T, N - 1, MACLAURIN_EXPANSION_REPEAT_NUMBER>::compute(y,
                                                                    x_wrapped);
    y = y * x_wrapped +
        static_cast<T>(
            COS_MACLAURIN_DOUBLEANGLE_FACTOR[MACLAURIN_EXPANSION_REPEAT_NUMBER -
                                             1 - N]);
  }
};

template <typename T, std::size_t MACLAURIN_EXPANSION_REPEAT_NUMBER>
struct FirstLoop<T, 0, MACLAURIN_EXPANSION_REPEAT_NUMBER> {
  /**
   * @brief Computes a value using the provided references y and x_wrapped.
   *
   * This static function updates the value of y by multiplying it with
   * x_wrapped and adding a constant factor from the
   * COS_MACLAURIN_DOUBLEANGLE_FACTOR array. The constant is selected based on
   * the MACLAURIN_EXPANSION_REPEAT_NUMBER.
   *
   * @tparam T The numeric type of the input and output variables.
   * @param y Reference to the value to be updated.
   * @param x_wrapped Reference to the value to be used in the computation.
   */
  static void compute(T &y, T &x_wrapped) {
    y = y * x_wrapped +
        static_cast<T>(
            COS_MACLAURIN_DOUBLEANGLE_FACTOR[MACLAURIN_EXPANSION_REPEAT_NUMBER -
                                             1]);
  }
};

template <typename T, std::size_t N> struct SecondLoop {
  /**
   * @brief Performs a recursive computation on the input value y.
   *
   * This static function multiplies y by (4 - y) and then recursively calls
   * itself via SecondLoop<T, N - 1>::compute(y) to continue the computation for
   * N - 1 iterations.
   *
   * @tparam T The numeric type of the input value y.
   * @tparam N The number of recursive iterations to perform.
   * @param y Reference to the value to be computed and updated in place.
   */
  static void compute(T &y) {
    y = y * (static_cast<T>(4) - y);
    SecondLoop<T, N - 1>::compute(y);
  }
};

template <typename T> struct SecondLoop<T, 0> {
  /**
   * @brief Computes a value and assigns it to the given parameter.
   *
   * This is a default implementation that performs no operation on the input
   * parameter. It is intended to be specialized or overridden for specific
   * computations.
   *
   * @tparam T The type of the parameter.
   * @param y Reference to the value to be computed or modified.
   */
  static void compute(T &y) {
    static_cast<void>(y);
    /* Do Nothing. */
  }
};

} // namespace CosMaclaurinExpansion

/**
 * @brief Computes the cosine of a value using the Maclaurin series expansion
 *        with the Double Angle Formula for improved accuracy and efficiency.
 *
 * This function wraps the input value into the range [-π, π], then evaluates
 * the cosine using a Maclaurin series expansion up to the specified number of
 * terms, leveraging precomputed double angle factors and optimized loops for
 * performance.
 *
 * @tparam T The floating-point type (e.g., float, double) for computation.
 * @tparam MACLAURIN_EXPANSION_REPEAT_NUMBER The number of terms to use in the
 * Maclaurin expansion. Must be less than
 * COS_MACLAURIN_DOUBLEANGLE_FACTOR_MAX_SIZE.
 * @param x The input value (in radians) for which to compute the cosine.
 * @return The approximate value of cos(x) computed using the Maclaurin
 * expansion.
 *
 * @note Requires supporting definitions:
 *       - Base::Math::wrap_value_in_minus_pi_and_pi
 *       - Factorial_2
 *       - COS_MACLAURIN_DOUBLEANGLE_FACTOR
 *       - CosMaclaurinExpansion::FirstLoop and SecondLoop
 *
 * @throws static_assert if MACLAURIN_EXPANSION_REPEAT_NUMBER is too large.
 */
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

  CosMaclaurinExpansion::FirstLoop<
      T, (MACLAURIN_EXPANSION_REPEAT_NUMBER - 1),
      MACLAURIN_EXPANSION_REPEAT_NUMBER>::compute(y, x_wrapped);

  y = y * x_wrapped;

  CosMaclaurinExpansion::SecondLoop<
      T, MACLAURIN_EXPANSION_REPEAT_NUMBER>::compute(y);

  return static_cast<T>(1) - y * static_cast<T>(0.5);
}

/* sin cos Maclaurin expansion with DoubleAngleFormula */
namespace SinCosMcLoughlinExpansion {

template <typename T, std::size_t N, std::size_t I> struct FirstLoop {
  /**
   * @brief Recursive template struct to compute coefficients for Maclaurin
   * series expansion.
   *
   * This struct performs a compile-time loop to recursively compute the
   * coefficients for the cosine and sine Maclaurin series using double-angle
   * factors.
   *
   * @tparam T  The numeric type used for computation (e.g., float, double).
   * @tparam N  The total number of terms in the Maclaurin series.
   * @tparam I  The current index in the recursion (should start from N - 1).
   *
   * The static compute function updates the coefficients `c` and `s` by
   * recursively calling itself and applying the corresponding double-angle
   * factors from the COS_MACLAURIN_DOUBLEANGLE_FACTOR and
   * SIN_MACLAURIN_DOUBLEANGLE_FACTOR arrays.
   *
   * @param[out] c  Reference to the cosine coefficient to be updated.
   * @param[out] s  Reference to the sine coefficient to be updated.
   * @param[in]  z  The variable used in the Maclaurin series expansion.
   */
  static void compute(T &c, T &s, const T &z) {
    FirstLoop<T, N, I - 1>::compute(c, s, z);
    c = c * z + static_cast<T>(COS_MACLAURIN_DOUBLEANGLE_FACTOR[N - 1 - I]);
    s = s * z + static_cast<T>(SIN_MACLAURIN_DOUBLEANGLE_FACTOR[N - 1 - I]);
  }
};

template <typename T, std::size_t N> struct FirstLoop<T, N, 0> {
  /**
   * @brief Computes the next terms in the Maclaurin series for cosine and sine
   * using the double-angle formula.
   *
   * This static function updates the values of `c` (cosine) and `s` (sine) by
   * multiplying them with `z` and adding the corresponding precomputed
   * Maclaurin series coefficients for the current term.
   *
   * @tparam T The numeric type used for computation (e.g., float, double).
   * @tparam N The current term index in the Maclaurin series.
   * @param[out] c Reference to the current cosine value, updated in-place.
   * @param[out] s Reference to the current sine value, updated in-place.
   * @param[in] z The variable used in the Maclaurin series expansion.
   */
  static void compute(T &c, T &s, const T &z) {
    c = c * z + static_cast<T>(COS_MACLAURIN_DOUBLEANGLE_FACTOR[N - 1]);
    s = s * z + static_cast<T>(SIN_MACLAURIN_DOUBLEANGLE_FACTOR[N - 1]);
  }
};

template <typename T, std::size_t N> struct SecondLoop {
  /**
   * @brief Performs a single iteration of a specialized trigonometric
   * computation on the provided values.
   *
   * This static function updates the values of \p s and \p c using the
   * following formulas:
   *   - s = s * (2 - c)
   *   - c = c * (4 - c)
   * After updating, it recursively calls SecondLoop<T, N - 1>::compute(s, c).
   *
   * @tparam T The numeric type of the input values (e.g., float, double).
   * @tparam N The current iteration count for recursion.
   * @param s Reference to the first value to be updated.
   * @param c Reference to the second value to be updated.
   */
  static void compute(T &s, T &c) {
    s = s * (static_cast<T>(2) - c);
    c = c * (static_cast<T>(4) - c);
    SecondLoop<T, N - 1>::compute(s, c);
  }
};

template <typename T> struct SecondLoop<T, 0> {
  /**
   * @brief Computes values for s and c.
   *
   * This static function is intended to compute or assign values to the
   * parameters `s` and `c`. In its current implementation, it does nothing and
   * simply suppresses unused parameter warnings.
   *
   * @tparam T The type of the parameters.
   * @param s Reference to a variable of type T, intended to receive a computed
   * value.
   * @param c Reference to a variable of type T, intended to receive a computed
   * value.
   */
  static void compute(T &s, T &c) {
    static_cast<void>(s);
    static_cast<void>(c);
    /* Do Nothing. */
  }
};

} // namespace SinCosMcLoughlinExpansion

/**
 * @brief Computes the sine and cosine of an angle using the Maclaurin series
 * expansion with the double angle formula for improved accuracy and efficiency.
 *
 * This function calculates the sine and cosine values of the input angle
 * `theta` using a Maclaurin series expansion, leveraging the double angle
 * formula to enhance numerical stability and convergence. The number of terms
 * in the expansion is controlled by the template parameter
 * `MACLAURIN_EXPANSION_REPEAT_NUMBER`.
 *
 * @tparam T The floating-point type for computation (e.g., float, double).
 * @tparam MACLAURIN_EXPANSION_REPEAT_NUMBER The number of terms to use in the
 * Maclaurin expansion. Must be less than
 * `COS_MACLAURIN_DOUBLEANGLE_FACTOR_MAX_SIZE`.
 * @param[in]  theta      The input angle (in radians).
 * @param[out] sin_value  The computed sine value of `theta`.
 * @param[out] cos_value  The computed cosine value of `theta`.
 *
 * @note The input angle is first wrapped to the range [-π, π] for improved
 * accuracy.
 * @note The function uses precomputed Maclaurin expansion factors and
 * factorials.
 * @note The function relies on helper classes and constants such as
 *       `Base::Math::wrap_value_in_minus_pi_and_pi`,
 * `COS_MACLAURIN_DOUBLEANGLE_FACTOR`, `SIN_MACLAURIN_DOUBLEANGLE_FACTOR`,
 * `Factorial_2`, and `SinCosMcLoughlinExpansion::FirstLoop`/`SecondLoop`.
 * @throws static_assert if `MACLAURIN_EXPANSION_REPEAT_NUMBER` is too large.
 */
template <typename T, std::size_t MACLAURIN_EXPANSION_REPEAT_NUMBER>
inline void sincos_maclaurin_expansion_with_DoubleAngleFormula(const T &theta,
                                                               T &sin_value,
                                                               T &cos_value) {
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

  SinCosMcLoughlinExpansion::FirstLoop<T, MACLAURIN_EXPANSION_REPEAT_NUMBER,
                                       (MACLAURIN_EXPANSION_REPEAT_NUMBER -
                                        1)>::compute(c, s, z);

  c = c * z;
  s = s * theta_wrapped;

  SinCosMcLoughlinExpansion::SecondLoop<
      T, MACLAURIN_EXPANSION_REPEAT_NUMBER>::compute(s, c);

  cos_value = static_cast<T>(1) - c * static_cast<T>(0.5);
  sin_value = s;
}

/* sin maclaurin expansion with DoubleAngleFormula */

/**
 * @brief Computes the sine of a value using the Maclaurin series expansion and
 * the double angle formula.
 *
 * This function calculates sin(x) by leveraging the cosine Maclaurin expansion
 * with the double angle formula, applying the identity sin(x) = cos(x - π/2).
 * The number of terms used in the Maclaurin expansion is controlled by the
 * template parameter MACLAURIN_EXPANSION_REPEAT_NUMBER.
 *
 * @tparam T The floating-point type for computation (e.g., float, double).
 * @tparam MACLAURIN_EXPANSION_REPEAT_NUMBER The number of terms to use in the
 * Maclaurin series expansion.
 * @param x The input value (in radians) for which to compute the sine.
 * @return The approximate value of sin(x).
 */
template <typename T, std::size_t MACLAURIN_EXPANSION_REPEAT_NUMBER>
inline T sin_maclaurin_expansion_with_DoubleAngleFormula(const T &x) {

  return Base::Math::cos_maclaurin_expansion_with_DoubleAngleFormula<
      T, MACLAURIN_EXPANSION_REPEAT_NUMBER>(
      x - static_cast<T>(Base::Math::HALF_PI));
}

/* tan maclaurin expansion with DoubleAngleFormula */

/**
 * @brief Computes the tangent of a value using the Maclaurin series expansion
 * and the double angle formula.
 *
 * This function calculates the tangent of the input value `x` by first
 * computing the sine and cosine using a Maclaurin series expansion with the
 * double angle formula for improved accuracy and efficiency. The tangent is
 * then obtained as the ratio of sine to cosine, with a safeguard to avoid
 * division by zero.
 *
 * @tparam T The floating-point type (e.g., float, double) for computation.
 * @tparam MACLAURIN_EXPANSION_REPEAT_NUMBER The number of terms to use in the
 * Maclaurin series expansion.
 * @param x The input value (in radians) for which to compute the tangent.
 * @return The tangent of the input value `x`.
 */
template <typename T, std::size_t MACLAURIN_EXPANSION_REPEAT_NUMBER>
inline T tan_maclaurin_expansion_with_DoubleAngleFormula(const T &x) {
  T cos_value = static_cast<T>(0);
  T sin_value = static_cast<T>(0);

  Base::Math::sincos_maclaurin_expansion_with_DoubleAngleFormula<
      T, MACLAURIN_EXPANSION_REPEAT_NUMBER>(x, sin_value, cos_value);

  return sin_value /
         Base::Utility::avoid_zero_divide(
             cos_value, static_cast<T>(Base::Math::TRIGONOMETRIC_DIVISION_MIN));
}

/* sincos */

/**
 * @brief Computes the sine and cosine of the given value.
 *
 * This function calculates both the sine and cosine of the input value `x` and
 * stores the results in `sin_value` and `cos_value`, respectively. The
 * implementation can use either the standard library functions (`std::sin`,
 * `std::cos`) or a custom Maclaurin expansion with the double angle formula,
 * depending on whether `__BASE_MATH_USE_STD_MATH__` is defined.
 *
 * @tparam T The numeric type of the input and output values.
 * @param x The input value (angle in radians) for which to compute sine and
 * cosine.
 * @param sin_value Reference to a variable where the computed sine value will
 * be stored.
 * @param cos_value Reference to a variable where the computed cosine value will
 * be stored.
 */
template <typename T>
inline void sincos(const T &x, T &sin_value, T &cos_value) {

#ifdef __BASE_MATH_USE_STD_MATH__
  sin_value = std::sin(x);
  cos_value = std::cos(x);
#else // __BASE_MATH_USE_STD_MATH__

  Base::Math::sincos_maclaurin_expansion_with_DoubleAngleFormula<
      T, Base::Math::SINCOS_MACLAURIN_DOUBLEANGLE_REPEAT_NUMBER>(x, sin_value,
                                                                 cos_value);

#endif // __BASE_MATH_USE_STD_MATH__
}

/* sin maclaurin expansion */
namespace SinMaclaurin {

template <typename T, std::size_t LOOP_MAX, std::size_t N> struct Loop {
  /**
   * @brief Computes the next term in a trigonometric series expansion and
   * updates the result.
   *
   * This static function multiplies the current term by a coefficient based on
   * the squared input value and the current loop index, then adds the term to
   * the result. It recursively calls itself to process the next term in the
   * series.
   *
   * @tparam T        Numeric type used for computation (e.g., float, double).
   * @tparam LOOP_MAX Maximum number of terms in the series expansion.
   * @tparam N        Current term index (decremented with each recursion).
   * @param x_squared The squared value of the input variable (e.g., x^2).
   * @param term      Reference to the current term in the series (updated
   * in-place).
   * @param result    Reference to the accumulated result of the series (updated
   * in-place).
   */
  static void compute(const T &x_squared, T &term, T &result) {
    term *=
        -x_squared * static_cast<T>(static_cast<T>(1) /
                                    static_cast<T>((2 * (LOOP_MAX - N)) *
                                                   (2 * (LOOP_MAX - N) + 1)));
    result += term;

    Loop<T, LOOP_MAX, N - 1>::compute(x_squared, term, result);
  }
};

template <typename T, std::size_t LOOP_MAX> struct Loop<T, LOOP_MAX, 0> {
  /**
   * @brief Computes a value based on the squared input, updating the term and
   * result.
   *
   * This is a placeholder implementation that performs no computation.
   * The function is intended to be specialized or overridden for specific
   * behavior.
   *
   * @tparam T The numeric type of the input and output parameters.
   * @param x_squared The squared value of the input variable.
   * @param term Reference to the term variable to be updated.
   * @param result Reference to the result variable to be updated.
   */
  static void compute(const T &x_squared, T &term, T &result) {
    /* Do Nothing. */
    static_cast<void>(x_squared);
    static_cast<void>(term);
    static_cast<void>(result);
  }
};

} // namespace SinMaclaurin

/**
 * @brief Computes the sine of a value using the Maclaurin series expansion.
 *
 * This function calculates an approximation of the sine of the input value `x`
 * using a Maclaurin series expansion up to a specified number of terms
 * (`LOOP_NUMBER`). The input value is first wrapped to the range [-π, π] for
 * improved numerical stability.
 *
 * @tparam T The floating-point type of the input and result (e.g., float,
 * double).
 * @tparam LOOP_NUMBER The number of terms to use in the Maclaurin series
 * expansion.
 * @param x The input value (in radians) for which to compute the sine.
 * @return The approximate sine of the input value.
 */
template <typename T, std::size_t LOOP_NUMBER>
inline T sin_maclaurin_expansion(const T &x) {

  T x_wrapped = Base::Math::wrap_value_in_minus_pi_and_pi(x);

  T term = x_wrapped;
  T result = x_wrapped;
  T x_squared = x_wrapped * x_wrapped;

  SinMaclaurin::Loop<T, LOOP_NUMBER, LOOP_NUMBER - 1>::compute(x_squared, term,
                                                               result);

  return result;
}

/* sin */

/**
 * @brief Computes the sine of the given value.
 *
 * This function calculates the sine of the input value `x`. If the macro
 * `__BASE_MATH_USE_STD_MATH__` is defined, it uses the standard library's
 * `std::sin` function. Otherwise, it computes the sine using a Maclaurin
 * series expansion with the double angle formula, as implemented in
 * `Base::Math::sin_maclaurin_expansion_with_DoubleAngleFormula`.
 *
 * @tparam T The type of the input value (e.g., float, double).
 * @param x The value (in radians) for which to compute the sine.
 * @return The sine of the input value.
 */
template <typename T> inline T sin(const T &x) {

#ifdef __BASE_MATH_USE_STD_MATH__
  return std::sin(x);
#else // __BASE_MATH_USE_STD_MATH__

  return Base::Math::sin_maclaurin_expansion_with_DoubleAngleFormula<
      T, Base::Math::SIN_MACLAURIN_DOUBLEANGLE_REPEAT_NUMBER>(x);

#endif // __BASE_MATH_USE_STD_MATH__
}

/* cos maclaurin expansion */
namespace CosMaclaurin {

template <typename T, std::size_t LOOP_MAX, std::size_t N> struct Loop {
  /**
   * @brief Computes and accumulates a term in a trigonometric series expansion.
   *
   * This static function updates the current term and result for a
   * trigonometric series (such as sine or cosine) using the provided squared
   * input value. The function multiplies the current term by a calculated
   * coefficient based on the template parameters and adds it to the result. It
   * then recursively calls the next iteration of the loop.
   *
   * @tparam T        Numeric type used for computation (e.g., float, double).
   * @tparam LOOP_MAX Maximum number of terms in the series expansion.
   * @tparam N        Current iteration index (decremented in recursion).
   * @param x_squared The squared value of the input variable (e.g., x^2).
   * @param term      Reference to the current term in the series (updated in
   * place).
   * @param result    Reference to the accumulated result of the series (updated
   * in place).
   */
  static void compute(const T &x_squared, T &term, T &result) {
    term *=
        -x_squared * static_cast<T>(static_cast<T>(1) /
                                    static_cast<T>((2 * (LOOP_MAX - N) - 1) *
                                                   (2 * (LOOP_MAX - N))));
    result += term;

    Loop<T, LOOP_MAX, N - 1>::compute(x_squared, term, result);
  }
};

template <typename T, std::size_t LOOP_MAX> struct Loop<T, LOOP_MAX, 0> {
  /**
   * @brief Computes a value based on the squared input, updating the term and
   * result.
   *
   * This is a placeholder implementation that performs no computation.
   * The function is intended to be specialized or overridden for specific
   * behavior.
   *
   * @tparam T The numeric type of the input and output parameters.
   * @param x_squared The squared value of the input variable.
   * @param term Reference to the term variable to be updated.
   * @param result Reference to the result variable to be updated.
   */
  static void compute(const T &x_squared, T &term, T &result) {
    /* Do Nothing. */
    static_cast<void>(x_squared);
    static_cast<void>(term);
    static_cast<void>(result);
  }
};

} // namespace CosMaclaurin

/**
 * @brief Computes the cosine of a value using the Maclaurin series expansion.
 *
 * This function evaluates the cosine of the input value `x` by wrapping it into
 * the range [-π, π] and then applying the Maclaurin series expansion up to the
 * specified number of terms (`LOOP_NUMBER`). The computation is performed using
 * a compile-time loop unrolling technique for efficiency.
 *
 * @tparam T The floating-point type to use for computation (e.g., float,
 * double).
 * @tparam LOOP_NUMBER The number of terms to use in the Maclaurin series
 * expansion.
 * @param x The input value (in radians) for which to compute the cosine.
 * @return The approximate value of cos(x) computed using the Maclaurin series.
 */
template <typename T, std::size_t LOOP_NUMBER>
inline T cos_maclaurin_expansion(const T &x) {

  T x_wrapped = wrap_value_in_minus_pi_and_pi(x);

  T term = static_cast<T>(1);
  T result = static_cast<T>(1);
  T x_squared = x_wrapped * x_wrapped;

  CosMaclaurin::Loop<T, LOOP_NUMBER, LOOP_NUMBER - 1>::compute(x_squared, term,
                                                               result);

  return result;
}

/* cos */

/**
 * @brief Computes the cosine of the given value.
 *
 * This function computes the cosine of the input value `x`. Depending on the
 * compilation settings, it either uses the standard library's `std::cos`
 * function or a custom implementation based on the Maclaurin series expansion
 * with the double angle formula.
 *
 * @tparam T The type of the input value (e.g., float, double).
 * @param x The value (in radians) for which to compute the cosine.
 * @return The cosine of the input value `x`.
 */
template <typename T> inline T cos(const T &x) {

#ifdef __BASE_MATH_USE_STD_MATH__
  return std::cos(x);
#else // __BASE_MATH_USE_STD_MATH__

  return Base::Math::cos_maclaurin_expansion_with_DoubleAngleFormula<
      T, Base::Math::COS_MACLAURIN_DOUBLEANGLE_REPEAT_NUMBER>(x);

#endif // __BASE_MATH_USE_STD_MATH__
}

/* tan maclaurin expansion */

/**
 * @brief Computes the tangent of x using the Maclaurin series expansions for
 * sine and cosine.
 *
 * This function calculates tan(x) by dividing the Maclaurin series expansion of
 * sine by that of cosine. To prevent division by zero, the denominator is
 * passed through a utility function that ensures it is not too close to zero.
 *
 * @tparam T Numeric type (e.g., float, double).
 * @tparam SIN_LOOP_NUMBER Number of terms to use in the sine Maclaurin
 * expansion.
 * @tparam COS_LOOP_NUMBER Number of terms to use in the cosine Maclaurin
 * expansion.
 * @param x The input angle in radians.
 * @return The approximate value of tan(x) computed via Maclaurin series.
 */
template <typename T, std::size_t SIN_LOOP_NUMBER, std::size_t COS_LOOP_NUMBER>
inline T tan_maclaurin_expansion(const T &x) {

  return Base::Math::sin_maclaurin_expansion<T, SIN_LOOP_NUMBER>(x) /
         Base::Utility::avoid_zero_divide(
             Base::Math::cos_maclaurin_expansion<T, COS_LOOP_NUMBER>(x),
             static_cast<T>(Base::Math::TRIGONOMETRIC_DIVISION_MIN));
}

/* tan */

/**
 * @brief Computes the tangent of the given value.
 *
 * This function calculates the tangent of the input value `x` of type `T`.
 * If the macro `__BASE_MATH_USE_STD_MATH__` is defined, it uses the standard
 * library's `std::tan` function. Otherwise, it computes the tangent using
 * a Maclaurin series expansion with a double angle formula, as implemented in
 * `Base::Math::tan_maclaurin_expansion_with_DoubleAngleFormula`.
 *
 * @tparam T The numeric type of the input value.
 * @param x The value (in radians) for which to compute the tangent.
 * @return The tangent of `x`.
 */
template <typename T> inline T tan(const T &x) {

#ifdef __BASE_MATH_USE_STD_MATH__
  return std::tan(x);
#else // __BASE_MATH_USE_STD_MATH__

  return Base::Math::tan_maclaurin_expansion_with_DoubleAngleFormula<
      T, Base::Math::SINCOS_MACLAURIN_DOUBLEANGLE_REPEAT_NUMBER>(x);

#endif // __BASE_MATH_USE_STD_MATH__
}

/* atan maclaurin expansion */

/**
 * @brief Computes the arctangent (atan) of a value using the Maclaurin series
 * expansion.
 *
 * This function calculates the arctangent of the input value `x` using a
 * Maclaurin series expansion, with optimizations for values of `x` outside the
 * interval [-1, 1] and for values with large magnitude. The number of terms in
 * the series is controlled by the template parameter `LOOP_NUMBER`.
 *
 * Special handling is performed for:
 * - |x| > 1: Uses the identity atan(x) = sign(x) * pi/2 - atan(1/x) for better
 * convergence.
 * - |x| > 0.5: Splits the computation to improve convergence using the double
 * angle formula.
 * - |x| <= 0.5: Directly applies the Maclaurin series for atan(x).
 *
 * @tparam T Numeric type (e.g., float, double).
 * @tparam LOOP_NUMBER Number of terms to use in the Maclaurin series expansion.
 * @param x The value for which to compute the arctangent.
 * @return The computed arctangent of `x`.
 */
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
namespace AtanChebyshev {

template <typename T, int N> struct Loop {
  /**
   * @brief Recursively computes a value using Chebyshev coefficients for
   * arctangent approximation.
   *
   * This static function applies a recursive loop to compute a value based on
   * the Chebyshev polynomial coefficients for the arctangent function. It takes
   * the squared input value and an accumulator, updating the accumulator at
   * each recursion step.
   *
   * @tparam T The numeric type used for computation.
   * @tparam N The current recursion depth or index for the Chebyshev
   * coefficient.
   * @param x_squared The squared value of the input variable (typically x^2).
   * @param y The accumulator value, updated at each recursion step.
   * @return The computed value after applying the recursive Chebyshev
   * polynomial evaluation.
   */
  static T compute(const T &x_squared, const T &y) {
    return Loop<T, N - 1>::compute(
        x_squared,
        y * x_squared +
            static_cast<T>(Base::Math::CHEBYSHEV_COEFFICIENT_FOR_ATAN[N]));
  }
};

template <typename T> struct Loop<T, 0> {
  /**
   * @brief Computes a value based on the given squared input and another
   * parameter using a Chebyshev coefficient.
   *
   * This static function calculates the result of the expression: y * x_squared
   * + CHEBYSHEV_COEFFICIENT_FOR_ATAN[0], where CHEBYSHEV_COEFFICIENT_FOR_ATAN
   * is a predefined array of coefficients used for Chebyshev polynomial
   * approximation of the arctangent function.
   *
   * @tparam T The numeric type of the input and output values.
   * @param x_squared The squared value of the input variable (typically x^2).
   * @param y An additional parameter to be multiplied with x_squared.
   * @return The computed value as per the described formula.
   */
  static T compute(const T &x_squared, const T &y) {
    return y * x_squared +
           static_cast<T>(Base::Math::CHEBYSHEV_COEFFICIENT_FOR_ATAN[0]);
  }
};

} // namespace AtanChebyshev

/**
 * @brief Computes an approximation of the arctangent function using a Chebyshev
 * polynomial expansion.
 *
 * This function evaluates the arctangent (atan) of the input value `x` by
 * applying a Chebyshev polynomial approximation. The computation is performed
 * using a recursive loop structure defined in the `AtanChebyshev::Loop` class
 * template, which iterates `LOOP_NUMBER - 1` times to accumulate the result.
 *
 * @tparam T           The numeric type of the input and output (e.g., float,
 * double).
 * @tparam LOOP_NUMBER The number of terms (iterations) to use in the Chebyshev
 * polynomial expansion.
 * @param x            The input value for which to compute the arctangent
 * approximation.
 * @return T           The approximated arctangent of `x`.
 */
template <typename T, std::size_t LOOP_NUMBER>
inline T atan_chebyshev_core(const T &x) {

  T x_squared = x * x;
  T y = static_cast<T>(0);

  y = AtanChebyshev::Loop<T, LOOP_NUMBER - 1>::compute(x_squared, y);

  return y * x;
}

/**
 * @brief Computes an approximation of the arctangent (atan) function using a
 * Chebyshev polynomial expansion.
 *
 * This function improves the accuracy of the Chebyshev approximation for values
 * of x outside the interval [-0.5, 0.5] by recursively reducing the argument.
 * For |x| > 0.5, it splits the computation into two parts using the half-angle
 * identity to ensure better convergence of the Chebyshev series. For |x| <=
 * 0.5, it directly applies the Chebyshev core.
 *
 * @tparam T           The floating-point type (e.g., float, double).
 * @tparam LOOP_NUMBER The number of terms to use in the Chebyshev expansion.
 * @param x            The input value for which to compute atan(x).
 * @return T           The approximated arctangent of x.
 */
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

/**
 * @brief Computes the arctangent (atan) of a value using a Chebyshev polynomial
 * approximation.
 *
 * This function evaluates the arctangent of the input value `x` using a
 * Chebyshev polynomial expansion with a specified number of terms
 * (`LOOP_NUMBER`). The implementation handles wide input ranges by transforming
 * the argument when |x| > 1, ensuring numerical stability.
 *
 * @tparam T           The floating-point type (e.g., float, double).
 * @tparam LOOP_NUMBER The number of Chebyshev coefficients/terms to use in the
 * approximation. Must not exceed
 * Base::Math::CHEBYSHEV_COEFFICIENT_FOR_ATAN_SIZE.
 * @param x            The input value for which to compute the arctangent.
 * @return T           The approximated arctangent of `x`.
 *
 * @note For |x| > 1, the function uses the identity atan(x) = sign(x) * HALF_PI
 * - atan(1/x) to improve accuracy.
 * @note Requires that `Base::Math::atan_chebyshev_core_wide` and related
 * constants are defined.
 */
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

/**
 * @brief Computes the arctangent (inverse tangent) of the given value.
 *
 * This function returns the arctangent of the input value \p x.
 * Depending on the compilation flag `__BASE_MATH_USE_STD_MATH__`, it either:
 * - Uses the standard library implementation (`std::atan`), or
 * - Uses a custom Chebyshev polynomial approximation
 * (`Base::Math::atan_chebyshev`).
 *
 * @tparam T Numeric type of the input value.
 * @param x The value whose arctangent is to be computed.
 * @return The arctangent of \p x, in radians.
 */
template <typename T> inline T atan(const T &x) {

#ifdef __BASE_MATH_USE_STD_MATH__
  return std::atan(x);
#else // __BASE_MATH_USE_STD_MATH__

  return Base::Math::atan_chebyshev<T, Base::Math::ATAN_REPEAT_NUMBER>(x);

#endif // __BASE_MATH_USE_STD_MATH__
}

/* atan2 Chebyshev */

/**
 * @brief Computes the arctangent of y/x using a Chebyshev polynomial
 * approximation.
 *
 * This function provides an efficient and accurate approximation of the atan2
 * function, which computes the angle (in radians) between the positive x-axis
 * and the point (x, y). The approximation is performed using Chebyshev
 * polynomials, with the number of terms controlled by the template parameter
 * LOOP_NUMBER.
 *
 * Special cases are handled for near-zero values of x and y to avoid division
 * by zero and to ensure correct quadrant determination:
 *   - If both x and y are near zero, returns 0.
 *   - If x > 0, returns atan_chebyshev(y / x).
 *   - If x < 0 and y >= 0, returns atan_chebyshev(y / x) + PI.
 *   - If x < 0 and y < 0, returns atan_chebyshev(y / x) - PI.
 *   - If x is near zero and y > 0, returns HALF_PI.
 *   - If x is near zero and y < 0, returns -HALF_PI.
 *
 * @tparam T Numeric type (e.g., float, double).
 * @tparam LOOP_NUMBER Number of terms in the Chebyshev polynomial
 * approximation.
 * @param y The y-coordinate.
 * @param x The x-coordinate.
 * @return The angle in radians between the positive x-axis and the point (x,
 * y).
 */
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

/**
 * @brief Computes the arc tangent of y/x, using either the standard library or
 * a custom Chebyshev approximation.
 *
 * This function calculates the angle (in radians) whose tangent is the quotient
 * of two specified values (y/x). If the macro __BASE_MATH_USE_STD_MATH__ is
 * defined, it uses std::atan2 from the standard library. Otherwise, it uses a
 * custom Chebyshev polynomial approximation via Base::Math::atan2_chebyshev.
 *
 * @tparam T Numeric type of the input arguments.
 * @param y The ordinate coordinate.
 * @param x The abscissa coordinate.
 * @return The arc tangent of y/x, in radians.
 */
template <typename T> inline T atan2(const T &y, const T &x) {

#ifdef __BASE_MATH_USE_STD_MATH__
  return std::atan2(y, x);
#else // __BASE_MATH_USE_STD_MATH__

  return Base::Math::atan2_chebyshev<T, Base::Math::ATAN_REPEAT_NUMBER>(y, x);

#endif // __BASE_MATH_USE_STD_MATH__
}

/* asin Chebyshev */

/**
 * @brief Computes the arcsine (inverse sine) of a value using Chebyshev
 * polynomial approximation.
 *
 * This function calculates the arcsine of the input value `x` using a Chebyshev
 * polynomial-based approximation for improved numerical stability and
 * performance. It handles edge cases where `x` is outside the domain [-1, 1] by
 * clamping the result to ±π/2.
 *
 * @tparam T                Numeric type (e.g., float, double).
 * @tparam ATAN_LOOP_NUMBER Number of iterations or terms for the atan Chebyshev
 * approximation.
 * @tparam SQRT_LOOP_NUMBER Number of iterations for the Newton method in square
 * root calculation.
 * @param x                 The input value for which to compute arcsine.
 * @return T                The approximated arcsine of `x` in radians.
 *
 * @note Requires `Base::Math::atan_chebyshev` and
 * `Base::Math::sqrt_newton_method` implementations.
 */
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

/**
 * @brief Computes the arcsine (inverse sine) of the given value.
 *
 * This function returns the arcsine of the input value \p x.
 * If the macro __BASE_MATH_USE_STD_MATH__ is defined, it uses the standard
 * library's std::asin function for values in the range (-1, 1), and clamps
 * the result to ±HALF_PI for values outside this range.
 * Otherwise, it uses a Chebyshev polynomial approximation via
 * Base::Math::asin_chebyshev.
 *
 * @tparam T Numeric type of the input and output (e.g., float, double).
 * @param x The value whose arcsine is to be computed.
 * @return The arcsine of \p x, in radians.
 */
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

/**
 * @brief Computes the arc cosine (inverse cosine) of a value using Chebyshev
 * polynomial approximation.
 *
 * This function calculates the arc cosine of the input value `x` by utilizing
 * the relationship between acos and asin, specifically: acos(x) = π/2 -
 * asin(x). The asin is computed using a Chebyshev polynomial approximation,
 * which is controlled by the template parameters.
 *
 * @tparam T The floating-point type (e.g., float, double) for the computation.
 * @tparam ATAN_LOOP_NUMBER The number of iterations or terms used in the
 * arctangent approximation (used internally by asin_chebyshev).
 * @tparam SQRT_LOOP_NUMBER The number of iterations or terms used in the square
 * root approximation (used internally by asin_chebyshev).
 * @param x The input value for which to compute the arc cosine. Should be in
 * the range [-1, 1].
 * @return The arc cosine of `x`, in radians.
 */
template <typename T, std::size_t ATAN_LOOP_NUMBER,
          std::size_t SQRT_LOOP_NUMBER>
inline T acos_chebyshev(const T &x) {

  return static_cast<T>(Base::Math::HALF_PI) -
         Base::Math::asin_chebyshev<T, ATAN_LOOP_NUMBER, SQRT_LOOP_NUMBER>(x);
}

/* acos */

/**
 * @brief Computes the arc cosine (inverse cosine) of the given value.
 *
 * This function returns the principal value of the arc cosine of x, expressed
 * in radians. The implementation uses either the standard library's std::acos
 * or a custom Chebyshev polynomial approximation, depending on the compilation
 * flag __BASE_MATH_USE_STD_MATH__.
 *
 * Special cases:
 * - If x > 1, returns 0.
 * - If x < -1, returns π.
 *
 * @tparam T Numeric type of the input value.
 * @param x The value whose arc cosine is to be computed.
 * @return The arc cosine of x, in radians.
 */
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

/**
 * @brief Computes the hyperbolic sine (sinh) of a value using the Maclaurin
 * series expansion with a precomputed table for exponentials.
 *
 * This function calculates sinh(x) using the formula:
 *   sinh(x) = (exp(x) - exp(-x)) / 2
 * The exponential function is approximated using a Maclaurin series expansion
 * with a lookup table for improved performance.
 *
 * @tparam T The numeric type of the input and output (e.g., float, double).
 * @tparam LOOP_NUMBER The number of terms to use in the Maclaurin series
 * expansion.
 * @param x The value for which to compute the hyperbolic sine.
 * @return The approximate value of sinh(x).
 */
template <typename T, std::size_t LOOP_NUMBER>
inline T sinh_maclaurin_expansion_with_table(const T &x) {
  return (Base::Math::exp_maclaurin_expansion_with_table<T, LOOP_NUMBER>(x) -
          Base::Math::exp_maclaurin_expansion_with_table<T, LOOP_NUMBER>(-x)) *
         static_cast<T>(0.5);
}

/* sinh maclaurin expansion */

/**
 * @brief Computes the hyperbolic sine (sinh) of a value using the Maclaurin
 * series expansion for the exponential function.
 *
 * This function calculates sinh(x) by evaluating the difference between the
 * Maclaurin expansions of exp(x) and exp(-x), divided by 2. The number of terms
 * used in the expansion is controlled by LOOP_NUMBER.
 *
 * @tparam T The numeric type of the input and output (e.g., float, double).
 * @tparam LOOP_NUMBER The number of terms to use in the Maclaurin series
 * expansion.
 * @param x The value for which to compute the hyperbolic sine.
 * @return The approximate value of sinh(x) computed using the Maclaurin series.
 */
template <typename T, std::size_t LOOP_NUMBER>
inline T sinh_maclaurin_expansion(const T &x) {
  return (Base::Math::exp_maclaurin_expansion<T, LOOP_NUMBER>(x) -
          Base::Math::exp_maclaurin_expansion<T, LOOP_NUMBER>(-x)) *
         static_cast<T>(0.5);
}

/* sinh */

/**
 * @brief Computes the hyperbolic sine of the given value.
 *
 * This function calculates the hyperbolic sine (sinh) of the input value `x`.
 * Depending on the compilation flag `__BASE_MATH_USE_STD_MATH__`, it either
 * uses the standard library implementation (`std::sinh`) or a custom Maclaurin
 * expansion provided by `Base::Math::sinh_maclaurin_expansion`.
 *
 * @tparam T The numeric type of the input value.
 * @param x The value for which to compute the hyperbolic sine.
 * @return The hyperbolic sine of `x`.
 */
template <typename T> inline T sinh(const T &x) {

#ifdef __BASE_MATH_USE_STD_MATH__
  return std::sinh(x);
#else // __BASE_MATH_USE_STD_MATH__

  return Base::Math::sinh_maclaurin_expansion<T, Base::Math::EXP_REPEAT_NUMBER>(
      x);

#endif // __BASE_MATH_USE_STD_MATH__
}

/* cosh Maclaurin Expansion with table */

/**
 * @brief Computes the hyperbolic cosine (cosh) of a value using the Maclaurin
 * expansion for the exponential function.
 *
 * This function calculates cosh(x) by evaluating the Maclaurin series expansion
 * for exp(x) and exp(-x), then averaging the results as per the mathematical
 * definition: cosh(x) = (exp(x) + exp(-x)) / 2
 *
 * The expansion is computed up to LOOP_NUMBER terms for improved accuracy.
 * A lookup table may be used internally for optimization, depending on the
 * implementation of Base::Math::exp_maclaurin_expansion_with_table.
 *
 * @tparam T           The numeric type of the input and output (e.g., float,
 * double).
 * @tparam LOOP_NUMBER The number of terms to use in the Maclaurin expansion.
 * @param x            The value for which to compute the hyperbolic cosine.
 * @return T           The computed cosh(x) value.
 */
template <typename T, std::size_t LOOP_NUMBER>
inline T cosh_maclaurin_expansion_with_table(const T &x) {
  return (Base::Math::exp_maclaurin_expansion_with_table<T, LOOP_NUMBER>(x) +
          Base::Math::exp_maclaurin_expansion_with_table<T, LOOP_NUMBER>(-x)) *
         static_cast<T>(0.5);
}

/* cosh maclaurin expansion */

/**
 * @brief Computes the hyperbolic cosine (cosh) of a value using the Maclaurin
 * series expansion.
 *
 * This function calculates cosh(x) by evaluating the Maclaurin series expansion
 * for exp(x) and exp(-x), and then averaging the results as per the
 * mathematical definition: cosh(x) = (exp(x) + exp(-x)) / 2
 *
 * @tparam T The floating-point type of the input and output (e.g., float,
 * double).
 * @tparam LOOP_NUMBER The number of terms to use in the Maclaurin series
 * expansion for exp(x).
 * @param x The value for which to compute the hyperbolic cosine.
 * @return The approximate value of cosh(x) computed using the Maclaurin series.
 */
template <typename T, std::size_t LOOP_NUMBER>
inline T cosh_maclaurin_expansion(const T &x) {
  return (Base::Math::exp_maclaurin_expansion<T, LOOP_NUMBER>(x) +
          Base::Math::exp_maclaurin_expansion<T, LOOP_NUMBER>(-x)) *
         static_cast<T>(0.5);
}

/* cosh */

/**
 * @brief Computes the hyperbolic cosine of the given value.
 *
 * This function calculates the hyperbolic cosine (cosh) of the input value `x`.
 * If the macro `__BASE_MATH_USE_STD_MATH__` is defined, it uses the standard
 * library's `std::cosh` function. Otherwise, it computes the result using a
 * Maclaurin series expansion implemented in
 * `Base::Math::cosh_maclaurin_expansion`.
 *
 * @tparam T Numeric type of the input value.
 * @param x The value for which to compute the hyperbolic cosine.
 * @return The hyperbolic cosine of `x`.
 */
template <typename T> inline T cosh(const T &x) {

#ifdef __BASE_MATH_USE_STD_MATH__
  return std::cosh(x);
#else // __BASE_MATH_USE_STD_MATH__

  return Base::Math::cosh_maclaurin_expansion<T, Base::Math::EXP_REPEAT_NUMBER>(
      x);

#endif // __BASE_MATH_USE_STD_MATH__
}

/* tanh maclaurin expansion with table */

/**
 * @brief Computes the hyperbolic tangent (tanh) of a value using the Maclaurin
 * series expansion for the exponential function, with table-based optimization.
 *
 * This function calculates tanh(x) by evaluating the Maclaurin series expansion
 * for exp(x) and exp(-x) up to a specified number of terms (LOOP_NUMBER), and
 * then applies the tanh identity: tanh(x) = (exp(x) - exp(-x)) / (exp(x) +
 * exp(-x)). To ensure numerical stability, the denominator is passed through a
 * utility function to avoid division by zero.
 *
 * @tparam T Numeric type (e.g., float, double).
 * @tparam LOOP_NUMBER Number of terms to use in the Maclaurin series expansion.
 * @param x The input value for which to compute tanh(x).
 * @return The approximate value of tanh(x).
 */
template <typename T, std::size_t LOOP_NUMBER>
inline T tanh_maclaurin_expansion_with_table(const T &x) {

  T a = Base::Math::exp_maclaurin_expansion_with_table<T, LOOP_NUMBER>(x);
  T b = Base::Math::exp_maclaurin_expansion_with_table<T, LOOP_NUMBER>(-x);

  return (a - b) / Base::Utility::avoid_zero_divide(
                       a + b, static_cast<T>(TRIGONOMETRIC_DIVISION_MIN));
}

/* tanh maclaurin expansion */

/**
 * @brief Computes the hyperbolic tangent (tanh) of a value using the Maclaurin
 * series expansion for the exponential function.
 *
 * This function calculates tanh(x) by evaluating the Maclaurin series expansion
 * for exp(x) and exp(-x), then applying the identity: tanh(x) = (exp(x) -
 * exp(-x)) / (exp(x) + exp(-x)). To avoid division by zero, a utility function
 * is used to ensure the denominator is not too close to zero.
 *
 * @tparam T           The numeric type of the input and output (e.g., float,
 * double).
 * @tparam LOOP_NUMBER The number of terms to use in the Maclaurin series
 * expansion for exp(x).
 * @param x            The input value for which to compute tanh(x).
 * @return             The approximate value of tanh(x) computed using the
 * Maclaurin series.
 */
template <typename T, std::size_t LOOP_NUMBER>
inline T tanh_maclaurin_expansion(const T &x) {

  T a = Base::Math::exp_maclaurin_expansion<T, LOOP_NUMBER>(x);
  T b = Base::Math::exp_maclaurin_expansion<T, LOOP_NUMBER>(-x);

  return (a - b) / Base::Utility::avoid_zero_divide(
                       a + b, static_cast<T>(TRIGONOMETRIC_DIVISION_MIN));
}

/* tanh */

/**
 * @brief Computes the hyperbolic tangent of the given value.
 *
 * This function calculates the hyperbolic tangent (tanh) of the input value
 * `x`. If the macro `__BASE_MATH_USE_STD_MATH__` is defined, it uses the
 * standard library's `std::tanh` function. Otherwise, it computes the result
 * using a custom Maclaurin series expansion implementation.
 *
 * @tparam T The type of the input value (e.g., float, double).
 * @param x The value for which to compute the hyperbolic tangent.
 * @return The hyperbolic tangent of `x`.
 */
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
