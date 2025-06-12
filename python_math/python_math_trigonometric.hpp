/**
 * @file python_math_trigonometric.hpp
 * @brief Provides Python-like trigonometric and hyperbolic math functions for
 * scalars, std::vector, and std::array types.
 *
 * This header defines a collection of templated functions within the PythonMath
 * namespace that mimic the behavior of Python's math and numpy trigonometric
 * and hyperbolic functions. The functions are overloaded to operate on single
 * values, std::vector, and std::array, enabling element-wise operations similar
 * to Python's math and numpy modules. All functions delegate the actual
 * computation to the corresponding functions in the Base::Math namespace.
 *
 * @namespace PythonMath
 * @brief Contains templated trigonometric and hyperbolic math functions for
 * scalar and container types.
 *
 * The PythonMath namespace provides the following functions:
 * - sin, cos, tan, asin, acos, atan, atan2: Standard trigonometric functions.
 * - sinh, cosh, tanh: Hyperbolic trigonometric functions.
 * Each function is overloaded to support:
 *   - Scalar values (e.g., double, float)
 *   - std::vector<T>
 *   - std::array<T, N>
 * For atan2, overloads are provided for all combinations of scalar and
 * container arguments.
 */
#ifndef __PYTHON_MATH_TRIGONOMETRIC_HPP__
#define __PYTHON_MATH_TRIGONOMETRIC_HPP__

#include "base_math.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace PythonMath {

/* sin */

/**
 * @brief Computes the sine of the given value.
 *
 * This function is a wrapper around Base::Math::sin, providing a generic
 * interface for calculating the sine of a value of type T.
 *
 * @tparam T The numeric type of the input value.
 * @param x The value (in radians) for which to compute the sine.
 * @return The sine of the input value.
 */
template <typename T> inline T sin(const T &x) { return Base::Math::sin(x); }

/**
 * @brief Computes the sine of each element in the input vector.
 *
 * This function takes a vector of type T and returns a new vector where each
 * element is the result of applying Base::Math::sin to the corresponding
 * element of the input vector.
 *
 * @tparam T The numeric type of the vector elements.
 * @param vector The input vector containing elements to compute the sine of.
 * @return std::vector<T> A vector containing the sine of each input element.
 */
template <typename T> inline std::vector<T> sin(const std::vector<T> &vector) {
  std::vector<T> result;
  result.reserve(vector.size());
  for (const auto &element : vector) {
    result.push_back(Base::Math::sin(element));
  }
  return result;
}

/**
 * @brief Computes the element-wise sine of the input array.
 *
 * This function takes a std::array of type T and size N, and returns a new
 * array where each element is the sine of the corresponding element in the
 * input array. The sine computation is performed using Base::Math::sin for each
 * element.
 *
 * @tparam T The type of the elements in the array (e.g., float, double).
 * @tparam N The number of elements in the array.
 * @param array The input array whose elements' sines are to be computed.
 * @return std::array<T, N> An array containing the sine of each input element.
 */
template <typename T, std::size_t N>
inline std::array<T, N> sin(const std::array<T, N> &array) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::sin(array[i]);
  }
  return result;
}

/* cos */

/**
 * @brief Computes the cosine of the given value.
 *
 * This function is a wrapper that calls Base::Math::cos to compute the cosine
 * of the input value.
 *
 * @tparam T The type of the input value (e.g., float, double).
 * @param x The value (in radians) for which to compute the cosine.
 * @return The cosine of the input value.
 */
template <typename T> inline T cos(const T &x) { return Base::Math::cos(x); }

/**
 * @brief Computes the cosine of each element in the input vector.
 *
 * This function takes a vector of type T and returns a new vector where each
 * element is the result of applying the cosine function (Base::Math::cos) to
 * the corresponding element of the input vector.
 *
 * @tparam T The numeric type of the elements in the vector.
 * @param vector The input vector containing elements to compute the cosine of.
 * @return std::vector<T> A vector containing the cosine of each input element.
 */
template <typename T> inline std::vector<T> cos(const std::vector<T> &vector) {
  std::vector<T> result;
  result.reserve(vector.size());
  for (const auto &element : vector) {
    result.push_back(Base::Math::cos(element));
  }
  return result;
}

/**
 * @brief Computes the cosine of each element in the input array.
 *
 * This function takes a std::array of type T and size N, and returns a new
 * array where each element is the cosine of the corresponding element in the
 * input array. The cosine computation is performed using Base::Math::cos.
 *
 * @tparam T The type of the elements in the array (e.g., float, double).
 * @tparam N The size of the array.
 * @param array The input array whose elements' cosines are to be computed.
 * @return std::array<T, N> An array containing the cosine of each input
 * element.
 */
template <typename T, std::size_t N>
inline std::array<T, N> cos(const std::array<T, N> &array) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::cos(array[i]);
  }
  return result;
}

/* tan */

/**
 * @brief Computes the tangent of the given value.
 *
 * This function is a template wrapper that calls Base::Math::tan to compute the
 * tangent of the input value x.
 *
 * @tparam T The type of the input value (e.g., float, double).
 * @param x The value (in radians) for which to compute the tangent.
 * @return The tangent of x.
 */
template <typename T> inline T tan(const T &x) { return Base::Math::tan(x); }

/**
 * @brief Applies the tangent function to each element of the input vector.
 *
 * This function takes a vector of type T and returns a new vector where each
 * element is the result of applying Base::Math::tan to the corresponding
 * element of the input vector.
 *
 * @tparam T The type of the elements in the input vector.
 * @param vector The input vector whose elements will be transformed.
 * @return std::vector<T> A vector containing the tangent of each input element.
 */
template <typename T> inline std::vector<T> tan(const std::vector<T> &vector) {
  std::vector<T> result;
  result.reserve(vector.size());
  for (const auto &element : vector) {
    result.push_back(Base::Math::tan(element));
  }
  return result;
}

/**
 * @brief Applies the tangent function element-wise to a std::array.
 *
 * This function takes a std::array of type T and size N, computes the tangent
 * of each element using Base::Math::tan, and returns a new std::array
 * containing the results.
 *
 * @tparam T The type of the elements in the array.
 * @tparam N The size of the array.
 * @param array The input array whose elements will be transformed.
 * @return std::array<T, N> A new array where each element is the tangent of the
 * corresponding input element.
 */
template <typename T, std::size_t N>
inline std::array<T, N> tan(const std::array<T, N> &array) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::tan(array[i]);
  }
  return result;
}

/* atan */

/**
 * @brief Computes the arc tangent (inverse tangent) of the given value.
 *
 * This function returns the principal value of the arc tangent of x, expressed
 * in radians. It is a wrapper around Base::Math::atan.
 *
 * @tparam T Numeric type of the input value.
 * @param x The value whose arc tangent is to be computed.
 * @return The arc tangent of x, in radians.
 */
template <typename T> inline T atan(const T &x) { return Base::Math::atan(x); }

/**
 * @brief Computes the element-wise arctangent (inverse tangent) of a vector.
 *
 * This function takes a vector of type T and returns a new vector where each
 * element is the result of applying Base::Math::atan to the corresponding
 * element of the input vector.
 *
 * @tparam T The type of the elements in the input vector.
 * @param vector The input vector whose elements will be processed.
 * @return std::vector<T> A vector containing the arctangent of each input
 * element.
 */
template <typename T> inline std::vector<T> atan(const std::vector<T> &vector) {
  std::vector<T> result;
  result.reserve(vector.size());
  for (const auto &element : vector) {
    result.push_back(Base::Math::atan(element));
  }
  return result;
}

/**
 * @brief Applies the arctangent (atan) function element-wise to a std::array.
 *
 * This function takes a std::array of type T and size N, and returns a new
 * array where each element is the result of applying Base::Math::atan to the
 * corresponding element in the input array.
 *
 * @tparam T The type of the elements in the array.
 * @tparam N The size of the array.
 * @param array The input array whose elements will be transformed.
 * @return std::array<T, N> An array containing the arctangent of each input
 * element.
 */
template <typename T, std::size_t N>
inline std::array<T, N> atan(const std::array<T, N> &array) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::atan(array[i]);
  }
  return result;
}

/* atan2 */

/**
 * @brief Computes the arc tangent of y/x using the signs of both arguments to
 * determine the correct quadrant.
 *
 * This function returns the angle in radians between the positive x-axis and
 * the point (x, y). It is a wrapper around Base::Math::atan2, providing
 * type-generic support.
 *
 * @tparam T Numeric type of the input arguments (e.g., float, double).
 * @param y The y-coordinate.
 * @param x The x-coordinate.
 * @return The angle in radians between the positive x-axis and the point (x,
 * y).
 */
template <typename T> inline T atan2(const T &y, const T &x) {
  return Base::Math::atan2(y, x);
}

/**
 * @brief Computes the element-wise arc tangent of the quotient of each element
 * in the input vector and a scalar.
 *
 * This function applies the two-argument arctangent (atan2) operation to each
 * element of the input vector `vector_y` with respect to the scalar `x`,
 * returning a vector of results. The operation is equivalent to calling
 * `Base::Math::atan2(y_i, x)` for each element `y_i` in `vector_y`.
 *
 * @tparam T Numeric type of the vector elements and the scalar.
 * @param vector_y A vector of values representing the numerator in the atan2
 * operation.
 * @param x A scalar value representing the denominator in the atan2 operation.
 * @return std::vector<T> A vector containing the result of atan2 for each
 * element in `vector_y` with respect to `x`.
 */
template <typename T>
inline std::vector<T> atan2(const std::vector<T> &vector_y, const T &x) {
  std::vector<T> result;
  result.reserve(vector_y.size());
  for (const auto &element : vector_y) {
    result.push_back(Base::Math::atan2(element, x));
  }
  return result;
}

/**
 * @brief Computes the element-wise arc tangent of two variables (y, x) for a
 * scalar y and a vector x.
 *
 * This function takes a scalar value `y` and a vector of values `vector_x`, and
 * computes the arc tangent of each pair (y, x_i) using `Base::Math::atan2`. The
 * result is a vector containing the computed values.
 *
 * @tparam T The numeric type of the input and output values.
 * @param y The scalar value representing the numerator for the atan2
 * computation.
 * @param vector_x The vector of denominator values for the atan2 computation.
 * @return std::vector<T> A vector containing the result of atan2(y, x_i) for
 * each element x_i in `vector_x`.
 */
template <typename T>
inline std::vector<T> atan2(const T &y, const std::vector<T> &vector_x) {
  std::vector<T> result;
  result.reserve(vector_x.size());
  for (const auto &element : vector_x) {
    result.push_back(Base::Math::atan2(y, element));
  }
  return result;
}

/**
 * @brief Computes the element-wise arc tangent of two vectors (atan2).
 *
 * This function takes two vectors of the same size, `vector_y` and `vector_x`,
 * and computes the arc tangent of the quotient of their corresponding elements
 * using `Base::Math::atan2`. The result is a vector containing the computed
 * values for each pair of elements.
 *
 * @tparam T The numeric type of the input vectors (e.g., float, double).
 * @param vector_y The vector containing the y-coordinates.
 * @param vector_x The vector containing the x-coordinates.
 * @return std::vector<T> A vector containing the element-wise atan2 results.
 *
 * @note Both input vectors must have the same size.
 */
template <typename T>
inline std::vector<T> atan2(const std::vector<T> &vector_y,
                            const std::vector<T> &vector_x) {
  std::vector<T> result;
  result.reserve(vector_y.size());
  for (std::size_t i = 0; i < vector_y.size(); ++i) {
    result.push_back(Base::Math::atan2(vector_y[i], vector_x[i]));
  }
  return result;
}

/**
 * @brief Computes the element-wise arc tangent of each element in the input
 * array and a scalar value.
 *
 * This function applies the two-argument arctangent (atan2) operation to each
 * element of the input array `array_y` and the scalar value `x`, returning a
 * new array where each element is the result of `atan2(array_y[i], x)`.
 *
 * @tparam T The type of the elements in the input array and the scalar.
 * @tparam N The size of the input array.
 * @param array_y The input array containing the y-coordinates.
 * @param x The scalar x-coordinate to be used for all elements.
 * @return std::array<T, N> An array containing the result of `atan2(array_y[i],
 * x)` for each element.
 */
template <typename T, std::size_t N>
inline std::array<T, N> atan2(const std::array<T, N> &array_y, const T &x) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::atan2(array_y[i], x);
  }
  return result;
}

/**
 * @brief Computes the element-wise arc tangent of y and each element in the
 * input array x.
 *
 * This function takes a scalar value y and an array of values array_x, and
 * returns a new array where each element is the result of Base::Math::atan2(y,
 * array_x[i]). The atan2 function computes the angle (in radians) whose tangent
 * is the quotient of its arguments, handling the correct quadrant.
 *
 * @tparam T The type of the elements (e.g., float, double).
 * @tparam N The size of the input and output arrays.
 * @param y The scalar value to use as the first argument for atan2.
 * @param array_x The input array whose elements are used as the second argument
 * for atan2.
 * @return std::array<T, N> An array containing the result of atan2(y,
 * array_x[i]) for each element.
 */
template <typename T, std::size_t N>
inline std::array<T, N> atan2(const T &y, const std::array<T, N> &array_x) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::atan2(y, array_x[i]);
  }
  return result;
}

/**
 * @brief Computes the element-wise arc tangent of two arrays (atan2).
 *
 * This function takes two input arrays, `array_y` and `array_x`, each of size
 * `N`, and computes the arc tangent of the quotient of their corresponding
 * elements, storing the result in a new array. The computation is performed
 * using `Base::Math::atan2` for each element.
 *
 * @tparam T The type of the elements in the arrays (e.g., float, double).
 * @tparam N The size of the input arrays.
 * @param array_y The array containing the y-coordinates.
 * @param array_x The array containing the x-coordinates.
 * @return std::array<T, N> An array containing the element-wise atan2 results.
 */
template <typename T, std::size_t N>
inline std::array<T, N> atan2(const std::array<T, N> &array_y,
                              const std::array<T, N> &array_x) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::atan2(array_y[i], array_x[i]);
  }
  return result;
}

/* asin */

/**
 * @brief Computes the arc sine (inverse sine) of the given value.
 *
 * This function returns the angle whose sine is the specified value.
 * The result is expressed in radians and lies in the range [-π/2, π/2].
 *
 * @tparam T Numeric type of the input value.
 * @param x The value whose arc sine is to be computed.
 * @return The arc sine of x, in radians.
 *
 * @note The input value x should be in the range [-1, 1]. If x is outside this
 * range, the result is undefined.
 */
template <typename T> inline T asin(const T &x) { return Base::Math::asin(x); }

/**
 * @brief Computes the arc sine (inverse sine) of each element in the input
 * vector.
 *
 * This function takes a vector of elements and returns a new vector where each
 * element is the result of applying Base::Math::asin to the corresponding
 * element of the input vector.
 *
 * @tparam T The type of the elements in the input vector.
 * @param vector The input vector containing elements to compute the arc sine
 * for.
 * @return std::vector<T> A vector containing the arc sine of each input
 * element.
 */
template <typename T> inline std::vector<T> asin(const std::vector<T> &vector) {
  std::vector<T> result;
  result.reserve(vector.size());
  for (const auto &element : vector) {
    result.push_back(Base::Math::asin(element));
  }
  return result;
}

/**
 * @brief Computes the element-wise arcsine (inverse sine) of the input array.
 *
 * This function takes a std::array of type T and size N, and returns a new
 * array where each element is the result of applying Base::Math::asin to the
 * corresponding element of the input array.
 *
 * @tparam T The type of the elements in the array (e.g., float, double).
 * @tparam N The size of the array.
 * @param array The input array whose elements will be processed.
 * @return std::array<T, N> An array containing the arcsine of each input
 * element.
 */
template <typename T, std::size_t N>
inline std::array<T, N> asin(const std::array<T, N> &array) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::asin(array[i]);
  }
  return result;
}

/* acos */

/**
 * @brief Computes the arc cosine (inverse cosine) of the given value.
 *
 * This function returns the principal value of the arc cosine of x, expressed
 * in radians. The input value x must be in the range [-1, 1]. The result is in
 * the range [0, π].
 *
 * @tparam T Numeric type of the input value (e.g., float, double).
 * @param x The value whose arc cosine is to be computed.
 * @return The arc cosine of x, in radians.
 */
template <typename T> inline T acos(const T &x) { return Base::Math::acos(x); }

/**
 * @brief Computes the arc cosine (inverse cosine) of each element in the input
 * vector.
 *
 * This function takes a vector of elements and applies the arc cosine function
 * (acos) to each element using Base::Math::acos, returning a new vector with
 * the results.
 *
 * @tparam T The type of the elements in the input vector.
 * @param vector The input vector containing elements to compute the arc cosine
 * for.
 * @return std::vector<T> A vector containing the arc cosine of each input
 * element.
 */
template <typename T> inline std::vector<T> acos(const std::vector<T> &vector) {
  std::vector<T> result;
  result.reserve(vector.size());
  for (const auto &element : vector) {
    result.push_back(Base::Math::acos(element));
  }
  return result;
}

/**
 * @brief Computes the arc cosine (inverse cosine) of each element in the input
 * array.
 *
 * This function applies the Base::Math::acos function to each element of the
 * input std::array and returns a new std::array containing the results.
 *
 * @tparam T The type of the elements in the array (e.g., float, double).
 * @tparam N The size of the array.
 * @param array The input array whose elements will be processed.
 * @return std::array<T, N> An array containing the arc cosine of each input
 * element.
 */
template <typename T, std::size_t N>
inline std::array<T, N> acos(const std::array<T, N> &array) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::acos(array[i]);
  }
  return result;
}

/* sinh */

/**
 * @brief Computes the hyperbolic sine of the given value.
 *
 * This function is a wrapper that calls Base::Math::sinh to compute the
 * hyperbolic sine (sinh) of the input value.
 *
 * @tparam T The type of the input value (e.g., float, double).
 * @param x The value for which to compute the hyperbolic sine.
 * @return The hyperbolic sine of x.
 */
template <typename T> inline T sinh(const T &x) { return Base::Math::sinh(x); }

/**
 * @brief Computes the hyperbolic sine (sinh) of each element in the input
 * vector.
 *
 * This function applies the Base::Math::sinh operation to every element of the
 * input vector and returns a new vector containing the results.
 *
 * @tparam T The type of the elements in the input vector.
 * @param vector The input vector whose elements will be processed.
 * @return std::vector<T> A vector containing the hyperbolic sine of each input
 * element.
 */
template <typename T> inline std::vector<T> sinh(const std::vector<T> &vector) {
  std::vector<T> result;
  result.reserve(vector.size());
  for (const auto &element : vector) {
    result.push_back(Base::Math::sinh(element));
  }
  return result;
}

/**
 * @brief Computes the hyperbolic sine (sinh) of each element in the input
 * array.
 *
 * This function applies the Base::Math::sinh operation to each element of the
 * input std::array and returns a new array containing the results.
 *
 * @tparam T The type of the elements in the array.
 * @tparam N The size of the array.
 * @param array The input array whose elements will be transformed.
 * @return std::array<T, N> An array where each element is the hyperbolic sine
 * of the corresponding input element.
 */
template <typename T, std::size_t N>
inline std::array<T, N> sinh(const std::array<T, N> &array) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::sinh(array[i]);
  }
  return result;
}

/* cosh */

/**
 * @brief Computes the hyperbolic cosine of the given value.
 *
 * This function template forwards the computation to Base::Math::cosh.
 *
 * @tparam T The type of the input value (e.g., float, double).
 * @param x The value for which to compute the hyperbolic cosine.
 * @return The hyperbolic cosine of x.
 */
template <typename T> inline T cosh(const T &x) { return Base::Math::cosh(x); }

/**
 * @brief Computes the hyperbolic cosine (cosh) of each element in the input
 * vector.
 *
 * This function applies the Base::Math::cosh function to each element of the
 * input std::vector<T> and returns a new vector containing the results.
 *
 * @tparam T The type of the elements in the input vector.
 * @param vector The input vector whose elements will be transformed.
 * @return std::vector<T> A vector containing the hyperbolic cosine of each
 * input element.
 */
template <typename T> inline std::vector<T> cosh(const std::vector<T> &vector) {
  std::vector<T> result;
  result.reserve(vector.size());
  for (const auto &element : vector) {
    result.push_back(Base::Math::cosh(element));
  }
  return result;
}

/**
 * @brief Computes the hyperbolic cosine (cosh) of each element in the input
 * array.
 *
 * This function applies the Base::Math::cosh function to each element of the
 * input std::array, returning a new std::array containing the results.
 *
 * @tparam T The type of the elements in the array (e.g., float, double).
 * @tparam N The number of elements in the array.
 * @param array The input array whose elements will be transformed.
 * @return std::array<T, N> An array where each element is the hyperbolic cosine
 * of the corresponding input element.
 */
template <typename T, std::size_t N>
inline std::array<T, N> cosh(const std::array<T, N> &array) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::cosh(array[i]);
  }
  return result;
}

/* tanh */

/**
 * @brief Computes the hyperbolic tangent of the given value.
 *
 * This function is a wrapper around Base::Math::tanh, providing a generic
 * interface for computing the hyperbolic tangent for any type T supported by
 * Base::Math::tanh.
 *
 * @tparam T The type of the input value (e.g., float, double).
 * @param x The value for which to compute the hyperbolic tangent.
 * @return The hyperbolic tangent of x.
 */
template <typename T> inline T tanh(const T &x) { return Base::Math::tanh(x); }

/**
 * @brief Applies the hyperbolic tangent function element-wise to a vector.
 *
 * This function computes the hyperbolic tangent (tanh) of each element in the
 * input vector and returns a new vector containing the results.
 *
 * @tparam T The type of the elements in the input vector.
 * @param vector The input vector whose elements will be transformed.
 * @return std::vector<T> A vector containing the tanh of each input element.
 */
template <typename T> inline std::vector<T> tanh(const std::vector<T> &vector) {
  std::vector<T> result;
  result.reserve(vector.size());
  for (const auto &element : vector) {
    result.push_back(Base::Math::tanh(element));
  }
  return result;
}

/**
 * @brief Applies the hyperbolic tangent function element-wise to a std::array.
 *
 * This function takes a std::array of type T and size N, and returns a new
 * std::array where each element is the result of applying Base::Math::tanh to
 * the corresponding element in the input array.
 *
 * @tparam T The type of the elements in the array.
 * @tparam N The number of elements in the array.
 * @param array The input std::array whose elements will be transformed.
 * @return std::array<T, N> A new array with the hyperbolic tangent applied to
 * each element.
 */
template <typename T, std::size_t N>
inline std::array<T, N> tanh(const std::array<T, N> &array) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::tanh(array[i]);
  }
  return result;
}

} // namespace PythonMath

#endif // __PYTHON_MATH_TRIGONOMETRIC_HPP__
