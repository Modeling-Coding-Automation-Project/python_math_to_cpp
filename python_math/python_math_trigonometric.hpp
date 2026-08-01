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
#ifndef PYTHON_MATH_TRIGONOMETRIC_HPP_
#define PYTHON_MATH_TRIGONOMETRIC_HPP_

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
 * element is the result of applying PythonMath::sin to the corresponding
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
    result.push_back(PythonMath::sin(element));
  }
  return result;
}

namespace SinAction {

template <typename T, std::size_t N, std::size_t Index> struct SinCore {
  /**
   * @brief Recursively computes the element-wise sine of a std::array.
   *
   * @param result The array to store the results.
   * @param array The input array whose elements' sines are to be computed.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array) {
    result[Index] = PythonMath::sin(array[Index]);
    SinCore<T, N, Index - 1>::compute(result, array);
  }
};

// Specialization to end the recursion
template <typename T, std::size_t N> struct SinCore<T, N, 0> {
  /**
   * @brief Base case for the recursive computation of sines in a std::array.
   *
   * @param result The array to store the results.
   * @param array The input array whose elements' sines are to be computed.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array) {
    result[0] = PythonMath::sin(array[0]);
  }
};

/**
 * @brief Computes the element-wise sine of a std::array.
 *
 * This function serves as an entry point for computing the sines of elements
 * in a std::array. It initializes the recursive computation by calling SinCore
 * with the appropriate parameters.
 *
 * @tparam T The type of the elements in the array.
 * @tparam N The size of the array.
 * @param result The array to store the results.
 * @param array The input array whose elements' sines are to be computed.
 */
template <typename T, std::size_t N>
inline void compute(std::array<T, N> &result, const std::array<T, N> &array) {
  SinCore<T, N, N - 1>::compute(result, array);
}

} // namespace SinAction

/**
 * @brief Computes the element-wise sine of the input array.
 *
 * This function takes a std::array of type T and size N, and returns a new
 * array where each element is the sine of the corresponding element in the
 * input array. The sine computation is performed using PythonMath::sin for each
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
  SinAction::compute(result, array);
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
 * element is the result of applying the cosine function (PythonMath::cos) to
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
    result.push_back(PythonMath::cos(element));
  }
  return result;
}

namespace CosAction {

template <typename T, std::size_t N, std::size_t Index> struct CosCore {
  /**
   * @brief Recursively computes the element-wise cosine of a std::array.
   *
   * @param result The array to store the results.
   * @param array The input array whose elements' cosines are to be computed.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array) {
    result[Index] = PythonMath::cos(array[Index]);
    CosCore<T, N, Index - 1>::compute(result, array);
  }
};

// Specialization to end the recursion
template <typename T, std::size_t N> struct CosCore<T, N, 0> {
  /**
   * @brief Base case for the recursive computation of cosines in a std::array.
   *
   * @param result The array to store the results.
   * @param array The input array whose elements' cosines are to be computed.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array) {
    result[0] = PythonMath::cos(array[0]);
  }
};

/**
 * @brief Computes the element-wise cosine of a std::array.
 *
 * This function serves as an entry point for computing the cosines of elements
 * in a std::array. It initializes the recursive computation by calling CosCore
 * with the appropriate parameters.
 *
 * @tparam T The type of the elements in the array.
 * @tparam N The size of the array.
 * @param result The array to store the results.
 * @param array The input array whose elements' cosines are to be computed.
 */
template <typename T, std::size_t N>
inline void compute(std::array<T, N> &result, const std::array<T, N> &array) {
  CosCore<T, N, N - 1>::compute(result, array);
}

} // namespace CosAction

/**
 * @brief Computes the cosine of each element in the input array.
 *
 * This function takes a std::array of type T and size N, and returns a new
 * array where each element is the cosine of the corresponding element in the
 * input array. The cosine computation is performed using PythonMath::cos.
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
  CosAction::compute(result, array);
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
 * element is the result of applying PythonMath::tan to the corresponding
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
    result.push_back(PythonMath::tan(element));
  }
  return result;
}

namespace TanAction {

template <typename T, std::size_t N, std::size_t Index> struct TanCore {
  /**
   * @brief Recursively computes the element-wise tangent of a std::array.
   *
   * @param result The array to store the results.
   * @param array The input array whose elements will be transformed.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array) {
    result[Index] = PythonMath::tan(array[Index]);
    TanCore<T, N, Index - 1>::compute(result, array);
  }
};

// Specialization to end the recursion
template <typename T, std::size_t N> struct TanCore<T, N, 0> {
  /**
   * @brief Base case for the recursive computation of tangents in a std::array.
   *
   * @param result The array to store the results.
   * @param array The input array whose elements will be transformed.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array) {
    result[0] = PythonMath::tan(array[0]);
  }
};

/**
 * @brief Computes the element-wise tangent of a std::array.
 *
 * This function serves as an entry point for computing the tangents of elements
 * in a std::array. It initializes the recursive computation by calling TanCore
 * with the appropriate parameters.
 *
 * @tparam T The type of the elements in the array.
 * @tparam N The size of the array.
 * @param result The array to store the results.
 * @param array The input array whose elements will be transformed.
 */
template <typename T, std::size_t N>
inline void compute(std::array<T, N> &result, const std::array<T, N> &array) {
  TanCore<T, N, N - 1>::compute(result, array);
}

} // namespace TanAction

/**
 * @brief Applies the tangent function element-wise to a std::array.
 *
 * This function takes a std::array of type T and size N, computes the tangent
 * of each element using PythonMath::tan, and returns a new std::array
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
  TanAction::compute(result, array);
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
 * element is the result of applying PythonMath::atan to the corresponding
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
    result.push_back(PythonMath::atan(element));
  }
  return result;
}

namespace AtanAction {

template <typename T, std::size_t N, std::size_t Index> struct AtanCore {
  /**
   * @brief Recursively computes the element-wise arctangent of a std::array.
   *
   * @param result The array to store the results.
   * @param array The input array whose elements will be transformed.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array) {
    result[Index] = PythonMath::atan(array[Index]);
    AtanCore<T, N, Index - 1>::compute(result, array);
  }
};

// Specialization to end the recursion
template <typename T, std::size_t N> struct AtanCore<T, N, 0> {
  /**
   * @brief Base case for the recursive computation of arctangents in a
   * std::array.
   *
   * @param result The array to store the results.
   * @param array The input array whose elements will be transformed.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array) {
    result[0] = PythonMath::atan(array[0]);
  }
};

/**
 * @brief Computes the element-wise arctangent of a std::array.
 *
 * This function serves as an entry point for computing the arctangents of
 * elements in a std::array. It initializes the recursive computation by calling
 * AtanCore with the appropriate parameters.
 *
 * @tparam T The type of the elements in the array.
 * @tparam N The size of the array.
 * @param result The array to store the results.
 * @param array The input array whose elements will be transformed.
 */
template <typename T, std::size_t N>
inline void compute(std::array<T, N> &result, const std::array<T, N> &array) {
  AtanCore<T, N, N - 1>::compute(result, array);
}

} // namespace AtanAction

/**
 * @brief Applies the arctangent (atan) function element-wise to a std::array.
 *
 * This function takes a std::array of type T and size N, and returns a new
 * array where each element is the result of applying PythonMath::atan to the
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
  AtanAction::compute(result, array);
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
 * `PythonMath::atan2(y_i, x)` for each element `y_i` in `vector_y`.
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
    result.push_back(PythonMath::atan2(element, x));
  }
  return result;
}

/**
 * @brief Computes the element-wise arc tangent of two variables (y, x) for a
 * scalar y and a vector x.
 *
 * This function takes a scalar value `y` and a vector of values `vector_x`, and
 * computes the arc tangent of each pair (y, x_i) using `PythonMath::atan2`. The
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
    result.push_back(PythonMath::atan2(y, element));
  }
  return result;
}

/**
 * @brief Computes the element-wise arc tangent of two vectors (atan2).
 *
 * This function takes two vectors of the same size, `vector_y` and `vector_x`,
 * and computes the arc tangent of the quotient of their corresponding elements
 * using `PythonMath::atan2`. The result is a vector containing the computed
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
    result.push_back(PythonMath::atan2(vector_y[i], vector_x[i]));
  }
  return result;
}

namespace Atan2Action1 {

template <typename T, std::size_t N, std::size_t Index> struct Atan2Core1 {
  /**
   * @brief Recursively computes element-wise atan2(array_y, scalar x).
   *
   * @param result The array to store the results.
   * @param array_y The y-coordinate array.
   * @param x The scalar x-coordinate.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array_y,
                      const T &x) {
    result[Index] = PythonMath::atan2(array_y[Index], x);
    Atan2Core1<T, N, Index - 1>::compute(result, array_y, x);
  }
};

// Specialization to end the recursion
template <typename T, std::size_t N> struct Atan2Core1<T, N, 0> {
  /**
   * @brief Base case for atan2(array_y, scalar x).
   *
   * @param result The array to store the results.
   * @param array_y The y-coordinate array.
   * @param x The scalar x-coordinate.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array_y,
                      const T &x) {
    result[0] = PythonMath::atan2(array_y[0], x);
  }
};

/**
 * @brief Entry point for atan2(array_y, scalar x).
 *
 * @tparam T The type of the elements.
 * @tparam N The size of the array.
 * @param result The array to store the results.
 * @param array_y The y-coordinate array.
 * @param x The scalar x-coordinate.
 */
template <typename T, std::size_t N>
inline void compute(std::array<T, N> &result, const std::array<T, N> &array_y,
                    const T &x) {
  Atan2Core1<T, N, N - 1>::compute(result, array_y, x);
}

} // namespace Atan2Action1

namespace Atan2Action2 {

template <typename T, std::size_t N, std::size_t Index> struct Atan2Core2 {
  /**
   * @brief Recursively computes element-wise atan2(scalar y, array_x).
   *
   * @param result The array to store the results.
   * @param y The scalar y-coordinate.
   * @param array_x The x-coordinate array.
   */
  static void compute(std::array<T, N> &result, const T &y,
                      const std::array<T, N> &array_x) {
    result[Index] = PythonMath::atan2(y, array_x[Index]);
    Atan2Core2<T, N, Index - 1>::compute(result, y, array_x);
  }
};

// Specialization to end the recursion
template <typename T, std::size_t N> struct Atan2Core2<T, N, 0> {
  /**
   * @brief Base case for atan2(scalar y, array_x).
   *
   * @param result The array to store the results.
   * @param y The scalar y-coordinate.
   * @param array_x The x-coordinate array.
   */
  static void compute(std::array<T, N> &result, const T &y,
                      const std::array<T, N> &array_x) {
    result[0] = PythonMath::atan2(y, array_x[0]);
  }
};

/**
 * @brief Entry point for atan2(scalar y, array_x).
 *
 * @tparam T The type of the elements.
 * @tparam N The size of the array.
 * @param result The array to store the results.
 * @param y The scalar y-coordinate.
 * @param array_x The x-coordinate array.
 */
template <typename T, std::size_t N>
inline void compute(std::array<T, N> &result, const T &y,
                    const std::array<T, N> &array_x) {
  Atan2Core2<T, N, N - 1>::compute(result, y, array_x);
}

} // namespace Atan2Action2

namespace Atan2Action3 {

template <typename T, std::size_t N, std::size_t Index> struct Atan2Core3 {
  /**
   * @brief Recursively computes element-wise atan2(array_y, array_x).
   *
   * @param result The array to store the results.
   * @param array_y The y-coordinate array.
   * @param array_x The x-coordinate array.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array_y,
                      const std::array<T, N> &array_x) {
    result[Index] = PythonMath::atan2(array_y[Index], array_x[Index]);
    Atan2Core3<T, N, Index - 1>::compute(result, array_y, array_x);
  }
};

// Specialization to end the recursion
template <typename T, std::size_t N> struct Atan2Core3<T, N, 0> {
  /**
   * @brief Base case for atan2(array_y, array_x).
   *
   * @param result The array to store the results.
   * @param array_y The y-coordinate array.
   * @param array_x The x-coordinate array.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array_y,
                      const std::array<T, N> &array_x) {
    result[0] = PythonMath::atan2(array_y[0], array_x[0]);
  }
};

/**
 * @brief Entry point for atan2(array_y, array_x).
 *
 * @tparam T The type of the elements.
 * @tparam N The size of the arrays.
 * @param result The array to store the results.
 * @param array_y The y-coordinate array.
 * @param array_x The x-coordinate array.
 */
template <typename T, std::size_t N>
inline void compute(std::array<T, N> &result, const std::array<T, N> &array_y,
                    const std::array<T, N> &array_x) {
  Atan2Core3<T, N, N - 1>::compute(result, array_y, array_x);
}

} // namespace Atan2Action3

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
  Atan2Action1::compute(result, array_y, x);
  return result;
}

/**
 * @brief Computes the element-wise arc tangent of y and each element in the
 * input array x.
 *
 * This function takes a scalar value y and an array of values array_x, and
 * returns a new array where each element is the result of PythonMath::atan2(y,
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
  Atan2Action2::compute(result, y, array_x);
  return result;
}

/**
 * @brief Computes the element-wise arc tangent of two arrays (atan2).
 *
 * This function takes two input arrays, `array_y` and `array_x`, each of size
 * `N`, and computes the arc tangent of the quotient of their corresponding
 * elements, storing the result in a new array. The computation is performed
 * using `PythonMath::atan2` for each element.
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
  Atan2Action3::compute(result, array_y, array_x);
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
 * element is the result of applying PythonMath::asin to the corresponding
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
    result.push_back(PythonMath::asin(element));
  }
  return result;
}

namespace AsinAction {

template <typename T, std::size_t N, std::size_t Index> struct AsinCore {
  /**
   * @brief Recursively computes the element-wise arcsine of a std::array.
   *
   * @param result The array to store the results.
   * @param array The input array whose elements will be processed.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array) {
    result[Index] = PythonMath::asin(array[Index]);
    AsinCore<T, N, Index - 1>::compute(result, array);
  }
};

// Specialization to end the recursion
template <typename T, std::size_t N> struct AsinCore<T, N, 0> {
  /**
   * @brief Base case for the recursive computation of arcsines in a std::array.
   *
   * @param result The array to store the results.
   * @param array The input array whose elements will be processed.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array) {
    result[0] = PythonMath::asin(array[0]);
  }
};

/**
 * @brief Computes the element-wise arcsine of a std::array.
 *
 * This function serves as an entry point for computing the arcsines of elements
 * in a std::array. It initializes the recursive computation by calling AsinCore
 * with the appropriate parameters.
 *
 * @tparam T The type of the elements in the array.
 * @tparam N The size of the array.
 * @param result The array to store the results.
 * @param array The input array whose elements will be processed.
 */
template <typename T, std::size_t N>
inline void compute(std::array<T, N> &result, const std::array<T, N> &array) {
  AsinCore<T, N, N - 1>::compute(result, array);
}

} // namespace AsinAction

/**
 * @brief Computes the element-wise arcsine (inverse sine) of the input array.
 *
 * This function takes a std::array of type T and size N, and returns a new
 * array where each element is the result of applying PythonMath::asin to the
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
  AsinAction::compute(result, array);
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
 * (acos) to each element using PythonMath::acos, returning a new vector with
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
    result.push_back(PythonMath::acos(element));
  }
  return result;
}

namespace AcosAction {

template <typename T, std::size_t N, std::size_t Index> struct AcosCore {
  /**
   * @brief Recursively computes the element-wise arccosine of a std::array.
   *
   * @param result The array to store the results.
   * @param array The input array whose elements will be processed.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array) {
    result[Index] = PythonMath::acos(array[Index]);
    AcosCore<T, N, Index - 1>::compute(result, array);
  }
};

// Specialization to end the recursion
template <typename T, std::size_t N> struct AcosCore<T, N, 0> {
  /**
   * @brief Base case for the recursive computation of arccosines in a
   * std::array.
   *
   * @param result The array to store the results.
   * @param array The input array whose elements will be processed.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array) {
    result[0] = PythonMath::acos(array[0]);
  }
};

/**
 * @brief Computes the element-wise arccosine of a std::array.
 *
 * This function serves as an entry point for computing the arccosines of
 * elements in a std::array. It initializes the recursive computation by calling
 * AcosCore with the appropriate parameters.
 *
 * @tparam T The type of the elements in the array.
 * @tparam N The size of the array.
 * @param result The array to store the results.
 * @param array The input array whose elements will be processed.
 */
template <typename T, std::size_t N>
inline void compute(std::array<T, N> &result, const std::array<T, N> &array) {
  AcosCore<T, N, N - 1>::compute(result, array);
}

} // namespace AcosAction

/**
 * @brief Computes the arc cosine (inverse cosine) of each element in the input
 * array.
 *
 * This function applies the PythonMath::acos function to each element of the
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
  AcosAction::compute(result, array);
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
 * This function applies the PythonMath::sinh operation to every element of the
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
    result.push_back(PythonMath::sinh(element));
  }
  return result;
}

namespace SinhAction {

template <typename T, std::size_t N, std::size_t Index> struct SinhCore {
  /**
   * @brief Recursively computes the element-wise hyperbolic sine of a
   * std::array.
   *
   * @param result The array to store the results.
   * @param array The input array whose elements will be transformed.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array) {
    result[Index] = PythonMath::sinh(array[Index]);
    SinhCore<T, N, Index - 1>::compute(result, array);
  }
};

// Specialization to end the recursion
template <typename T, std::size_t N> struct SinhCore<T, N, 0> {
  /**
   * @brief Base case for the recursive computation of hyperbolic sines in a
   * std::array.
   *
   * @param result The array to store the results.
   * @param array The input array whose elements will be transformed.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array) {
    result[0] = PythonMath::sinh(array[0]);
  }
};

/**
 * @brief Computes the element-wise hyperbolic sine of a std::array.
 *
 * This function serves as an entry point for computing the hyperbolic sines of
 * elements in a std::array. It initializes the recursive computation by calling
 * SinhCore with the appropriate parameters.
 *
 * @tparam T The type of the elements in the array.
 * @tparam N The size of the array.
 * @param result The array to store the results.
 * @param array The input array whose elements will be transformed.
 */
template <typename T, std::size_t N>
inline void compute(std::array<T, N> &result, const std::array<T, N> &array) {
  SinhCore<T, N, N - 1>::compute(result, array);
}

} // namespace SinhAction

/**
 * @brief Computes the hyperbolic sine (sinh) of each element in the input
 * array.
 *
 * This function applies the PythonMath::sinh operation to each element of the
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
  SinhAction::compute(result, array);
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
 * This function applies the PythonMath::cosh function to each element of the
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
    result.push_back(PythonMath::cosh(element));
  }
  return result;
}

namespace CoshAction {

template <typename T, std::size_t N, std::size_t Index> struct CoshCore {
  /**
   * @brief Recursively computes the element-wise hyperbolic cosine of a
   * std::array.
   *
   * @param result The array to store the results.
   * @param array The input array whose elements will be transformed.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array) {
    result[Index] = PythonMath::cosh(array[Index]);
    CoshCore<T, N, Index - 1>::compute(result, array);
  }
};

// Specialization to end the recursion
template <typename T, std::size_t N> struct CoshCore<T, N, 0> {
  /**
   * @brief Base case for the recursive computation of hyperbolic cosines in a
   * std::array.
   *
   * @param result The array to store the results.
   * @param array The input array whose elements will be transformed.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array) {
    result[0] = PythonMath::cosh(array[0]);
  }
};

/**
 * @brief Computes the element-wise hyperbolic cosine of a std::array.
 *
 * This function serves as an entry point for computing the hyperbolic cosines
 * of elements in a std::array. It initializes the recursive computation by
 * calling CoshCore with the appropriate parameters.
 *
 * @tparam T The type of the elements in the array.
 * @tparam N The size of the array.
 * @param result The array to store the results.
 * @param array The input array whose elements will be transformed.
 */
template <typename T, std::size_t N>
inline void compute(std::array<T, N> &result, const std::array<T, N> &array) {
  CoshCore<T, N, N - 1>::compute(result, array);
}

} // namespace CoshAction

/**
 * @brief Computes the hyperbolic cosine (cosh) of each element in the input
 * array.
 *
 * This function applies the PythonMath::cosh function to each element of the
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
  CoshAction::compute(result, array);
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
    result.push_back(PythonMath::tanh(element));
  }
  return result;
}

namespace TanhAction {

template <typename T, std::size_t N, std::size_t Index> struct TanhCore {
  /**
   * @brief Recursively computes the element-wise hyperbolic tangent of a
   * std::array.
   *
   * @param result The array to store the results.
   * @param array The input array whose elements will be transformed.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array) {
    result[Index] = PythonMath::tanh(array[Index]);
    TanhCore<T, N, Index - 1>::compute(result, array);
  }
};

// Specialization to end the recursion
template <typename T, std::size_t N> struct TanhCore<T, N, 0> {
  /**
   * @brief Base case for the recursive computation of hyperbolic tangents in a
   * std::array.
   *
   * @param result The array to store the results.
   * @param array The input array whose elements will be transformed.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array) {
    result[0] = PythonMath::tanh(array[0]);
  }
};

/**
 * @brief Computes the element-wise hyperbolic tangent of a std::array.
 *
 * This function serves as an entry point for computing the hyperbolic tangents
 * of elements in a std::array. It initializes the recursive computation by
 * calling TanhCore with the appropriate parameters.
 *
 * @tparam T The type of the elements in the array.
 * @tparam N The size of the array.
 * @param result The array to store the results.
 * @param array The input array whose elements will be transformed.
 */
template <typename T, std::size_t N>
inline void compute(std::array<T, N> &result, const std::array<T, N> &array) {
  TanhCore<T, N, N - 1>::compute(result, array);
}

} // namespace TanhAction

/**
 * @brief Applies the hyperbolic tangent function element-wise to a std::array.
 *
 * This function takes a std::array of type T and size N, and returns a new
 * std::array where each element is the result of applying PythonMath::tanh to
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
  TanhAction::compute(result, array);
  return result;
}

} // namespace PythonMath

#endif // PYTHON_MATH_TRIGONOMETRIC_HPP_
