#ifndef __PYTHON_MATH_EXPONENTIAL_LOGARITHMIC_HPP__
#define __PYTHON_MATH_EXPONENTIAL_LOGARITHMIC_HPP__

#include "base_math.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace PythonMath {

/* sqrt */

template <typename T> inline T sqrt(const T &x) { return Base::Math::sqrt(x); }

/**
 * @brief Computes the square root of each element in the input vector.
 *
 * This function takes a vector of type T and returns a new vector where each
 * element is the result of applying Base::Math::sqrt to the corresponding
 * element of the input vector.
 *
 * @tparam T The type of the elements in the input vector.
 * @param vector The input vector whose elements' square roots are to be
 * computed.
 * @return std::vector<T> A vector containing the square roots of the input
 * elements.
 */
template <typename T> inline std::vector<T> sqrt(const std::vector<T> &vector) {
  std::vector<T> result;
  result.reserve(vector.size());
  for (const auto &element : vector) {
    result.push_back(Base::Math::sqrt(element));
  }
  return result;
}

/**
 * @brief Computes the square root of each element in the input array.
 *
 * This function takes a std::array of type T and size N, and returns a new
 * array where each element is the square root of the corresponding element in
 * the input array. The square root operation is performed using
 * Base::Math::sqrt.
 *
 * @tparam T The type of the elements in the array.
 * @tparam N The size of the array.
 * @param array The input array whose elements' square roots are to be computed.
 * @return std::array<T, N> An array containing the square roots of the input
 * elements.
 */
template <typename T, std::size_t N>
inline std::array<T, N> sqrt(const std::array<T, N> &array) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::sqrt(array[i]);
  }
  return result;
}

/* exp */

/**
 * @brief Computes the exponential of the given value.
 *
 * This function returns the result of raising Euler's number (e) to the power
 * of the input value `x`. It delegates the computation to `Base::Math::exp`.
 *
 * @tparam T Numeric type of the input value.
 * @param x The exponent to raise e to.
 * @return The exponential of `x`.
 */
template <typename T> inline T exp(const T &x) { return Base::Math::exp(x); }

/**
 * @brief Applies the exponential function to each element of the input vector.
 *
 * This function takes a vector of type T and returns a new vector where each
 * element is the result of applying Base::Math::exp to the corresponding
 * element in the input vector.
 *
 * @tparam T The type of the elements in the input vector.
 * @param vector The input vector whose elements will be exponentiated.
 * @return std::vector<T> A vector containing the exponentials of the input
 * elements.
 */
template <typename T> inline std::vector<T> exp(const std::vector<T> &vector) {
  std::vector<T> result;
  result.reserve(vector.size());
  for (const auto &element : vector) {
    result.push_back(Base::Math::exp(element));
  }
  return result;
}

/**
 * @brief Applies the exponential function to each element of the input array.
 *
 * This function takes a std::array of type T and size N, and returns a new
 * array where each element is the result of applying Base::Math::exp to the
 * corresponding element of the input array.
 *
 * @tparam T The type of the elements in the array.
 * @tparam N The size of the array.
 * @param array The input array whose elements will be exponentiated.
 * @return std::array<T, N> An array containing the exponentials of the input
 * elements.
 */
template <typename T, std::size_t N>
inline std::array<T, N> exp(const std::array<T, N> &array) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::exp(array[i]);
  }
  return result;
}

/* exp2 */

/**
 * @brief Computes the base-2 exponential of the given value.
 *
 * This function returns 2 raised to the power of \p x, i.e., 2^x.
 * It forwards the computation to Base::Math::exp2.
 *
 * @tparam T Numeric type of the input value.
 * @param x The exponent to which 2 is raised.
 * @return The value of 2 raised to the power of \p x.
 */
template <typename T> inline T exp2(const T &x) { return Base::Math::exp2(x); }

/**
 * @brief Applies the base-2 exponential function to each element of the input
 * vector.
 *
 * This function computes 2 raised to the power of each element in the input
 * vector using Base::Math::exp2, and returns a new vector containing the
 * results.
 *
 * @tparam T The type of the elements in the input vector.
 * @param vector The input vector whose elements will be exponentiated.
 * @return std::vector<T> A vector containing the base-2 exponentials of the
 * input elements.
 */
template <typename T> inline std::vector<T> exp2(const std::vector<T> &vector) {
  std::vector<T> result;
  result.reserve(vector.size());
  for (const auto &element : vector) {
    result.push_back(Base::Math::exp2(element));
  }
  return result;
}

/**
 * @brief Applies the base-2 exponential function to each element of the input
 * array.
 *
 * This function computes 2 raised to the power of each element in the input
 * array, returning a new array containing the results.
 *
 * @tparam T The type of the elements in the array.
 * @tparam N The size of the array.
 * @param array The input array whose elements will be exponentiated.
 * @return std::array<T, N> An array where each element is 2 raised to the power
 * of the corresponding input element.
 */
template <typename T, std::size_t N>
inline std::array<T, N> exp2(const std::array<T, N> &array) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::exp2(array[i]);
  }
  return result;
}

/* log */

/**
 * @brief Computes the natural logarithm (base e) of the given value.
 *
 * This function is a wrapper that calls Base::Math::log to compute the natural
 * logarithm of the input value x. It is templated to support various numeric
 * types.
 *
 * @tparam T Numeric type of the input value.
 * @param x The value for which to compute the natural logarithm.
 * @return The natural logarithm of x.
 */
template <typename T> inline T log(const T &x) { return Base::Math::log(x); }

/**
 * @brief Applies the natural logarithm function to each element of the input
 * vector.
 *
 * This function takes a vector of elements of type T and returns a new vector
 * where each element is the result of applying Base::Math::log to the
 * corresponding element in the input vector.
 *
 * @tparam T The type of the elements in the input vector.
 * @param vector The input vector whose elements will be transformed.
 * @return std::vector<T> A vector containing the natural logarithm of each
 * input element.
 */
template <typename T> inline std::vector<T> log(const std::vector<T> &vector) {
  std::vector<T> result;
  result.reserve(vector.size());
  for (const auto &element : vector) {
    result.push_back(Base::Math::log(element));
  }
  return result;
}

/**
 * @brief Applies the natural logarithm element-wise to a std::array.
 *
 * This function takes a std::array of type T and size N, and returns a new
 * std::array where each element is the result of applying Base::Math::log to
 * the corresponding element of the input array.
 *
 * @tparam T The type of the elements in the array.
 * @tparam N The size of the array.
 * @param array The input array whose elements will be transformed.
 * @return std::array<T, N> A new array with the natural logarithm applied to
 * each element.
 */
template <typename T, std::size_t N>
inline std::array<T, N> log(const std::array<T, N> &array) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::log(array[i]);
  }
  return result;
}

/* log2 */

/**
 * @brief Computes the base-2 logarithm of the given value.
 *
 * This function returns the logarithm of \p x to base 2 by delegating
 * the computation to Base::Math::log2.
 *
 * @tparam T Numeric type of the input value.
 * @param x The value whose base-2 logarithm is to be computed.
 * @return The base-2 logarithm of \p x.
 */
template <typename T> inline T log2(const T &x) { return Base::Math::log2(x); }

/**
 * @brief Computes the base-2 logarithm of each element in the input vector.
 *
 * This function takes a vector of type T and returns a new vector where each
 * element is the result of applying the base-2 logarithm (log2) to the
 * corresponding element in the input vector. The log2 operation is performed
 * using Base::Math::log2.
 *
 * @tparam T The type of the elements in the input vector.
 * @param vector The input vector containing elements to compute the log2 for.
 * @return std::vector<T> A vector containing the base-2 logarithms of the input
 * elements.
 */
template <typename T> inline std::vector<T> log2(const std::vector<T> &vector) {
  std::vector<T> result;
  result.reserve(vector.size());
  for (const auto &element : vector) {
    result.push_back(Base::Math::log2(element));
  }
  return result;
}

/**
 * @brief Computes the base-2 logarithm of each element in the input array.
 *
 * This function takes a std::array of type T and size N, and returns a new
 * array where each element is the result of applying Base::Math::log2 to the
 * corresponding element in the input array.
 *
 * @tparam T The type of the elements in the array.
 * @tparam N The size of the array.
 * @param array The input array whose elements will be transformed.
 * @return std::array<T, N> An array containing the base-2 logarithms of the
 * input elements.
 */
template <typename T, std::size_t N>
inline std::array<T, N> log2(const std::array<T, N> &array) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::log2(array[i]);
  }
  return result;
}

/* log10 */

/**
 * @brief Computes the base-10 logarithm of the given value.
 *
 * This function returns the logarithm of the input value `x` to base 10.
 * It delegates the computation to `Base::Math::log10`.
 *
 * @tparam T The type of the input value. Must support logarithmic operations.
 * @param x The value whose base-10 logarithm is to be computed.
 * @return The base-10 logarithm of `x`.
 */
template <typename T> inline T log10(const T &x) {
  return Base::Math::log10(x);
}

/**
 * @brief Computes the base-10 logarithm of each element in the input vector.
 *
 * This function takes a vector of type T and returns a new vector where each
 * element is the result of applying Base::Math::log10 to the corresponding
 * element of the input vector.
 *
 * @tparam T The type of the elements in the input vector.
 * @param vector The input vector whose elements will be transformed.
 * @return std::vector<T> A vector containing the base-10 logarithms of the
 * input elements.
 */
template <typename T>
inline std::vector<T> log10(const std::vector<T> &vector) {
  std::vector<T> result;
  result.reserve(vector.size());
  for (const auto &element : vector) {
    result.push_back(Base::Math::log10(element));
  }
  return result;
}

/**
 * @brief Computes the base-10 logarithm of each element in the input array.
 *
 * This function takes a std::array of type T and size N, and returns a new
 * array where each element is the result of applying Base::Math::log10 to the
 * corresponding element of the input array.
 *
 * @tparam T The type of the elements in the array.
 * @tparam N The size of the array.
 * @param array The input array whose elements will be transformed.
 * @return std::array<T, N> An array containing the base-10 logarithm of each
 * input element.
 */
template <typename T, std::size_t N>
inline std::array<T, N> log10(const std::array<T, N> &array) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::log10(array[i]);
  }
  return result;
}

/* pow */

/**
 * @brief Computes the value of x raised to the power of y.
 *
 * This function is a wrapper around Base::Math::pow, providing
 * a generic interface for exponentiation with any type T that
 * supports the operation.
 *
 * @tparam T The type of the base and exponent.
 * @param x The base value.
 * @param y The exponent value.
 * @return The result of raising x to the power of y.
 */
template <typename T> inline T pow(const T &x, const T &y) {
  return Base::Math::pow(x, y);
}

/**
 * @brief Raises each element of the input vector to the given power.
 *
 * This function takes a vector of elements and a scalar exponent, and returns a
 * new vector where each element is the result of raising the corresponding
 * input element to the power of y. The exponentiation is performed using
 * Base::Math::pow.
 *
 * @tparam T The type of the elements in the vector.
 * @param vector_x The input vector whose elements will be exponentiated.
 * @param y The exponent to which each element of the vector will be raised.
 * @return std::vector<T> A vector containing the results of the exponentiation.
 */
template <typename T>
inline std::vector<T> pow(const std::vector<T> &vector_x, const T &y) {
  std::vector<T> result;
  result.reserve(vector_x.size());
  for (const auto &element : vector_x) {
    result.push_back(Base::Math::pow(element, y));
  }
  return result;
}

/**
 * @brief Raises a scalar value to the power of each element in a vector.
 *
 * This function computes the result of raising the scalar value `x` to the
 * power of each element in the input vector `vector_y`. The results are stored
 * in a new vector, where each element is calculated as `Base::Math::pow(x,
 * y_i)` for each `y_i` in `vector_y`.
 *
 * @tparam T The numeric type of the scalar and vector elements.
 * @param x The base scalar value to be raised to the power of each element in
 * `vector_y`.
 * @param vector_y A vector containing the exponents.
 * @return std::vector<T> A vector containing the results of the exponentiation.
 */
template <typename T>
inline std::vector<T> pow(const T &x, const std::vector<T> &vector_y) {
  std::vector<T> result;
  result.reserve(vector_y.size());
  for (const auto &element : vector_y) {
    result.push_back(Base::Math::pow(x, element));
  }
  return result;
}

/**
 * @brief Computes the element-wise power of two vectors.
 *
 * This function takes two vectors of the same size, `vector_x` and `vector_y`,
 * and returns a new vector where each element is the result of raising the
 * corresponding element in `vector_x` to the power of the corresponding element
 * in `vector_y`, using `Base::Math::pow`.
 *
 * @tparam T The type of the elements in the input vectors.
 * @param vector_x The base values as a vector.
 * @param vector_y The exponent values as a vector.
 * @return std::vector<T> A vector containing the element-wise powers.
 *
 * @note The input vectors must have the same size.
 */
template <typename T>
inline std::vector<T> pow(const std::vector<T> &vector_x,
                          const std::vector<T> &vector_y) {

  std::vector<T> result;
  result.reserve(vector_x.size());
  for (std::size_t i = 0; i < vector_x.size(); ++i) {
    result.push_back(Base::Math::pow(vector_x[i], vector_y[i]));
  }
  return result;
}

/**
 * @brief Raises each element of the input array to the given power.
 *
 * This function takes an input array of type T and size N, and returns a new
 * array where each element is raised to the power of y using Base::Math::pow.
 *
 * @tparam T The type of the elements in the array.
 * @tparam N The size of the array.
 * @param array_x The input array whose elements will be exponentiated.
 * @param y The exponent to which each element of the array will be raised.
 * @return std::array<T, N> A new array with each element as array_x[i] raised
 * to the power y.
 */
template <typename T, std::size_t N>
inline std::array<T, N> pow(const std::array<T, N> &array_x, const T &y) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::pow(array_x[i], y);
  }
  return result;
}

/**
 * @brief Raises a scalar value to the power of each element in an array.
 *
 * This function computes the power of a scalar base `x` raised to each exponent
 * in the input array `array_y`. The result is an array where each element is
 * calculated as `pow(x, array_y[i])`.
 *
 * @tparam T The type of the scalar and array elements.
 * @tparam N The size of the input and output arrays.
 * @param x The scalar base value.
 * @param array_y The array of exponents.
 * @return std::array<T, N> An array containing the results of raising `x` to
 * each exponent in `array_y`.
 */
template <typename T, std::size_t N>
inline std::array<T, N> pow(const T &x, const std::array<T, N> &array_y) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::pow(x, array_y[i]);
  }
  return result;
}

/**
 * @brief Computes the element-wise power of two arrays.
 *
 * This function takes two arrays of the same size and computes the power of
 * each corresponding element, i.e., result[i] = pow(array_x[i], array_y[i]),
 * using Base::Math::pow.
 *
 * @tparam T The type of the elements in the arrays.
 * @tparam N The size of the arrays.
 * @param array_x The base values as an array.
 * @param array_y The exponent values as an array.
 * @return std::array<T, N> An array containing the result of raising each
 * element of array_x to the power of the corresponding element in array_y.
 */
template <typename T, std::size_t N>
inline std::array<T, N> pow(const std::array<T, N> &array_x,
                            const std::array<T, N> &array_y) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::pow(array_x[i], array_y[i]);
  }
  return result;
}

} // namespace PythonMath

#endif // __PYTHON_MATH_EXPONENTIAL_LOGARITHMIC_HPP__
