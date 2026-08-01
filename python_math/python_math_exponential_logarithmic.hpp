#ifndef PYTHON_MATH_EXPONENTIAL_LOGARITHMIC_HPP_
#define PYTHON_MATH_EXPONENTIAL_LOGARITHMIC_HPP_

#include "base_math.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace PythonMath {

/* sqrt */

template <typename T, std::size_t N>
inline std::array<T, N> sqrt(const std::array<T, N> &array);

template <typename T> inline T sqrt(const T &x) { return Base::Math::sqrt(x); }

/**
 * @brief Computes the square root of each element in the input vector.
 *
 * This function takes a vector of type T and returns a new vector where each
 * element is the result of applying PythonMath::sqrt to the corresponding
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
    result.push_back(PythonMath::sqrt(element));
  }
  return result;
}

namespace SqrtAction {

template <typename T, std::size_t N, std::size_t Index> struct SqrtCore {
  /**
   * @brief Recursively computes the element-wise square root of a std::array.
   *
   * @param result The array to store the results.
   * @param array The input array whose elements' square roots are to be
   * computed.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array) {
    result[Index] = PythonMath::sqrt(array[Index]);
    SqrtCore<T, N, Index - 1>::compute(result, array);
  }
};

// Specialization to end the recursion
template <typename T, std::size_t N> struct SqrtCore<T, N, 0> {
  /**
   * @brief Base case for the recursive computation of square roots in a
   * std::array.
   *
   * @param result The array to store the results.
   * @param array The input array whose elements' square roots are to be
   * computed.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array) {
    result[0] = PythonMath::sqrt(array[0]);
  }
};

/**
 * @brief Computes the element-wise square root of a std::array.
 *
 * This function serves as an entry point for computing the square roots of
 * elements in a std::array. It initializes the recursive computation by calling
 * SqrtCore with the appropriate parameters.
 *
 * @tparam T The type of the elements in the array.
 * @tparam N The size of the array.
 * @param result The array to store the results.
 * @param array The input array whose elements' square roots are to be computed.
 */
template <typename T, std::size_t N>
inline void compute(std::array<T, N> &result, const std::array<T, N> &array) {
  SqrtCore<T, N, N - 1>::compute(result, array);
}

} // namespace SqrtAction

/**
 * @brief Computes the square root of each element in the input array.
 *
 * This function takes a std::array of type T and size N, and returns a new
 * array where each element is the square root of the corresponding element in
 * the input array. The square root operation is performed using
 * PythonMath::sqrt.
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
  SqrtAction::compute(result, array);
  return result;
}

/* exp */

template <typename T, std::size_t N>
inline std::array<T, N> exp(const std::array<T, N> &array);

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
 * element is the result of applying PythonMath::exp to the corresponding
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
    result.push_back(PythonMath::exp(element));
  }
  return result;
}

namespace ExpAction {

template <typename T, std::size_t N, std::size_t Index> struct ExpCore {
  /**
   * @brief Recursively computes the element-wise exponential of a std::array.
   *
   * @param result The array to store the results.
   * @param array The input array whose elements will be exponentiated.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array) {
    result[Index] = PythonMath::exp(array[Index]);
    ExpCore<T, N, Index - 1>::compute(result, array);
  }
};

// Specialization to end the recursion
template <typename T, std::size_t N> struct ExpCore<T, N, 0> {
  /**
   * @brief Base case for the recursive computation of exponentials in a
   * std::array.
   *
   * @param result The array to store the results.
   * @param array The input array whose elements will be exponentiated.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array) {
    result[0] = PythonMath::exp(array[0]);
  }
};

/**
 * @brief Computes the element-wise exponential of a std::array.
 *
 * This function serves as an entry point for computing the exponentials of
 * elements in a std::array. It initializes the recursive computation by calling
 * ExpCore with the appropriate parameters.
 *
 * @tparam T The type of the elements in the array.
 * @tparam N The size of the array.
 * @param result The array to store the results.
 * @param array The input array whose elements will be exponentiated.
 */
template <typename T, std::size_t N>
inline void compute(std::array<T, N> &result, const std::array<T, N> &array) {
  ExpCore<T, N, N - 1>::compute(result, array);
}

} // namespace ExpAction

/**
 * @brief Applies the exponential function to each element of the input array.
 *
 * This function takes a std::array of type T and size N, and returns a new
 * array where each element is the result of applying PythonMath::exp to the
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
  ExpAction::compute(result, array);
  return result;
}

/* exp2 */

template <typename T, std::size_t N>
inline std::array<T, N> exp2(const std::array<T, N> &array);

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
 * vector using PythonMath::exp2, and returns a new vector containing the
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
    result.push_back(PythonMath::exp2(element));
  }
  return result;
}

namespace Exp2Action {

template <typename T, std::size_t N, std::size_t Index> struct Exp2Core {
  /**
   * @brief Recursively computes the element-wise base-2 exponential of a
   * std::array.
   *
   * @param result The array to store the results.
   * @param array The input array whose elements will be exponentiated.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array) {
    result[Index] = PythonMath::exp2(array[Index]);
    Exp2Core<T, N, Index - 1>::compute(result, array);
  }
};

// Specialization to end the recursion
template <typename T, std::size_t N> struct Exp2Core<T, N, 0> {
  /**
   * @brief Base case for the recursive computation of base-2 exponentials in a
   * std::array.
   *
   * @param result The array to store the results.
   * @param array The input array whose elements will be exponentiated.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array) {
    result[0] = PythonMath::exp2(array[0]);
  }
};

/**
 * @brief Computes the element-wise base-2 exponential of a std::array.
 *
 * This function serves as an entry point for computing the base-2 exponentials
 * of elements in a std::array. It initializes the recursive computation by
 * calling Exp2Core with the appropriate parameters.
 *
 * @tparam T The type of the elements in the array.
 * @tparam N The size of the array.
 * @param result The array to store the results.
 * @param array The input array whose elements will be exponentiated.
 */
template <typename T, std::size_t N>
inline void compute(std::array<T, N> &result, const std::array<T, N> &array) {
  Exp2Core<T, N, N - 1>::compute(result, array);
}

} // namespace Exp2Action

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
  Exp2Action::compute(result, array);
  return result;
}

/* log */

template <typename T, std::size_t N>
inline void compute(std::array<T, N> &result, const std::array<T, N> &array);

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
 * where each element is the result of applying PythonMath::log to the
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
    result.push_back(PythonMath::log(element));
  }
  return result;
}

namespace LogAction {

template <typename T, std::size_t N, std::size_t Index> struct LogCore {
  /**
   * @brief Recursively computes the element-wise natural logarithm of a
   * std::array.
   *
   * @param result The array to store the results.
   * @param array The input array whose elements will be transformed.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array) {
    result[Index] = PythonMath::log(array[Index]);
    LogCore<T, N, Index - 1>::compute(result, array);
  }
};

// Specialization to end the recursion
template <typename T, std::size_t N> struct LogCore<T, N, 0> {
  /**
   * @brief Base case for the recursive computation of natural logarithms in a
   * std::array.
   *
   * @param result The array to store the results.
   * @param array The input array whose elements will be transformed.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array) {
    result[0] = PythonMath::log(array[0]);
  }
};

/**
 * @brief Computes the element-wise natural logarithm of a std::array.
 *
 * This function serves as an entry point for computing the natural logarithms
 * of elements in a std::array. It initializes the recursive computation by
 * calling LogCore with the appropriate parameters.
 *
 * @tparam T The type of the elements in the array.
 * @tparam N The size of the array.
 * @param result The array to store the results.
 * @param array The input array whose elements will be transformed.
 */
template <typename T, std::size_t N>
inline void compute(std::array<T, N> &result, const std::array<T, N> &array) {
  LogCore<T, N, N - 1>::compute(result, array);
}

} // namespace LogAction

/**
 * @brief Applies the natural logarithm element-wise to a std::array.
 *
 * This function takes a std::array of type T and size N, and returns a new
 * std::array where each element is the result of applying PythonMath::log to
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
  LogAction::compute(result, array);
  return result;
}

/* log2 */

template <typename T, std::size_t N>
inline std::array<T, N> log2(const std::array<T, N> &array);

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
 * using PythonMath::log2.
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
    result.push_back(PythonMath::log2(element));
  }
  return result;
}

namespace Log2Action {

template <typename T, std::size_t N, std::size_t Index> struct Log2Core {
  /**
   * @brief Recursively computes the element-wise base-2 logarithm of a
   * std::array.
   *
   * @param result The array to store the results.
   * @param array The input array whose elements will be transformed.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array) {
    result[Index] = PythonMath::log2(array[Index]);
    Log2Core<T, N, Index - 1>::compute(result, array);
  }
};

// Specialization to end the recursion
template <typename T, std::size_t N> struct Log2Core<T, N, 0> {
  /**
   * @brief Base case for the recursive computation of base-2 logarithms in a
   * std::array.
   *
   * @param result The array to store the results.
   * @param array The input array whose elements will be transformed.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array) {
    result[0] = PythonMath::log2(array[0]);
  }
};

/**
 * @brief Computes the element-wise base-2 logarithm of a std::array.
 *
 * This function serves as an entry point for computing the base-2 logarithms
 * of elements in a std::array. It initializes the recursive computation by
 * calling Log2Core with the appropriate parameters.
 *
 * @tparam T The type of the elements in the array.
 * @tparam N The size of the array.
 * @param result The array to store the results.
 * @param array The input array whose elements will be transformed.
 */
template <typename T, std::size_t N>
inline void compute(std::array<T, N> &result, const std::array<T, N> &array) {
  Log2Core<T, N, N - 1>::compute(result, array);
}

} // namespace Log2Action

/**
 * @brief Computes the base-2 logarithm of each element in the input array.
 *
 * This function takes a std::array of type T and size N, and returns a new
 * array where each element is the result of applying PythonMath::log2 to the
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
  Log2Action::compute(result, array);
  return result;
}

/* log10 */

template <typename T, std::size_t N>
inline std::array<T, N> log10(const std::array<T, N> &array);

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
 * element is the result of applying PythonMath::log10 to the corresponding
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
    result.push_back(PythonMath::log10(element));
  }
  return result;
}

namespace Log10Action {

template <typename T, std::size_t N, std::size_t Index> struct Log10Core {
  /**
   * @brief Recursively computes the element-wise base-10 logarithm of a
   * std::array.
   *
   * @param result The array to store the results.
   * @param array The input array whose elements will be transformed.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array) {
    result[Index] = PythonMath::log10(array[Index]);
    Log10Core<T, N, Index - 1>::compute(result, array);
  }
};

// Specialization to end the recursion
template <typename T, std::size_t N> struct Log10Core<T, N, 0> {
  /**
   * @brief Base case for the recursive computation of base-10 logarithms in a
   * std::array.
   *
   * @param result The array to store the results.
   * @param array The input array whose elements will be transformed.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array) {
    result[0] = PythonMath::log10(array[0]);
  }
};

/**
 * @brief Computes the element-wise base-10 logarithm of a std::array.
 *
 * This function serves as an entry point for computing the base-10 logarithms
 * of elements in a std::array. It initializes the recursive computation by
 * calling Log10Core with the appropriate parameters.
 *
 * @tparam T The type of the elements in the array.
 * @tparam N The size of the array.
 * @param result The array to store the results.
 * @param array The input array whose elements will be transformed.
 */
template <typename T, std::size_t N>
inline void compute(std::array<T, N> &result, const std::array<T, N> &array) {
  Log10Core<T, N, N - 1>::compute(result, array);
}

} // namespace Log10Action

/**
 * @brief Computes the base-10 logarithm of each element in the input array.
 *
 * This function takes a std::array of type T and size N, and returns a new
 * array where each element is the result of applying PythonMath::log10 to the
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
  Log10Action::compute(result, array);
  return result;
}

/* pow */

template <typename T, std::size_t N>
inline std::array<T, N> pow(const std::array<T, N> &array_x, const T &y);

template <typename T, std::size_t N>
inline std::array<T, N> pow(const T &x, const std::array<T, N> &array_y);

template <typename T, std::size_t N>
inline std::array<T, N> pow(const std::array<T, N> &array_x,
                            const std::array<T, N> &array_y);

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
 * PythonMath::pow.
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
    result.push_back(PythonMath::pow(element, y));
  }
  return result;
}

/**
 * @brief Raises a scalar value to the power of each element in a vector.
 *
 * This function computes the result of raising the scalar value `x` to the
 * power of each element in the input vector `vector_y`. The results are stored
 * in a new vector, where each element is calculated as `PythonMath::pow(x,
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
    result.push_back(PythonMath::pow(x, element));
  }
  return result;
}

/**
 * @brief Computes the element-wise power of two vectors.
 *
 * This function takes two vectors of the same size, `vector_x` and `vector_y`,
 * and returns a new vector where each element is the result of raising the
 * corresponding element in `vector_x` to the power of the corresponding element
 * in `vector_y`, using `PythonMath::pow`.
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
    result.push_back(PythonMath::pow(vector_x[i], vector_y[i]));
  }
  return result;
}

namespace PowAction1 {

template <typename T, std::size_t N, std::size_t Index> struct PowCore1 {
  /**
   * @brief Recursively computes the element-wise power (array base, scalar
   * exponent).
   *
   * @param result The array to store the results.
   * @param array_x The base values as an array.
   * @param y The scalar exponent.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array_x,
                      const T &y) {
    result[Index] = PythonMath::pow(array_x[Index], y);
    PowCore1<T, N, Index - 1>::compute(result, array_x, y);
  }
};

// Specialization to end the recursion
template <typename T, std::size_t N> struct PowCore1<T, N, 0> {
  /**
   * @brief Base case for the recursive computation of power (array base, scalar
   * exponent).
   *
   * @param result The array to store the results.
   * @param array_x The base values as an array.
   * @param y The scalar exponent.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array_x,
                      const T &y) {
    result[0] = PythonMath::pow(array_x[0], y);
  }
};

/**
 * @brief Entry point for computing element-wise power (array base, scalar
 * exponent).
 *
 * @tparam T The type of the elements in the array.
 * @tparam N The size of the array.
 * @param result The array to store the results.
 * @param array_x The base values as an array.
 * @param y The scalar exponent.
 */
template <typename T, std::size_t N>
inline void compute(std::array<T, N> &result, const std::array<T, N> &array_x,
                    const T &y) {
  PowCore1<T, N, N - 1>::compute(result, array_x, y);
}

} // namespace PowAction1

namespace PowAction2 {

template <typename T, std::size_t N, std::size_t Index> struct PowCore2 {
  /**
   * @brief Recursively computes the element-wise power (scalar base, array
   * exponent).
   *
   * @param result The array to store the results.
   * @param x The scalar base.
   * @param array_y The exponent values as an array.
   */
  static void compute(std::array<T, N> &result, const T &x,
                      const std::array<T, N> &array_y) {
    result[Index] = PythonMath::pow(x, array_y[Index]);
    PowCore2<T, N, Index - 1>::compute(result, x, array_y);
  }
};

// Specialization to end the recursion
template <typename T, std::size_t N> struct PowCore2<T, N, 0> {
  /**
   * @brief Base case for the recursive computation of power (scalar base, array
   * exponent).
   *
   * @param result The array to store the results.
   * @param x The scalar base.
   * @param array_y The exponent values as an array.
   */
  static void compute(std::array<T, N> &result, const T &x,
                      const std::array<T, N> &array_y) {
    result[0] = PythonMath::pow(x, array_y[0]);
  }
};

/**
 * @brief Entry point for computing element-wise power (scalar base, array
 * exponent).
 *
 * @tparam T The type of the elements in the array.
 * @tparam N The size of the array.
 * @param result The array to store the results.
 * @param x The scalar base.
 * @param array_y The exponent values as an array.
 */
template <typename T, std::size_t N>
inline void compute(std::array<T, N> &result, const T &x,
                    const std::array<T, N> &array_y) {
  PowCore2<T, N, N - 1>::compute(result, x, array_y);
}

} // namespace PowAction2

namespace PowAction3 {

template <typename T, std::size_t N, std::size_t Index> struct PowCore3 {
  /**
   * @brief Recursively computes the element-wise power (array base, array
   * exponent).
   *
   * @param result The array to store the results.
   * @param array_x The base values as an array.
   * @param array_y The exponent values as an array.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array_x,
                      const std::array<T, N> &array_y) {
    result[Index] = PythonMath::pow(array_x[Index], array_y[Index]);
    PowCore3<T, N, Index - 1>::compute(result, array_x, array_y);
  }
};

// Specialization to end the recursion
template <typename T, std::size_t N> struct PowCore3<T, N, 0> {
  /**
   * @brief Base case for the recursive computation of power (array base, array
   * exponent).
   *
   * @param result The array to store the results.
   * @param array_x The base values as an array.
   * @param array_y The exponent values as an array.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array_x,
                      const std::array<T, N> &array_y) {
    result[0] = PythonMath::pow(array_x[0], array_y[0]);
  }
};

/**
 * @brief Entry point for computing element-wise power (array base, array
 * exponent).
 *
 * @tparam T The type of the elements in the arrays.
 * @tparam N The size of the arrays.
 * @param result The array to store the results.
 * @param array_x The base values as an array.
 * @param array_y The exponent values as an array.
 */
template <typename T, std::size_t N>
inline void compute(std::array<T, N> &result, const std::array<T, N> &array_x,
                    const std::array<T, N> &array_y) {
  PowCore3<T, N, N - 1>::compute(result, array_x, array_y);
}

} // namespace PowAction3

/**
 * @brief Raises each element of the input array to the given power.
 *
 * This function takes an input array of type T and size N, and returns a new
 * array where each element is raised to the power of y using PythonMath::pow.
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
  PowAction1::compute(result, array_x, y);
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
  PowAction2::compute(result, x, array_y);
  return result;
}

/**
 * @brief Computes the element-wise power of two arrays.
 *
 * This function takes two arrays of the same size and computes the power of
 * each corresponding element, i.e., result[i] = pow(array_x[i], array_y[i]),
 * using PythonMath::pow.
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
  PowAction3::compute(result, array_x, array_y);
  return result;
}

} // namespace PythonMath

#endif // PYTHON_MATH_EXPONENTIAL_LOGARITHMIC_HPP_
