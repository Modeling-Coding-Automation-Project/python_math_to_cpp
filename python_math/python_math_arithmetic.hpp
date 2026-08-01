/**
 * @file python_math_arithmetic.hpp
 * @brief Provides Python-like arithmetic functions for C++ containers.
 *
 * This header defines the `PythonMath` namespace, which offers template
 * functions that mimic Python's arithmetic operations such as `abs` and `fmod`
 * (floating-point modulo) for scalar types as well as for standard C++
 * containers like `std::vector` and `std::array`. These functions internally
 * delegate to the corresponding implementations in the `Base::Math` namespace.
 *
 * @namespace PythonMath
 * @brief Contains template functions for element-wise arithmetic operations on
 * scalars and containers.
 */
#ifndef PYTHON_MATH_ARITHMETIC_HPP_
#define PYTHON_MATH_ARITHMETIC_HPP_

#include "base_math.hpp"

#include <array>
#include <vector>

namespace PythonMath {

/* abs */

template <typename T, std::size_t N>
inline std::array<T, N> abs(const std::array<T, N> &array);

template <typename T> inline std::vector<T> abs(const std::vector<T> &vector);

/**
 * @brief Returns the absolute value of the given input.
 *
 * This function is a template that computes the absolute value of its argument
 * by delegating to Base::Math::abs. It works for any type T for which
 * Base::Math::abs is defined.
 *
 * @tparam T The type of the input value.
 * @param x The value whose absolute value is to be computed.
 * @return The absolute value of x.
 */
template <typename T> inline T abs(const T &x) { return Base::Math::abs(x); }

/**
 * @brief Computes the element-wise absolute value of a vector.
 *
 * This function takes a constant reference to a std::vector of type T and
 * returns a new std::vector containing the absolute values of each element. The
 * function assumes that PythonMath::abs is defined for type T.
 *
 * @tparam T The type of the elements in the input vector.
 * @param vector The input vector whose elements' absolute values are to be
 * computed.
 * @return std::vector<T> A vector containing the absolute values of the input
 * vector's elements.
 */
template <typename T> inline std::vector<T> abs(const std::vector<T> &vector) {
  std::vector<T> result;
  result.reserve(vector.size());
  for (const auto &element : vector) {
    result.push_back(PythonMath::abs(element));
  }
  return result;
}

namespace AbsAction {

template <typename T, std::size_t N, std::size_t Index> struct AbsCore {
  /**
   * @brief Recursively computes the element-wise absolute value of a
   * std::array.
   *
   * @param result The array to store the results.
   * @param array The input array whose elements' absolute values are to be
   * computed.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array) {
    result[Index] = PythonMath::abs(array[Index]);
    AbsCore<T, N, Index - 1>::compute(result, array);
  }
};

// Specialization to end the recursion
template <typename T, std::size_t N> struct AbsCore<T, N, 0> {
  /**
   * @brief Base case for the recursive computation of absolute values in a
   * std::array.
   *
   * @param result The array to store the results.
   * @param array The input array whose elements' absolute values are to be
   * computed.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array) {
    result[0] = PythonMath::abs(array[0]);
  }
};

/**
 * @brief Computes the element-wise absolute value of a std::array.
 *
 * This function serves as an entry point for computing the absolute values of
 * elements in a std::array. It initializes the recursive computation by calling
 * AbsCore with the appropriate parameters.
 *
 * @tparam T The type of the elements in the array.
 * @tparam N The size of the array.
 * @param result The array to store the results.
 * @param array The input array whose elements' absolute values are to be
 * computed.
 */
template <typename T, std::size_t N>
inline void compute(std::array<T, N> &result, const std::array<T, N> &array) {
  AbsCore<T, N, N - 1>::compute(result, array);
}

} // namespace AbsAction

/**
 * @brief Computes the element-wise absolute value of a std::array.
 *
 * This function takes a constant reference to a std::array of type T and size
 * N, and returns a new std::array where each element is the absolute value of
 * the corresponding element in the input array. The absolute value is computed
 * using PythonMath::abs.
 *
 * @tparam T The type of the elements in the array.
 * @tparam N The size of the array.
 * @param array The input array whose elements' absolute values are to be
 * computed.
 * @return std::array<T, N> A new array containing the absolute values of the
 * input array's elements.
 */
template <typename T, std::size_t N>
inline std::array<T, N> abs(const std::array<T, N> &array) {
  std::array<T, N> result;

  AbsAction::compute(result, array);

  return result;
}

namespace FmodAction {

template <typename T, std::size_t N, std::size_t Index> struct FmodCore {
  /**
   * @brief Recursively computes the element-wise floating-point modulus of a
   * std::array.
   *
   * @param result The array to store the results.
   * @param array The input array whose elements will be processed.
   * @param y The divisor to use for the modulus operation.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array,
                      const T &y) {
    result[Index] = Base::Math::mod(array[Index], y);
    FmodCore<T, N, Index - 1>::compute(result, array, y);
  }
};

// Specialization to end the recursion
template <typename T, std::size_t N> struct FmodCore<T, N, 0> {
  /**
   * @brief Base case for the recursive computation of floating-point modulus
   * in a std::array.
   *
   * @param result The array to store the results.
   * @param array The input array whose elements will be processed.
   * @param y The divisor to use for the modulus operation.
   */
  static void compute(std::array<T, N> &result, const std::array<T, N> &array,
                      const T &y) {
    result[0] = Base::Math::mod(array[0], y);
  }
};

/**
 * @brief Computes the element-wise floating-point modulus of a std::array.
 *
 * This function serves as an entry point for computing the floating-point
 * modulus of elements in a std::array. It initializes the recursive computation
 * by calling FmodCore with the appropriate parameters.
 *
 * @tparam T The type of the elements in the array.
 * @tparam N The size of the array.
 * @param result The array to store the results.
 * @param array The input array whose elements will be processed.
 * @param y The divisor to use for the modulus operation.
 */
template <typename T, std::size_t N>
inline void compute(std::array<T, N> &result, const std::array<T, N> &array,
                    const T &y) {
  FmodCore<T, N, N - 1>::compute(result, array, y);
}

} // namespace FmodAction

/* fmod */

template <typename T, std::size_t N>
inline std::array<T, N> fmod(const std::array<T, N> &array, const T &y);

template <typename T>
inline std::vector<T> fmod(const std::vector<T> &vector, const T &y);

/**
 * @brief Computes the floating-point remainder of the division of two values.
 *
 * This function returns the result of the modulus operation (remainder after
 * division) for the given values x and y, using Base::Math::mod. It is a
 * generic template that works with any type T supported by the underlying mod
 * implementation.
 *
 * @tparam T The type of the input values (e.g., float, double, etc.).
 * @param x The dividend.
 * @param y The divisor.
 * @return The remainder of x divided by y.
 */
template <typename T> inline T fmod(const T &x, const T &y) {
  return Base::Math::mod(x, y);
}

/**
 * @brief Computes the element-wise floating-point modulus of a vector.
 *
 * This function applies the modulus operation to each element of the input
 * vector with respect to the given scalar value `y`, using `Base::Math::mod`.
 *
 * @tparam T The type of the elements in the vector.
 * @param vector The input vector whose elements will be processed.
 * @param y The scalar value to use as the modulus divisor.
 * @return std::vector<T> A vector containing the result of the modulus
 * operation for each element.
 */
template <typename T>
inline std::vector<T> fmod(const std::vector<T> &vector, const T &y) {
  std::vector<T> result;
  result.reserve(vector.size());
  for (const auto &element : vector) {
    result.push_back(Base::Math::mod(element, y));
  }
  return result;
}

/**
 * @brief Computes the element-wise floating-point modulus of an array.
 *
 * This function applies the modulus operation to each element of the input
 * array using the provided divisor `y`, and returns a new array containing the
 * results. The modulus operation is performed using `Base::Math::mod`.
 *
 * @tparam T The type of the elements in the array.
 * @tparam N The size of the array.
 * @param array The input array whose elements will be processed.
 * @param y The divisor to use for the modulus operation.
 * @return std::array<T, N> An array containing the result of the modulus
 * operation for each element.
 */
template <typename T, std::size_t N>
inline std::array<T, N> fmod(const std::array<T, N> &array, const T &y) {
  std::array<T, N> result;

  FmodAction::compute(result, array, y);

  return result;
}

} // namespace PythonMath

#endif // PYTHON_MATH_ARITHMETIC_HPP_
