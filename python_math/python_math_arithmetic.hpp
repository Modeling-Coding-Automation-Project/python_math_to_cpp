#ifndef PYTHON_MATH_ARITHMETIC_HPP
#define PYTHON_MATH_ARITHMETIC_HPP

#include "base_math.hpp"

#include <array>
#include <vector>

namespace PythonMath {

/* abs */
template <typename T> inline T abs(const T &x) { return Base::Math::abs(x); }

template <typename T> inline std::vector<T> abs(const std::vector<T> &vector) {
  std::vector<T> result;
  result.reserve(vector.size());
  for (const auto &element : vector) {
    result.push_back(Base::Math::abs(element));
  }
  return result;
}

template <typename T, std::size_t N>
inline std::array<T, N> abs(const std::array<T, N> &array) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::abs(array[i]);
  }
  return result;
}

/* fmod */
template <typename T> inline T fmod(const T &x, const T &y) {
  return Base::Math::mod(x, y);
}

template <typename T>
inline std::vector<T> fmod(const std::vector<T> &vector, const T &y) {
  std::vector<T> result;
  result.reserve(vector.size());
  for (const auto &element : vector) {
    result.push_back(Base::Math::mod(element, y));
  }
  return result;
}

template <typename T, std::size_t N>
inline std::array<T, N> fmod(const std::array<T, N> &array, const T &y) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::mod(array[i], y);
  }
  return result;
}

} // namespace PythonMath

#endif // PYTHON_MATH_ARITHMETIC_HPP
