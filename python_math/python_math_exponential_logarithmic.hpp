#ifndef PYTHON_MATH_EXPONENTIAL_LOGARITHMIC_HPP
#define PYTHON_MATH_EXPONENTIAL_LOGARITHMIC_HPP

#include "base_math.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace PythonMath {

/* sqrt */
template <typename T> inline T sqrt(const T &x) { return Base::Math::sqrt(x); }

template <typename T> inline std::vector<T> sqrt(const std::vector<T> &vector) {
  std::vector<T> result;
  result.reserve(vector.size());
  for (const auto &element : vector) {
    result.push_back(Base::Math::sqrt(element));
  }
  return result;
}

template <typename T, std::size_t N>
inline std::array<T, N> sqrt(const std::array<T, N> &array) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::sqrt(array[i]);
  }
  return result;
}

/* exp */
template <typename T> inline T exp(const T &x) { return Base::Math::exp(x); }

template <typename T> inline std::vector<T> exp(const std::vector<T> &vector) {
  std::vector<T> result;
  result.reserve(vector.size());
  for (const auto &element : vector) {
    result.push_back(Base::Math::exp(element));
  }
  return result;
}

template <typename T, std::size_t N>
inline std::array<T, N> exp(const std::array<T, N> &array) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::exp(array[i]);
  }
  return result;
}

/* exp2 */
template <typename T> inline T exp2(const T &x) { return Base::Math::exp2(x); }

template <typename T> inline std::vector<T> exp2(const std::vector<T> &vector) {
  std::vector<T> result;
  result.reserve(vector.size());
  for (const auto &element : vector) {
    result.push_back(Base::Math::exp2(element));
  }
  return result;
}

template <typename T, std::size_t N>
inline std::array<T, N> exp2(const std::array<T, N> &array) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::exp2(array[i]);
  }
  return result;
}

/* log */
template <typename T> inline T log(const T &x) { return Base::Math::log(x); }

template <typename T> inline std::vector<T> log(const std::vector<T> &vector) {
  std::vector<T> result;
  result.reserve(vector.size());
  for (const auto &element : vector) {
    result.push_back(Base::Math::log(element));
  }
  return result;
}

template <typename T, std::size_t N>
inline std::array<T, N> log(const std::array<T, N> &array) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::log(array[i]);
  }
  return result;
}

/* log2 */
template <typename T> inline T log2(const T &x) { return Base::Math::log2(x); }

template <typename T> inline std::vector<T> log2(const std::vector<T> &vector) {
  std::vector<T> result;
  result.reserve(vector.size());
  for (const auto &element : vector) {
    result.push_back(Base::Math::log2(element));
  }
  return result;
}

template <typename T, std::size_t N>
inline std::array<T, N> log2(const std::array<T, N> &array) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::log2(array[i]);
  }
  return result;
}

/* log10 */
template <typename T> inline T log10(const T &x) {
  return Base::Math::log10(x);
}

template <typename T>
inline std::vector<T> log10(const std::vector<T> &vector) {
  std::vector<T> result;
  result.reserve(vector.size());
  for (const auto &element : vector) {
    result.push_back(Base::Math::log10(element));
  }
  return result;
}

template <typename T, std::size_t N>
inline std::array<T, N> log10(const std::array<T, N> &array) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::log10(array[i]);
  }
  return result;
}

/* pow */
template <typename T> inline T pow(const T &x, const T &y) {
  return Base::Math::pow(x, y);
}

template <typename T>
inline std::vector<T> pow(const std::vector<T> &vector_x, const T &y) {
  std::vector<T> result;
  result.reserve(vector_x.size());
  for (const auto &element : vector_x) {
    result.push_back(Base::Math::pow(element, y));
  }
  return result;
}

template <typename T>
inline std::vector<T> pow(const T &x, const std::vector<T> &vector_y) {
  std::vector<T> result;
  result.reserve(vector_y.size());
  for (const auto &element : vector_y) {
    result.push_back(Base::Math::pow(x, element));
  }
  return result;
}

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

template <typename T, std::size_t N>
inline std::array<T, N> pow(const std::array<T, N> &array_x, const T &y) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::pow(array_x[i], y);
  }
  return result;
}

template <typename T, std::size_t N>
inline std::array<T, N> pow(const T &x, const std::array<T, N> &array_y) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = Base::Math::pow(x, array_y[i]);
  }
  return result;
}

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

#endif // PYTHON_MATH_EXPONENTIAL_LOGARITHMIC_HPP
