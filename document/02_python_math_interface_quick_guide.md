# 公開インターフェース（python_math）クイックガイド

この章は「使い方だけを短く」まとめます。
内部アルゴリズムは [内部アルゴリズム章](./03_algorithm_arithmetic_and_safety.md) 以降を参照してください。

## 1. インクルード

```cpp
#include "python_math/python_math.hpp"
```

## 2. 基本（スカラ）

```cpp
double a = PythonMath::sqrt(2.0);
double b = PythonMath::exp(1.0);
double c = PythonMath::sin(0.5);
double d = PythonMath::atan2(1.0, 1.0);
```

## 3. `std::vector` / `std::array` への要素演算

```cpp
std::vector<double> v{0.0, 0.5, 1.0};
auto y = PythonMath::exp(v);   // 各要素にexpを適用

std::array<float, 3> a{{1.0f, 4.0f, 9.0f}};
auto r = PythonMath::sqrt(a);  // 各要素にsqrtを適用
```

## 4. 提供される主な関数

| 分類 | 関数 |
|---|---|
| 算術 | `abs`, `fmod` |
| 指数・対数 | `sqrt`, `exp`, `exp2`, `log`, `log2`, `log10`, `pow` |
| 三角・双曲線 | `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `atan2`, `sinh`, `cosh`, `tanh` |
| 定数 | `PI`, `TWO_PI`, `HALF_PI`, `TAU`, `E`, `LN_2`, `LN_10` |

## 5. 補足

- `python_math` の実体計算は `Base::Math::*` に委譲されます。
- そのため性能・精度特性は `base_math` の実装モード（マクロ設定）に依存します。
