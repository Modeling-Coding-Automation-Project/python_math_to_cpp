# 内部アルゴリズム: 算術・安全化基盤

本章では、`base_math` が他の関数群の土台として使う、算術処理と安全化処理を解説します。
`base_math` の公開インターフェース説明ではなく、**実際の計算戦略**に焦点を当てます。

## 1. 全体フロー

```mermaid
flowchart LR
    A[入力値] --> B[abs/mod などの基礎演算]
    B --> C[bitcast/clz/frexp/ldexp]
    C --> D[sqrt/exp/log/trig の上位関数]
    D --> E[avoid_zero_divide / near_zero で保護]
```

## 2. `abs` と `mod` の方針

### `abs`
- `BASE_MATH_USE_STD_MATH_` 有効時: `std::abs`
- それ以外: `x < 0 ? -x : x` の軽量分岐

### `mod` (`fmod`相当)
- 標準関数を使わないモードでは、商を整数キャストして `x - quotient * |y|` で剰余を計算
- `division_min` 未満の除数を抑制して極端な不安定計算を回避

## 3. ビット演算ユーティリティ（高速演算の基盤）

`base_math_arithmetic.hpp` では次を自前実装しています。

- `bitcast_u32` / `bitcast_f32` / `bitcast_u64` / `bitcast_f64`
  - `memcpy` による安全なビット再解釈
- `clz32` / `clz64`
  - 先頭ゼロカウント（fallback実装あり）
- `fast_frexpf`, `fast_ldexp_float`, `fast_ldexp_double`
  - IEEE754表現を直接操作して正規化・指数調整

これらは、`sqrt`・`exp`・`log` の近似法で、桁スケーリングや初期値生成に使われます。

## 4. `base_utility_cpp` の依存要点（本リポジトリ内再掲）

`base_math` は `external_libraries/base_utility_cpp/base_utility/base_utility.hpp` の次関数に依存します。

| 関数 | 役割 | 使用箇所の例 |
|---|---|---|
| `avoid_zero_divide(in, division_min)` | 分母が0近傍になるのを防ぐ | `tan`, `tanh` などの分母保護 |
| `near_zero(in, division_min)` | 0近傍判定 | `atan2` の象限分岐（x/yとも0近傍時など） |
| `sign(x)` | 符号判定（±1） | 一部近似の符号処理 |

### `avoid_zero_divide` の考え方

- `in` が `division_min` 未満で:
  - `in >= 0` なら `division_min` を返す
  - `-division_min < in < 0` なら `-division_min` を返す
- それ以外は `in` を返す

つまり「符号を維持したまま、絶対値が小さすぎる値を閾値まで押し上げる」ことで、ゼロ除算・過大ゲインを抑えます。

## 5. 設計上の意味

- 高速近似そのものより、**前処理（正規化）と安全化（ゼロ近傍保護）**が全体安定性を支える構造
- そのため、上位関数の誤差・発散挙動を読むときは、まず本章の補助関数の影響を確認するのが有効です
