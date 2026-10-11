# 参考にした記事（添付PDF由来）

本章は、添付された `math_1.pdf` 〜 `math_4.pdf` の先頭情報から、記事名と本リポジトリ実装との対応を整理したものです。
（※ サブモジュール内ドキュメントへのリンクは使用していません）

## 1. 記事一覧

| 添付PDF | 記事タイトル（PDF抽出） | URL |
|---|---|---|
| `math_1.pdf` | 平方根の高速計算アルゴリズムを設計する【sqrt】 | https://note.com/claude_a/n/n10c31c91844a |
| `math_2.pdf` | 指数、対数関数の高速計算アルゴリズムを設計する【exp, log】 | https://note.com/claude_a/n/n6d056e242b4b?magazine_key=m439fa71b3eef |
| `math_3.pdf` | 三角関数の高速計算アルゴリズムを設計する Part.1【sin, cos, tan】 | https://note.com/claude_a/n/n59491b5c5088?magazine_key=m439fa71b3eef |
| `math_4.pdf` | 三角関数の高速計算アルゴリズムを設計する Part.2【atan, asin, acos】 | https://note.com/claude_a/n/na4fd8d232e28?magazine_key=m439fa71b3eef |

## 2. 本実装との対応

| 記事トピック | 本リポジトリ対応箇所 |
|---|---|
| sqrt 高速化 | `base_math_exponential_logarithmic.hpp` の `fast_inverse_square_root`, `rsqrt_newton_method`, `sqrt_extraction`, `sqrt_newton_method` |
| exp/log 高速化 | 同ヘッダの `exp_maclaurin_expansion*`, `log_newton_method`, `log_maclaurin_expansion_with_table`, `pow*` |
| sin/cos/tan | `base_math_trigonometric.hpp` の `wrap_value_in_minus_pi_and_pi`, `sincos_maclaurin_expansion_with_DoubleAngleFormula`, `tan*` |
| atan/asin/acos | 同ヘッダの `atan_chebyshev*`, `atan2_chebyshev`, `asin_chebyshev`, `acos_chebyshev` |

## 3. ドキュメント反映方針

上記記事で扱われる主要アイデア（反復法、テーブル活用、近似式、レンジ縮約）を、本リポジトリの実装コードに対応づけて次章で解説しています。

- [指数・対数・平方根](./04_algorithm_exponential_logarithmic.md)
- [三角・逆三角・双曲線関数](./05_algorithm_trigonometric.md)
