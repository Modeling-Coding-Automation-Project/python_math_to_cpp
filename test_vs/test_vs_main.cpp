﻿#include <type_traits>
#include <iostream>
#include <cmath>

#include "base_math.hpp"
#include "python_math.hpp"

#include "MCAP_tester.hpp"


using namespace Tester;


template<typename T>
void check_base_math_exponential_logarithmic(void) {
    MCAPTester<T> tester;

#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    const T NEAR_LIMIT_SOFT = 1.0e-2F;
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    constexpr T NEAR_LIMIT_STRICT = std::is_same<T, double>::value ? T(1.0e-6) : T(1.0e-5);
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__


    std::vector<T> test_values_sqrt({
        static_cast<T>(-1),
        static_cast<T>(0),
        static_cast<T>(0.25),
        static_cast<T>(0.5),
        static_cast<T>(0.7777),
        static_cast<T>(1),
        static_cast<T>(2),
        static_cast<T>(3.4),
        static_cast<T>(5),
        static_cast<T>(11),
        static_cast<T>(56),
        static_cast<T>(102),
        static_cast<T>(1070),
        static_cast<T>(39900),
        static_cast<T>(608800)
    });

    /* sqrt */
    for (std::size_t i = 0; i < test_values_sqrt.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T sqrt_value = Base::Math::sqrt(test_values_sqrt[i]);
        T sqrt_answer = static_cast<T>(0);
        if (test_values_sqrt[i] < static_cast<T>(Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN)) {
            sqrt_answer = std::sqrt(static_cast<T>(Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN));
        }
        else {
            sqrt_answer = std::sqrt(test_values_sqrt[i]);
        }

        tester.expect_near(sqrt_value, sqrt_answer,
            NEAR_LIMIT_SOFT * std::abs(sqrt_answer) * static_cast<T>(10),
            "check sqrt.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T sqrt_value = Base::Math::sqrt(test_values_sqrt[i]);
        T sqrt_answer = static_cast<T>(0);
        if (test_values_sqrt[i] < static_cast<T>(Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN)) {
            sqrt_answer = std::sqrt(static_cast<T>(Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN));
        }
        else {
            sqrt_answer = std::sqrt(test_values_sqrt[i]);
        }

        tester.expect_near(sqrt_value, sqrt_answer, NEAR_LIMIT_STRICT * std::abs(sqrt_answer),
            "check sqrt.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

    /* sqrt check avoid zero divide */
    {
        T sqrt_d_value = Base::Math::sqrt_newton_method<T, Base::Math::SQRT_REPEAT_NUMBER>(
            static_cast<T>(2), static_cast<T>(3.123));
        T sqrt_d_answer = std::sqrt(static_cast<T>(3.123));

        tester.expect_near(sqrt_d_value, sqrt_d_answer, static_cast<T>(1.0e-1),
            "check sqrt division_min argument.");
    }

    /* rsqrt */
    for (std::size_t i = 0; i < test_values_sqrt.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T rsqrt_value = Base::Math::rsqrt(test_values_sqrt[i]);
        T rsqrt_answer = static_cast<T>(0);
        if (test_values_sqrt[i] < static_cast<T>(Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN)) {
            rsqrt_answer = static_cast<T>(1) / 
                std::sqrt(static_cast<T>(Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN));
        }
        else {
            rsqrt_answer = static_cast<T>(1) / std::sqrt(test_values_sqrt[i]);
        }

        tester.expect_near(rsqrt_value, rsqrt_answer,
            NEAR_LIMIT_SOFT * std::abs(rsqrt_answer) * static_cast<T>(10),
            "check rsqrt.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T rsqrt_value = Base::Math::rsqrt(test_values_sqrt[i]);
        T rsqrt_answer = static_cast<T>(0);
        if (test_values_sqrt[i] < static_cast<T>(Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN)) {
            rsqrt_answer = static_cast<T>(1) / 
                std::sqrt(static_cast<T>(Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN));
        }
        else {
            rsqrt_answer = static_cast<T>(1) / std::sqrt(test_values_sqrt[i]);
        }

        tester.expect_near(rsqrt_value, rsqrt_answer, NEAR_LIMIT_STRICT * std::abs(rsqrt_answer),
            "check rsqrt.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
}

    /* sqrt check avoid zero divide */
    {
        T sqrt_d_value = Base::Math::sqrt_newton_method<T, Base::Math::SQRT_REPEAT_NUMBER>(
            static_cast<T>(2), static_cast<T>(3.123));
        T sqrt_d_answer = std::sqrt(static_cast<T>(3.123));

        tester.expect_near(sqrt_d_value, sqrt_d_answer, static_cast<T>(1.0e-1),
            "check sqrt division_min argument.");
    }


    /* sqrt extraction */
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    for (std::size_t i = 0; i < test_values_sqrt.size(); i++) {
        T sqrt_value = static_cast<T>(0);
        T sqrt_answer = static_cast<T>(0);
        if (test_values_sqrt[i] < static_cast<T>(Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN)) {
            sqrt_value = Base::Math::sqrt_extraction<T, Base::Math::SQRT_EXTRACTION_REPEAT_NUMBER>(
                    static_cast<T>(Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN));
            sqrt_answer = std::sqrt(static_cast<T>(Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN));
        }
        else {
            sqrt_value = Base::Math::sqrt_extraction<T, Base::Math::SQRT_EXTRACTION_REPEAT_NUMBER>(
                    test_values_sqrt[i]);
            sqrt_answer = std::sqrt(test_values_sqrt[i]);
        }

        tester.expect_near(sqrt_value, sqrt_answer, NEAR_LIMIT_SOFT * std::abs(sqrt_answer),
            "check sqrt extraction double.");
    }
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    for (std::size_t i = 0; i < test_values_sqrt.size(); i++) {
        T sqrt_value = static_cast<T>(0);
        T sqrt_answer = static_cast<T>(0);
        if (test_values_sqrt[i] < static_cast<T>(Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN)) {
            sqrt_value = Base::Math::sqrt_extraction<T, Base::Math::SQRT_EXTRACTION_REPEAT_NUMBER>(
                    static_cast<T>(Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN));
            sqrt_answer = std::sqrt(static_cast<T>(Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN));
        }
        else {
            sqrt_value = Base::Math::sqrt_extraction<T, Base::Math::SQRT_EXTRACTION_REPEAT_NUMBER>(
                    test_values_sqrt[i]);
            sqrt_answer = std::sqrt(test_values_sqrt[i]);
        }

        tester.expect_near(sqrt_value, sqrt_answer, NEAR_LIMIT_STRICT * std::abs(sqrt_answer),
            "check sqrt extraction double.");
    }
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__


    std::vector<T> test_values_exp({
        static_cast<T>(-100),
        static_cast<T>(-87),
        static_cast<T>(-50),
        static_cast<T>(-11.113),
        static_cast<T>(-10),
        static_cast<T>(-1),
        static_cast<T>(0),
        static_cast<T>(2),
        static_cast<T>(10),
        static_cast<T>(11.113),
        static_cast<T>(50),
        static_cast<T>(87),
        static_cast<T>(100)
        });

    /* exp */
    for (std::size_t i = 0; i < test_values_exp.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T exp_value = Base::Math::exp(test_values_exp[i]);
        T exp_answer = static_cast<T>(0);
        if (test_values_exp[i] > static_cast<T>(Base::Math::EXP_INPUT_MAX)) {
            exp_answer = std::exp(static_cast<T>(Base::Math::EXP_INPUT_MAX));
        }
        else if (test_values_exp[i] < static_cast<T>(Base::Math::EXP_INPUT_MIN)) {
            exp_answer = std::exp(static_cast<T>(Base::Math::EXP_INPUT_MIN));
        }
        else {
            exp_answer = std::exp(test_values_exp[i]);
        }

        tester.expect_near(exp_value, exp_answer, (NEAR_LIMIT_SOFT * std::abs(exp_answer)),
            "check exp.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T exp_value = Base::Math::exp(test_values_exp[i]);
        T exp_answer = static_cast<T>(0);
        if (test_values_exp[i] > static_cast<T>(Base::Math::EXP_INPUT_MAX)) {
            exp_answer = std::exp(static_cast<T>(Base::Math::EXP_INPUT_MAX));
        }
        else if (test_values_exp[i] < static_cast<T>(Base::Math::EXP_INPUT_MIN)) {
            exp_answer = std::exp(static_cast<T>(Base::Math::EXP_INPUT_MIN));
        }
        else {
            exp_answer = std::exp(test_values_exp[i]);
        }

        tester.expect_near(exp_value, exp_answer, (NEAR_LIMIT_STRICT * std::abs(exp_answer)),
            "check exp.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

    /* exp Maclaurin Expansion with table double */
    for (std::size_t i = 0; i < test_values_exp.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T exp_value = Base::Math::exp_maclaurin_expansion_with_table<
            T, Base::Math::EXP_MACLAURIN_WITH_TABLE_REPEAT_NUMBER>(test_values_exp[i]);
        T exp_answer = static_cast<T>(0);
        if (test_values_exp[i] > static_cast<T>(Base::Math::EXP_INPUT_MAX)) {
            exp_answer = std::exp(static_cast<T>(Base::Math::EXP_INPUT_MAX));
        }
        else if (test_values_exp[i] < static_cast<T>(Base::Math::EXP_INPUT_MIN)) {
            exp_answer = std::exp(static_cast<T>(Base::Math::EXP_INPUT_MIN));
        }
        else {
            exp_answer = std::exp(test_values_exp[i]);
        }

        tester.expect_near(exp_value, exp_answer, (NEAR_LIMIT_SOFT * std::abs(exp_answer)),
            "check exp Maclaurin Expansion with table.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T exp_value = Base::Math::exp_maclaurin_expansion_with_table<
            T, Base::Math::EXP_MACLAURIN_WITH_TABLE_REPEAT_NUMBER>(test_values_exp[i]);
        T exp_answer = static_cast<T>(0);
        if (test_values_exp[i] > static_cast<T>(Base::Math::EXP_INPUT_MAX)) {
            exp_answer = std::exp(static_cast<T>(Base::Math::EXP_INPUT_MAX));
        }
        else if (test_values_exp[i] < static_cast<T>(Base::Math::EXP_INPUT_MIN)) {
            exp_answer = std::exp(static_cast<T>(Base::Math::EXP_INPUT_MIN));
        }
        else {
            exp_answer = std::exp(test_values_exp[i]);
        }

        tester.expect_near(exp_value, exp_answer, (NEAR_LIMIT_STRICT * std::abs(exp_answer)),
            "check exp Maclaurin Expansion with table.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }


    std::vector<T> test_values_exp2({
    static_cast<T>(-200),
    static_cast<T>(-126),
    static_cast<T>(-50),
    static_cast<T>(-8.445),
    static_cast<T>(-1),
    static_cast<T>(-0.777),
    static_cast<T>(0),
    static_cast<T>(0.777),
    static_cast<T>(1),
    static_cast<T>(8.445),
    static_cast<T>(50),
    static_cast<T>(126),
    static_cast<T>(200)
        });


    /* exp2 */
    for (std::size_t i = 0; i < test_values_exp2.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T exp2_value = Base::Math::exp2(test_values_exp2[i]);
        T exp2_answer = static_cast<T>(0);
        if (test_values_exp2[i] > static_cast<T>(Base::Math::EXP2_INPUT_MAX)) {
            exp2_answer = std::exp2(static_cast<T>(Base::Math::EXP2_INPUT_MAX));
        }
        else if (test_values_exp2[i] < static_cast<T>(Base::Math::EXP2_INPUT_MIN)) {
            exp2_answer = std::exp2(static_cast<T>(Base::Math::EXP2_INPUT_MIN));
        }
        else {
            exp2_answer = std::exp2(test_values_exp2[i]);
        }

        tester.expect_near(exp2_value, exp2_answer, (NEAR_LIMIT_SOFT * std::abs(exp2_answer)),
            "check exp2.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T exp2_value = Base::Math::exp2(test_values_exp2[i]);
        T exp2_answer = static_cast<T>(0);
        if (test_values_exp2[i] > static_cast<T>(Base::Math::EXP2_INPUT_MAX)) {
            exp2_answer = std::exp2(static_cast<T>(Base::Math::EXP2_INPUT_MAX));
        }
        else if (test_values_exp2[i] < static_cast<T>(Base::Math::EXP2_INPUT_MIN)) {
            exp2_answer = std::exp2(static_cast<T>(Base::Math::EXP2_INPUT_MIN));
        }
        else {
            exp2_answer = std::exp2(test_values_exp2[i]);
        }

        tester.expect_near(exp2_value, exp2_answer, (NEAR_LIMIT_STRICT * std::abs(exp2_answer)),
            "check exp2.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }


    std::vector<T> test_values_log({
    static_cast<T>(-1),
    static_cast<T>(0.001),
    static_cast<T>(0.1),
    static_cast<T>(0.5),
    static_cast<T>(0.77777),
    static_cast<T>(1),
    static_cast<T>(2),
    static_cast<T>(5),
    static_cast<T>(7.77777),
    static_cast<T>(10),
    static_cast<T>(50),
    static_cast<T>(111),
    static_cast<T>(500),
    static_cast<T>(1020),
    static_cast<T>(14350)
        });

    /* log */
    for (std::size_t i = 0; i < test_values_log.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T log_value = Base::Math::log(test_values_log[i]);
        T log_answer = static_cast<T>(0);
        if (test_values_log[i] <= static_cast<T>(0)) {
            log_value = static_cast<T>(Base::Math::LOG_OUTPUT_MIN);
            log_answer = static_cast<T>(Base::Math::LOG_OUTPUT_MIN);
        }
        else {
            log_answer = std::log(test_values_log[i]);
        }

        tester.expect_near(log_value, log_answer, NEAR_LIMIT_SOFT * std::abs(log_answer),
            "check log.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T log_value = Base::Math::log(test_values_log[i]);
        T log_answer = static_cast<T>(0);
        if (test_values_log[i] <= static_cast<T>(0)) {
            log_answer = static_cast<T>(Base::Math::LOG_OUTPUT_MIN);
        }
        else {
            log_answer = std::log(test_values_log[i]);
        }

        tester.expect_near(log_value, log_answer, NEAR_LIMIT_STRICT * std::abs(log_answer),
            "check log.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

    /* log maclaurin expansion with table */
    for (std::size_t i = 0; i < test_values_log.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T log_value = Base::Math::log_maclaurin_expansion_with_table(test_values_log[i]);
        T log_answer = static_cast<T>(0);
        if (test_values_log[i] <= static_cast<T>(0)) {
            log_answer = static_cast<T>(Base::Math::LOG_OUTPUT_MIN);
        }
        else {
            log_answer = std::log(test_values_log[i]);
        }

        tester.expect_near(log_value, log_answer, NEAR_LIMIT_SOFT * std::abs(log_answer),
            "check log maclaurin expansion with table.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T log_value = Base::Math::log_maclaurin_expansion_with_table(test_values_log[i]);
        T log_answer = static_cast<T>(0);
        if (test_values_log[i] <= static_cast<T>(0)) {
            log_answer = static_cast<T>(Base::Math::LOG_OUTPUT_MIN);
        }
        else {
            log_answer = std::log(test_values_log[i]);
        }

        tester.expect_near(log_value, log_answer, NEAR_LIMIT_STRICT * std::abs(log_answer),
            "check log maclaurin expansion with table.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

    /* log2 */
    for (std::size_t i = 0; i < test_values_log.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T log2_value = Base::Math::log2(test_values_log[i]);
        T log2_answer = static_cast<T>(0);
        if (test_values_log[i] <= static_cast<T>(0)) {
            log2_value = static_cast<T>(Base::Math::LOG_OUTPUT_MIN)
                / static_cast<T>(Base::Math::LN_2);
            log2_answer = static_cast<T>(Base::Math::LOG_OUTPUT_MIN)
                / static_cast<T>(Base::Math::LN_2);
        }
        else {
            log2_answer = std::log2(test_values_log[i]);
        }

        tester.expect_near(log2_value, log2_answer, NEAR_LIMIT_SOFT * std::abs(log2_answer),
            "check log2.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T log2_value = Base::Math::log2(test_values_log[i]);
        T log2_answer = static_cast<T>(0);
        if (test_values_log[i] <= static_cast<T>(0)) {
            log2_answer = static_cast<T>(Base::Math::LOG_OUTPUT_MIN)
                / static_cast<T>(Base::Math::LN_2);
        }
        else {
            log2_answer = std::log2(test_values_log[i]);
        }

        tester.expect_near(log2_value, log2_answer, NEAR_LIMIT_STRICT * std::abs(log2_answer),
            "check log2.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

    /* log10 */
    for (std::size_t i = 0; i < test_values_log.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T log10_value = Base::Math::log10(test_values_log[i]);
        T log10_answer = static_cast<T>(0);
        if (test_values_log[i] <= static_cast<T>(0)) {
            log10_value = static_cast<T>(Base::Math::LOG_OUTPUT_MIN)
                / static_cast<T>(Base::Math::LN_10);
            log10_answer = static_cast<T>(Base::Math::LOG_OUTPUT_MIN)
                / static_cast<T>(Base::Math::LN_10);
        }
        else {
            log10_answer = std::log10(test_values_log[i]);
        }

        tester.expect_near(log10_value, log10_answer, NEAR_LIMIT_SOFT * std::abs(log10_answer),
            "check log10.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T log10_value = Base::Math::log10(test_values_log[i]);
        T log10_answer = static_cast<T>(0);
        if (test_values_log[i] <= static_cast<T>(0)) {
            log10_answer = static_cast<T>(Base::Math::LOG_OUTPUT_MIN)
                / static_cast<T>(Base::Math::LN_10);
        }
        else {
            log10_answer = std::log10(test_values_log[i]);
        }

        tester.expect_near(log10_value, log10_answer, NEAR_LIMIT_STRICT * std::abs(log10_answer),
            "check log10.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

    std::vector<std::pair<T, T>> test_values_pow({
        {static_cast<T>(0.1), static_cast<T>(0.1)},
        {static_cast<T>(0.1), static_cast<T>(0.5)},
        {static_cast<T>(0.1), static_cast<T>(1)},
        {static_cast<T>(0.1), static_cast<T>(2)},
        {static_cast<T>(0.5), static_cast<T>(0.1)},
        {static_cast<T>(0.5), static_cast<T>(0.5)},
        {static_cast<T>(0.5), static_cast<T>(1)},
        {static_cast<T>(0.5), static_cast<T>(2)},
        {static_cast<T>(1), static_cast<T>(0.1)},
        {static_cast<T>(1), static_cast<T>(0.5)},
        {static_cast<T>(1), static_cast<T>(1)},
        {static_cast<T>(1), static_cast<T>(2)},
        {static_cast<T>(2), static_cast<T>(0.1)},
        {static_cast<T>(2), static_cast<T>(0.5)},
        {static_cast<T>(2), static_cast<T>(1)},
        {static_cast<T>(2), static_cast<T>(2)},
        {static_cast<T>(5), static_cast<T>(0.1)},
        {static_cast<T>(5), static_cast<T>(0.5)},
        {static_cast<T>(5), static_cast<T>(1)},
        {static_cast<T>(5), static_cast<T>(2)},
        {static_cast<T>(10), static_cast<T>(0.1)},
        {static_cast<T>(10), static_cast<T>(0.5)},
        {static_cast<T>(10), static_cast<T>(1)},
        {static_cast<T>(10), static_cast<T>(2)},
        {static_cast<T>(50), static_cast<T>(0.1)},
        {static_cast<T>(50), static_cast<T>(0.5)},
        {static_cast<T>(50), static_cast<T>(1)},
        {static_cast<T>(50), static_cast<T>(2)},
        {static_cast<T>(100), static_cast<T>(0.1)},
        {static_cast<T>(100), static_cast<T>(0.5)},
        {static_cast<T>(100), static_cast<T>(1)},
        {static_cast<T>(100), static_cast<T>(2)},
    });

    /* pow */
    for (std::size_t i = 0; i < test_values_pow.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T pow_value = Base::Math::pow(test_values_pow[i].first, test_values_pow[i].second);
        T pow_answer = std::pow(test_values_pow[i].first, test_values_pow[i].second);

        tester.expect_near(pow_value, pow_answer, NEAR_LIMIT_SOFT * std::abs(pow_answer),
            "check pow.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T pow_value = Base::Math::pow(test_values_pow[i].first, test_values_pow[i].second);
        T pow_answer = std::pow(test_values_pow[i].first, test_values_pow[i].second);

        tester.expect_near(pow_value, pow_answer, NEAR_LIMIT_STRICT * std::abs(pow_answer),
            "check pow.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }


    tester.throw_error_if_test_failed();
}


template<typename T>
void check_base_math_trigonometric(void) {

    MCAPTester<T> tester;

#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    const T NEAR_LIMIT_SOFT = 1.0e-2F;
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    constexpr T NEAR_LIMIT_STRICT = std::is_same<T, double>::value ? T(1.0e-6) : T(1.0e-5);
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__


    std::vector<T> test_values_sin = {
        static_cast<T>(-10 * Base::Math::PI),
        static_cast<T>(-3 * Base::Math::PI),
        static_cast<T>(-2 * Base::Math::PI),
        static_cast<T>(-3 * Base::Math::PI / 2),
        static_cast<T>(-Base::Math::PI),
        static_cast<T>(-Base::Math::PI / 2),
        static_cast<T>(0),
        static_cast<T>(Base::Math::PI / 2),
        static_cast<T>(Base::Math::PI),
        static_cast<T>(3 * Base::Math::PI / 2),
        static_cast<T>(2 * Base::Math::PI),
        static_cast<T>(-5 * Base::Math::PI / 4),
        static_cast<T>(1),
        static_cast<T>(5 * Base::Math::PI / 4),
        static_cast<T>(-7 * Base::Math::PI / 4),
        static_cast<T>(7 * Base::Math::PI / 4),
        static_cast<T>(3 * Base::Math::PI),
        static_cast<T>(10 * Base::Math::PI)
    };

    /* sin */
    for (std::size_t i = 0; i < test_values_sin.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T sin_value = Base::Math::sin(test_values_sin[i]);
        T sin_answer = std::sin(test_values_sin[i]);

        tester.expect_near(sin_value, sin_answer, NEAR_LIMIT_SOFT,
            "check sin.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T sin_value = Base::Math::sin(test_values_sin[i]);
        T sin_answer = std::sin(test_values_sin[i]);

        tester.expect_near(sin_value, sin_answer, NEAR_LIMIT_STRICT,
            "check sin.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

    /* cos */
    for (std::size_t i = 0; i < test_values_sin.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T cos_value = Base::Math::cos(test_values_sin[i]);
        T cos_answer = std::cos(test_values_sin[i]);

        tester.expect_near(cos_value, cos_answer, NEAR_LIMIT_SOFT,
            "check cos.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T cos_value = Base::Math::cos(test_values_sin[i]);
        T cos_answer = std::cos(test_values_sin[i]);

        tester.expect_near(cos_value, cos_answer, NEAR_LIMIT_STRICT,
            "check cos.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

    /* cos maclaurin expansion with DoubleAngleFormula */
    for (std::size_t i = 0; i < test_values_sin.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T cos_value = Base::Math::cos_maclaurin_expansion_with_DoubleAngleFormula<
            T, Base::Math::COS_MACLAURIN_DOUBLEANGLE_REPEAT_NUMBER>(test_values_sin[i]);
        T cos_answer = std::cos(test_values_sin[i]);

        tester.expect_near(cos_value, cos_answer, NEAR_LIMIT_SOFT,
            "check cos maclaurin expansion with DoubleAngleFormula.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T cos_value = Base::Math::cos_maclaurin_expansion_with_DoubleAngleFormula<
            T, Base::Math::COS_MACLAURIN_DOUBLEANGLE_REPEAT_NUMBER>(test_values_sin[i]);
        T cos_answer = std::cos(test_values_sin[i]);

        tester.expect_near(cos_value, cos_answer, NEAR_LIMIT_STRICT,
            "check cos maclaurin expansion with DoubleAngleFormula.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

    /* sin cos Maclaurin expansion with DoubleAngleFormula */
    for (std::size_t i = 0; i < test_values_sin.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T sin_value = static_cast<T>(0);
        T cos_value = static_cast<T>(0);
        Base::Math::sincos_maclaurin_expansion_with_DoubleAngleFormula<
            T, Base::Math::SINCOS_MACLAURIN_DOUBLEANGLE_REPEAT_NUMBER>(
            test_values_sin[i], sin_value, cos_value);
    
        T sin_answer = std::sin(test_values_sin[i]);
        T cos_answer = std::cos(test_values_sin[i]);

        tester.expect_near(sin_value, sin_answer, NEAR_LIMIT_SOFT,
            "check sin CORDIC for floating point number.");
        tester.expect_near(cos_value, cos_answer, NEAR_LIMIT_SOFT,
            "check cos CORDIC for floating point number.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T sin_value = static_cast<T>(0);
        T cos_value = static_cast<T>(0);
        Base::Math::sincos_maclaurin_expansion_with_DoubleAngleFormula<
            T, Base::Math::SINCOS_MACLAURIN_DOUBLEANGLE_REPEAT_NUMBER>(
                test_values_sin[i], sin_value, cos_value);

        T sin_answer = std::sin(test_values_sin[i]);
        T cos_answer = std::cos(test_values_sin[i]);

        tester.expect_near(sin_value, sin_answer, NEAR_LIMIT_STRICT,
            "check sincos maclaurin expansion with DoubleAngleFormula.");
        tester.expect_near(cos_value, cos_answer, NEAR_LIMIT_STRICT,
            "check sincos maclaurin expansion with DoubleAngleFormula.");

        sin_value = static_cast<T>(0);
        cos_value = static_cast<T>(0);

        Base::Math::sincos(test_values_sin[i], sin_value, cos_value);

        tester.expect_near(sin_value, sin_answer, NEAR_LIMIT_STRICT,
            "check sincos.");

#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

    /* sin maclaurin expansion with DoubleAngleFormula */
    for (std::size_t i = 0; i < test_values_sin.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T sin_value = Base::Math::sin_maclaurin_expansion_with_DoubleAngleFormula<
            T, Base::Math::SIN_MACLAURIN_DOUBLEANGLE_REPEAT_NUMBER>(test_values_sin[i]);
        T sin_answer = std::sin(test_values_sin[i]);

        tester.expect_near(sin_value, sin_answer, NEAR_LIMIT_SOFT,
            "check sin maclaurin expansion with DoubleAngleFormula.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T sin_value = Base::Math::sin_maclaurin_expansion_with_DoubleAngleFormula<
            T, Base::Math::SIN_MACLAURIN_DOUBLEANGLE_REPEAT_NUMBER>(test_values_sin[i]);
        T sin_answer = std::sin(test_values_sin[i]);

        tester.expect_near(sin_value, sin_answer, NEAR_LIMIT_STRICT,
            "check sin maclaurin expansion with DoubleAngleFormula.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

    /* tan maclaurin expansion with DoubleAngleFormula */
    for (std::size_t i = 0; i < test_values_sin.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T tan_value = Base::Math::tan_maclaurin_expansion_with_DoubleAngleFormula<
            T, Base::Math::SINCOS_MACLAURIN_DOUBLEANGLE_REPEAT_NUMBER>(test_values_sin[i]);
        T tan_answer = std::tan(test_values_sin[i]);

        if (std::abs(tan_value) > static_cast<T>(1 / NEAR_LIMIT_SOFT) &&
            std::abs(tan_answer) > static_cast<T>(1 / NEAR_LIMIT_SOFT)) {
            /* Do Nothing. */
    }
        else {
            tester.expect_near(tan_value, tan_answer, NEAR_LIMIT_SOFT,
                "check tan maclaurin expansion with DoubleAngleFormula.");
        }
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T tan_value = Base::Math::tan_maclaurin_expansion_with_DoubleAngleFormula<
            T, Base::Math::SINCOS_MACLAURIN_DOUBLEANGLE_REPEAT_NUMBER>(test_values_sin[i]);
        T tan_answer = std::tan(test_values_sin[i]);

        if (std::abs(tan_value) > static_cast<T>(1 / NEAR_LIMIT_STRICT) &&
            std::abs(tan_answer) > static_cast<T>(1 / NEAR_LIMIT_STRICT)) {
            /* Do Nothing. */
        }
        else {
            tester.expect_near(tan_value, tan_answer, NEAR_LIMIT_STRICT,
                "check tan maclaurin expansion with DoubleAngleFormula.");
        }
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

    /* tan */
    for (std::size_t i = 0; i < test_values_sin.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T tan_value = Base::Math::tan(test_values_sin[i]);
        T tan_answer = std::tan(test_values_sin[i]);

        if (std::abs(tan_value) > static_cast<T>(1 / NEAR_LIMIT_SOFT) &&
            std::abs(tan_answer) > static_cast<T>(1 / NEAR_LIMIT_SOFT)) {
            /* Do Nothing. */
        }
        else {
            tester.expect_near(tan_value, tan_answer, NEAR_LIMIT_SOFT,
                "check tan.");
        }
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T tan_value = Base::Math::tan(test_values_sin[i]);
        T tan_answer = std::tan(test_values_sin[i]);

        if (std::abs(tan_value) > static_cast<T>(1 / NEAR_LIMIT_STRICT) &&
            std::abs(tan_answer) > static_cast<T>(1 / NEAR_LIMIT_STRICT)) {
            /* Do Nothing. */
        }
        else {
            tester.expect_near(tan_value, tan_answer, NEAR_LIMIT_STRICT,
                "check tan.");
        }
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }


    std::vector<T> test_values_atan = {
    static_cast<T>(-10),
    static_cast<T>(-3),
    static_cast<T>(-1.1),
    static_cast<T>(-1),
    static_cast<T>(-0.9),
    static_cast<T>(-0.5),
    static_cast<T>(0),
    static_cast<T>(0.5),
    static_cast<T>(0.9),
    static_cast<T>(1),
    static_cast<T>(1.1),
    static_cast<T>(3),
    static_cast<T>(10)
    };

    /* atan */
    for (std::size_t i = 0; i < test_values_atan.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T atan_value = Base::Math::atan(test_values_atan[i]);
        T atan_answer = std::atan(test_values_atan[i]);

        tester.expect_near(atan_value, atan_answer, NEAR_LIMIT_SOFT,
            "check atan.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T atan_value = Base::Math::atan(test_values_atan[i]);
        T atan_answer = std::atan(test_values_atan[i]);

        tester.expect_near(atan_value, atan_answer, NEAR_LIMIT_STRICT,
            "check atan.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

    std::vector<std::pair<T, T>> test_values_atan2({
            {static_cast<T>(0), static_cast<T>(0)},    /* 原点 */
            {static_cast<T>(1), static_cast<T>(0)},    /* 正のx軸 */
            {static_cast<T>(0), static_cast<T>(1)},    /* 正のy軸 */
            {static_cast<T>(-1), static_cast<T>(0)},   /* 負のx軸 */
            {static_cast<T>(0), static_cast<T>(-1)},   /* 負のy軸 */
            {static_cast<T>(1), static_cast<T>(1)},    /* 第1象限 */
            {static_cast<T>(-1), static_cast<T>(1)},   /* 第2象限 */
            {static_cast<T>(-1), static_cast<T>(-1)},  /* 第3象限 */
            {static_cast<T>(1), static_cast<T>(-1)},   /* 第4象限 */
            {static_cast<T>(3), static_cast<T>(4)},    /* 一般的な正のケース */
            {static_cast<T>(-3), static_cast<T>(4)},   /* 第2象限 */
            {static_cast<T>(-3), static_cast<T>(-4)},  /* 第3象限 */
            {static_cast<T>(3), static_cast<T>(-4)},   /* 第4象限 */
            {static_cast<T>(0.5), static_cast<T>(0.5)}, /* 小さい正の値 */
            {static_cast<T>(-0.5), static_cast<T>(0.5)}, /* 小さい負の値 */
            {static_cast<T>(1e-6), static_cast<T>(1e-6)}, /* 非常に小さい正の値 */
            {static_cast<T>(-1e-6), static_cast<T>(1e-6)}, /* 非常に小さい負の値 */
            {static_cast<T>(1), static_cast<T>(1e6)},  /* 大きな差があるケース */
            {static_cast<T>(1e6), static_cast<T>(1)},  /* 大きな差が逆になるケース */
            {static_cast<T>(-1), static_cast<T>(1e6)}, /* 異なる象限のケース */
        });

    /* atan2 */
    for (std::size_t i = 0; i < test_values_atan2.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T atan2_value = Base::Math::atan2(test_values_atan2[i].second, test_values_atan2[i].first);
        T atan2_answer = std::atan2(test_values_atan2[i].second, test_values_atan2[i].first);

        tester.expect_near(atan2_value, atan2_answer, NEAR_LIMIT_SOFT,
            "check atan2.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T atan2_value = Base::Math::atan2(test_values_atan2[i].second, test_values_atan2[i].first);
        T atan2_answer = std::atan2(test_values_atan2[i].second, test_values_atan2[i].first);

        tester.expect_near(atan2_value, atan2_answer, NEAR_LIMIT_STRICT,
            "check atan2.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

    std::vector<T> test_values_asin_acos = {
        static_cast<T>(-2),
        static_cast<T>(-1),
        static_cast<T>(-0.86602540378), // -sqrt(3)/2
        static_cast<T>(-0.70710678118), // -sqrt(2)/2
        static_cast<T>(-0.5),
        static_cast<T>(-0.25881904510), // -1/2 * sqrt(3)/2
        static_cast<T>(0),
        static_cast<T>(0.25881904510),  // 1/2 * sqrt(3)/2
        static_cast<T>(0.5),
        static_cast<T>(0.70710678118),  // sqrt(2)/2
        static_cast<T>(0.86602540378),  // sqrt(3)/2
        static_cast<T>(1),
        static_cast<T>(2)
    };

    /* asin */
    for (std::size_t i = 0; i < test_values_asin_acos.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T asin_value = Base::Math::asin(test_values_asin_acos[i]);

        T asin_answer = static_cast<T>(0);
        if (test_values_asin_acos[i] < static_cast<T>(-1)) {
            asin_answer = std::asin(static_cast<T>(-1));
        }
        else if (test_values_asin_acos[i] > static_cast<T>(1)) {
            asin_answer = std::asin(static_cast<T>(1));
        }
        else
        {
            asin_answer = std::asin(test_values_asin_acos[i]);
        }

        tester.expect_near(asin_value, asin_answer, NEAR_LIMIT_SOFT,
            "check asin.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T asin_value = Base::Math::asin(test_values_asin_acos[i]);

        T asin_answer = static_cast<T>(0);
        if (test_values_asin_acos[i] < static_cast<T>(-1)) {
            asin_answer = std::asin(static_cast<T>(-1));
        }
        else if (test_values_asin_acos[i] > static_cast<T>(1)) {
            asin_answer = std::asin(static_cast<T>(1));
        }
        else
        {
            asin_answer = std::asin(test_values_asin_acos[i]);
        }

        tester.expect_near(asin_value, asin_answer, NEAR_LIMIT_STRICT,
            "check asin.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

    /* acos */
    for (std::size_t i = 0; i < test_values_asin_acos.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T acos_value = Base::Math::acos(test_values_asin_acos[i]);

        T acos_answer = static_cast<T>(0);
        if (test_values_asin_acos[i] < static_cast<T>(-1)) {
            acos_answer = std::acos(static_cast<T>(-1));
        }
        else if (test_values_asin_acos[i] > static_cast<T>(1)) {
            acos_answer = std::acos(static_cast<T>(1));
        }
        else
        {
            acos_answer = std::acos(test_values_asin_acos[i]);
        }

        tester.expect_near(acos_value, acos_answer, NEAR_LIMIT_SOFT,
            "check acos.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T acos_value = Base::Math::acos(test_values_asin_acos[i]);

        T acos_answer = static_cast<T>(0);
        if (test_values_asin_acos[i] < static_cast<T>(-1)) {
            acos_answer = std::acos(static_cast<T>(-1));
        }
        else if (test_values_asin_acos[i] > static_cast<T>(1)) {
            acos_answer = std::acos(static_cast<T>(1));
        }
        else
        {
            acos_answer = std::acos(test_values_asin_acos[i]);
        }

        tester.expect_near(acos_value, acos_answer, NEAR_LIMIT_STRICT,
            "check acos.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

    std::vector<T> test_values_hyperbolic({
        static_cast<T>(-10),
        static_cast<T>(-8.6),
        static_cast<T>(-5),
        static_cast<T>(-1),
        static_cast<T>(-0.1),
        static_cast<T>(0),
        static_cast<T>(0.1),
        static_cast<T>(1),
        static_cast<T>(5),
        static_cast<T>(8.7),
        static_cast<T>(10)
        });

    /* sinh */
    for (std::size_t i = 0; i < test_values_hyperbolic.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T sinh_value = Base::Math::sinh(test_values_hyperbolic[i]);
        T sinh_answer = std::sinh(test_values_hyperbolic[i]);

        tester.expect_near(sinh_value, sinh_answer, NEAR_LIMIT_SOFT * std::abs(sinh_answer),
            "check sinh.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T sinh_value = Base::Math::sinh(test_values_hyperbolic[i]);
        T sinh_answer = std::sinh(test_values_hyperbolic[i]);

        tester.expect_near(sinh_value, sinh_answer, NEAR_LIMIT_STRICT * std::abs(sinh_answer),
            "check sinh.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

    /* sinh maclaurin expansion */
    for (std::size_t i = 0; i < test_values_hyperbolic.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T sinh_value = Base::Math::sinh_maclaurin_expansion<
            T, Base::Math::EXP_REPEAT_NUMBER>(test_values_hyperbolic[i]);
        T sinh_answer = std::sinh(test_values_hyperbolic[i]);

        tester.expect_near(sinh_value, sinh_answer, NEAR_LIMIT_SOFT * std::abs(sinh_answer),
            "check sinh maclaurin expansion.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T sinh_value = Base::Math::sinh_maclaurin_expansion<
            T, Base::Math::EXP_REPEAT_NUMBER>(test_values_hyperbolic[i]);
        T sinh_answer = std::sinh(test_values_hyperbolic[i]);

        tester.expect_near(sinh_value, sinh_answer, NEAR_LIMIT_STRICT * std::abs(sinh_answer),
            "check sinh maclaurin expansion.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

    /* sinh maclaurin expansion with table */
    for (std::size_t i = 0; i < test_values_hyperbolic.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T sinh_value = Base::Math::sinh_maclaurin_expansion_with_table<
            T, Base::Math::EXP_MACLAURIN_WITH_TABLE_REPEAT_NUMBER>(test_values_hyperbolic[i]);
        T sinh_answer = std::sinh(test_values_hyperbolic[i]);

        tester.expect_near(sinh_value, sinh_answer, NEAR_LIMIT_SOFT * std::abs(sinh_answer),
            "check sinh maclaurin expansion with table.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T sinh_value = Base::Math::sinh_maclaurin_expansion_with_table<
            T, Base::Math::EXP_MACLAURIN_WITH_TABLE_REPEAT_NUMBER>(test_values_hyperbolic[i]);
        T sinh_answer = std::sinh(test_values_hyperbolic[i]);

        tester.expect_near(sinh_value, sinh_answer, NEAR_LIMIT_STRICT * std::abs(sinh_answer),
            "check sinh maclaurin expansion with table.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

    /* cosh */
    for (std::size_t i = 0; i < test_values_hyperbolic.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T cosh_value = Base::Math::cosh(test_values_hyperbolic[i]);
        T cosh_answer = std::cosh(test_values_hyperbolic[i]);

        tester.expect_near(cosh_value, cosh_answer, NEAR_LIMIT_SOFT * std::abs(cosh_answer),
            "check cosh.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T cosh_value = Base::Math::cosh(test_values_hyperbolic[i]);
        T cosh_answer = std::cosh(test_values_hyperbolic[i]);

        tester.expect_near(cosh_value, cosh_answer, NEAR_LIMIT_STRICT * std::abs(cosh_answer),
            "check cosh.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

    /* cosh maclaurin expansion */
    for (std::size_t i = 0; i < test_values_hyperbolic.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T cosh_value = Base::Math::cosh_maclaurin_expansion<
            T, Base::Math::EXP_REPEAT_NUMBER>(test_values_hyperbolic[i]);
        T cosh_answer = std::cosh(test_values_hyperbolic[i]);

        tester.expect_near(cosh_value, cosh_answer, NEAR_LIMIT_SOFT * std::abs(cosh_answer),
            "check cosh maclaurin expansion.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T cosh_value = Base::Math::cosh_maclaurin_expansion<
            T, Base::Math::EXP_REPEAT_NUMBER>(test_values_hyperbolic[i]);
        T cosh_answer = std::cosh(test_values_hyperbolic[i]);

        tester.expect_near(cosh_value, cosh_answer, NEAR_LIMIT_STRICT * std::abs(cosh_answer),
            "check cosh maclaurin expansion.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

    /* cosh maclaurin expansion with table */
    for (std::size_t i = 0; i < test_values_hyperbolic.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T cosh_value = Base::Math::cosh_maclaurin_expansion_with_table<
            T, Base::Math::EXP_MACLAURIN_WITH_TABLE_REPEAT_NUMBER>(test_values_hyperbolic[i]);
        T cosh_answer = std::cosh(test_values_hyperbolic[i]);

        tester.expect_near(cosh_value, cosh_answer, NEAR_LIMIT_SOFT * std::abs(cosh_answer),
            "check cosh maclaurin expansion with table.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T cosh_value = Base::Math::cosh_maclaurin_expansion_with_table<
            T, Base::Math::EXP_MACLAURIN_WITH_TABLE_REPEAT_NUMBER>(test_values_hyperbolic[i]);
        T cosh_answer = std::cosh(test_values_hyperbolic[i]);

        tester.expect_near(cosh_value, cosh_answer, NEAR_LIMIT_STRICT * std::abs(cosh_answer),
            "check cosh maclaurin expansion with table.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

    /* tanh */
    for (std::size_t i = 0; i < test_values_hyperbolic.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T tanh_value = Base::Math::tanh(test_values_hyperbolic[i]);
        T tanh_answer = std::tanh(test_values_hyperbolic[i]);

        tester.expect_near(tanh_value, tanh_answer, NEAR_LIMIT_SOFT,
            "check tanh.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T tanh_value = Base::Math::tanh(test_values_hyperbolic[i]);
        T tanh_answer = std::tanh(test_values_hyperbolic[i]);

        tester.expect_near(tanh_value, tanh_answer, NEAR_LIMIT_STRICT,
            "check tanh.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

    /* tanh maclaurin expansion */
    for (std::size_t i = 0; i < test_values_hyperbolic.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T tanh_value = Base::Math::tanh_maclaurin_expansion<
            T, Base::Math::EXP_REPEAT_NUMBER>(test_values_hyperbolic[i]);
        T tanh_answer = std::tanh(test_values_hyperbolic[i]);

        tester.expect_near(tanh_value, tanh_answer, NEAR_LIMIT_SOFT * std::abs(tanh_answer),
            "check tanh maclaurin expansion.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T tanh_value = Base::Math::tanh_maclaurin_expansion<
            T, Base::Math::EXP_REPEAT_NUMBER>(test_values_hyperbolic[i]);
        T tanh_answer = std::tanh(test_values_hyperbolic[i]);

        tester.expect_near(tanh_value, tanh_answer, NEAR_LIMIT_STRICT * std::abs(tanh_answer),
            "check tanh maclaurin expansion.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

    /* tanh maclaurin expansion with table */
    for (std::size_t i = 0; i < test_values_hyperbolic.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T tanh_value = Base::Math::tanh_maclaurin_expansion_with_table<
            T, Base::Math::EXP_MACLAURIN_WITH_TABLE_REPEAT_NUMBER>(test_values_hyperbolic[i]);
        T tanh_answer = std::tanh(test_values_hyperbolic[i]);

        tester.expect_near(tanh_value, tanh_answer, NEAR_LIMIT_SOFT * std::abs(tanh_answer),
            "check tanh maclaurin expansion with table.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T tanh_value = Base::Math::tanh_maclaurin_expansion_with_table<
            T, Base::Math::EXP_MACLAURIN_WITH_TABLE_REPEAT_NUMBER>(test_values_hyperbolic[i]);
        T tanh_answer = std::tanh(test_values_hyperbolic[i]);

        tester.expect_near(tanh_value, tanh_answer, NEAR_LIMIT_STRICT * std::abs(tanh_answer),
            "check tanh maclaurin expansion with table.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }


    tester.throw_error_if_test_failed();
}

template<typename T>
void check_base_math_calc(void) {

    check_base_math_exponential_logarithmic<T>();

    check_base_math_trigonometric<T>();
}

template<typename T>
void check_python_math_arithmetic(void) {

    MCAPTester<T> tester;

    constexpr T NEAR_LIMIT_STRICT = std::is_same<T, double>::value ? T(1.0e-6) : T(1.0e-5);


    std::vector<T> test_values_abs({
        static_cast<T>(-10),
        static_cast<T>(-8.6),
        static_cast<T>(-5),
        static_cast<T>(-1),
        static_cast<T>(-0.1),
        static_cast<T>(0),
        static_cast<T>(0.1),
        static_cast<T>(1),
        static_cast<T>(5),
        static_cast<T>(8.7),
        static_cast<T>(10)
        });

    /* abs */
    for (std::size_t i = 0; i < test_values_abs.size(); i++) {
        T abs_value = PythonMath::abs(test_values_abs[i]);
        T abs_answer = std::abs(test_values_abs[i]);

        tester.expect_near(abs_value, abs_answer, NEAR_LIMIT_STRICT,
            "check PythonMath abs.");
    }

    std::vector<T> abs_value_vector = PythonMath::abs(test_values_abs);
    for (std::size_t i = 0; i < test_values_abs.size(); i++) {
        T abs_answer = std::abs(test_values_abs[i]);

        tester.expect_near(abs_value_vector[i], abs_answer, NEAR_LIMIT_STRICT,
            "check PythonMath abs vector.");
    }

    std::array<T, 11> test_values_abs_array = {
        test_values_abs[0],
        test_values_abs[1],
        test_values_abs[2],
        test_values_abs[3],
        test_values_abs[4],
        test_values_abs[5],
        test_values_abs[6],
        test_values_abs[7],
        test_values_abs[8],
        test_values_abs[9],
        test_values_abs[10]
    };

    std::array<T, 11> abs_value_array = PythonMath::abs(test_values_abs_array);
    for (std::size_t i = 0; i < test_values_abs.size(); i++) {
        T abs_answer = std::abs(test_values_abs[i]);

        tester.expect_near(abs_value_array[i], abs_answer, NEAR_LIMIT_STRICT,
            "check PythonMath abs array.");
    }

    std::vector<T> test_values_fmod({
        static_cast<T>(-10),
        static_cast<T>(-8.6),
        static_cast<T>(-5),
        static_cast<T>(-1),
        static_cast<T>(-0.1),
        static_cast<T>(0),
        static_cast<T>(0.1),
        static_cast<T>(1),
        static_cast<T>(5),
        static_cast<T>(8.7),
        static_cast<T>(10)
        });

    T divisor_for_fmod = static_cast<T>(3);

    /* fmod */
    for (std::size_t i = 0; i < test_values_fmod.size(); i++) {
        T fmod_value = PythonMath::fmod(test_values_fmod[i], divisor_for_fmod);
        T fmod_answer = std::fmod(test_values_fmod[i], divisor_for_fmod);

        tester.expect_near(fmod_value, fmod_answer, NEAR_LIMIT_STRICT,
            "check PythonMath fmod.");
    }

    std::vector<T> fmod_value_vector = PythonMath::fmod(test_values_fmod, divisor_for_fmod);
    for (std::size_t i = 0; i < test_values_fmod.size(); i++) {
        T fmod_answer = std::fmod(test_values_fmod[i], divisor_for_fmod);

        tester.expect_near(fmod_value_vector[i], fmod_answer, NEAR_LIMIT_STRICT,
            "check PythonMath fmod vector.");
    }

    std::array<T, 11> test_values_fmod_array = {
        test_values_fmod[0],
        test_values_fmod[1],
        test_values_fmod[2],
        test_values_fmod[3],
        test_values_fmod[4],
        test_values_fmod[5],
        test_values_fmod[6],
        test_values_fmod[7],
        test_values_fmod[8],
        test_values_fmod[9],
        test_values_fmod[10]
    };

    std::array<T, 11> fmod_value_array = PythonMath::fmod(test_values_fmod_array, divisor_for_fmod);
    for (std::size_t i = 0; i < test_values_fmod.size(); i++) {
        T fmod_answer = std::fmod(test_values_fmod[i], divisor_for_fmod);

        tester.expect_near(fmod_value_array[i], fmod_answer, NEAR_LIMIT_STRICT,
            "check PythonMath fmod array.");
    }


    tester.throw_error_if_test_failed();
}

template<typename T>
void check_python_math_exponential_logarithmic(void) {

    MCAPTester<T> tester;

#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    const T NEAR_LIMIT_SOFT = 1.0e-2F;
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    constexpr T NEAR_LIMIT_STRICT = std::is_same<T, double>::value ? T(1.0e-6) : T(1.0e-5);
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__


    std::vector<T> test_values_sqrt({
        static_cast<T>(-1),
        static_cast<T>(0),
        static_cast<T>(0.25),
        static_cast<T>(0.5),
        static_cast<T>(0.75),
        static_cast<T>(1),
        static_cast<T>(2),
        static_cast<T>(3),
        static_cast<T>(4),
        static_cast<T>(5),
        static_cast<T>(10)
        });

    /* sqrt */
    for (std::size_t i = 0; i < test_values_sqrt.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T sqrt_value = PythonMath::sqrt(test_values_sqrt[i]);
        T sqrt_answer = static_cast<T>(0);
        if (test_values_sqrt[i] < static_cast<T>(Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN)) {
            sqrt_answer = std::sqrt(static_cast<T>(Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN));
        }
        else {
            sqrt_answer = std::sqrt(test_values_sqrt[i]);
        }

        tester.expect_near(sqrt_value, sqrt_answer, NEAR_LIMIT_SOFT * std::abs(sqrt_answer),
            "check sqrt.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T sqrt_value = PythonMath::sqrt(test_values_sqrt[i]);
        T sqrt_answer = static_cast<T>(0);
        if (test_values_sqrt[i] < static_cast<T>(Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN)) {
            sqrt_answer = std::sqrt(static_cast<T>(Base::Math::EXPONENTIAL_LOGARITHMIC_DIVISION_MIN));
        }
        else {
            sqrt_answer = std::sqrt(test_values_sqrt[i]);
        }

        tester.expect_near(sqrt_value, sqrt_answer, NEAR_LIMIT_STRICT,
            "check sqrt.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

#ifndef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    std::vector<T> test_vector_sqrt({
        static_cast<T>(0.25),
        static_cast<T>(2),
        static_cast<T>(3)
    });

    std::vector<T> sqrt_value_vector = PythonMath::sqrt(test_vector_sqrt);
    for (std::size_t i = 0; i < test_vector_sqrt.size(); i++) {
        T sqrt_answer = std::sqrt(test_vector_sqrt[i]);

        tester.expect_near(sqrt_value_vector[i], sqrt_answer, NEAR_LIMIT_STRICT,
            "check sqrt vector.");
    }

    std::array<T, 3> test_array_sqrt = {
        test_vector_sqrt[0],
        test_vector_sqrt[1],
        test_vector_sqrt[2]
    };

    std::array<T, 3> sqrt_value_array = PythonMath::sqrt(test_array_sqrt);
    for (std::size_t i = 0; i < test_vector_sqrt.size(); i++) {
        T sqrt_answer = std::sqrt(test_vector_sqrt[i]);

        tester.expect_near(sqrt_value_array[i], sqrt_answer, NEAR_LIMIT_STRICT,
            "check sqrt array.");
    }
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__

    std::vector<T> test_values_exp({
        static_cast<T>(-100),
        static_cast<T>(-87),
        static_cast<T>(-50),
        static_cast<T>(-10),
        static_cast<T>(-1),
        static_cast<T>(0),
        static_cast<T>(1),
        static_cast<T>(2),
        static_cast<T>(50),
        static_cast<T>(87),
        static_cast<T>(100)
        });

    /* exp */
    for (std::size_t i = 0; i < test_values_exp.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T exp_value = PythonMath::exp(test_values_exp[i]);
        T exp_answer = static_cast<T>(0);
        if (test_values_exp[i] > static_cast<T>(Base::Math::EXP_INPUT_MAX)) {
            exp_answer = std::exp(static_cast<T>(Base::Math::EXP_INPUT_MAX));
        }
        else if (test_values_exp[i] < static_cast<T>(Base::Math::EXP_INPUT_MIN)) {
            exp_answer = std::exp(static_cast<T>(Base::Math::EXP_INPUT_MIN));
        }
        else {
            exp_answer = std::exp(test_values_exp[i]);
        }

        tester.expect_near(exp_value, exp_answer, (NEAR_LIMIT_SOFT * std::abs(exp_answer)),
            "check exp.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T exp_value = PythonMath::exp(test_values_exp[i]);
        T exp_answer = static_cast<T>(0);
        if (test_values_exp[i] > static_cast<T>(Base::Math::EXP_INPUT_MAX)) {
            exp_answer = std::exp(static_cast<T>(Base::Math::EXP_INPUT_MAX));
        }
        else if (test_values_exp[i] < static_cast<T>(Base::Math::EXP_INPUT_MIN)) {
            exp_answer = std::exp(static_cast<T>(Base::Math::EXP_INPUT_MIN));
        }
        else {
            exp_answer = std::exp(test_values_exp[i]);
        }

        tester.expect_near(exp_value, exp_answer, (NEAR_LIMIT_STRICT * std::abs(exp_answer)),
            "check exp.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

#ifndef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    std::vector<T> test_vector_exp({
        static_cast<T>(0),
        static_cast<T>(0.333),
        static_cast<T>(0.777)
        });

    std::vector<T> exp_value_vector = PythonMath::exp(test_vector_exp);
    for (std::size_t i = 0; i < test_vector_exp.size(); i++) {
        T exp_answer = std::exp(test_vector_exp[i]);

        tester.expect_near(exp_value_vector[i], exp_answer, NEAR_LIMIT_STRICT,
            "check exp vector.");
    }

    std::array<T, 3> test_array_exp = {
        test_vector_exp[0],
        test_vector_exp[1],
        test_vector_exp[2]
    };

    std::array<T, 3> exp_value_array = PythonMath::exp(test_array_exp);
    for (std::size_t i = 0; i < test_vector_exp.size(); i++) {
        T exp_answer = std::exp(test_vector_exp[i]);

        tester.expect_near(exp_value_array[i], exp_answer, NEAR_LIMIT_STRICT,
            "check exp array.");
    }
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__


    std::vector<T> test_values_exp2({
    static_cast<T>(-200),
    static_cast<T>(-126),
    static_cast<T>(-50),
    static_cast<T>(-8.445),
    static_cast<T>(-1),
    static_cast<T>(-0.777),
    static_cast<T>(0),
    static_cast<T>(0.777),
    static_cast<T>(1),
    static_cast<T>(8.445),
    static_cast<T>(50),
    static_cast<T>(126),
    static_cast<T>(200)
        });


    /* exp2 */
    for (std::size_t i = 0; i < test_values_exp2.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T exp2_value = PythonMath::exp2(test_values_exp2[i]);
        T exp2_answer = static_cast<T>(0);
        if (test_values_exp2[i] > static_cast<T>(Base::Math::EXP2_INPUT_MAX)) {
            exp2_answer = std::exp2(static_cast<T>(Base::Math::EXP2_INPUT_MAX));
        }
        else if (test_values_exp2[i] < static_cast<T>(Base::Math::EXP2_INPUT_MIN)) {
            exp2_answer = std::exp2(static_cast<T>(Base::Math::EXP2_INPUT_MIN));
        }
        else {
            exp2_answer = std::exp2(test_values_exp2[i]);
        }

        tester.expect_near(exp2_value, exp2_answer, (NEAR_LIMIT_SOFT * std::abs(exp2_answer)),
            "check exp2.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T exp2_value = PythonMath::exp2(test_values_exp2[i]);
        T exp2_answer = static_cast<T>(0);
        if (test_values_exp2[i] > static_cast<T>(Base::Math::EXP2_INPUT_MAX)) {
            exp2_answer = std::exp2(static_cast<T>(Base::Math::EXP2_INPUT_MAX));
        }
        else if (test_values_exp2[i] < static_cast<T>(Base::Math::EXP2_INPUT_MIN)) {
            exp2_answer = std::exp2(static_cast<T>(Base::Math::EXP2_INPUT_MIN));
        }
        else {
            exp2_answer = std::exp2(test_values_exp2[i]);
        }

        tester.expect_near(exp2_value, exp2_answer, (NEAR_LIMIT_STRICT* std::abs(exp2_answer)),
            "check exp2.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

#ifndef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    std::vector<T> test_vector_exp2({
        static_cast<T>(0),
        static_cast<T>(0.333),
        static_cast<T>(0.777)
        });

    std::vector<T> exp2_value_vector = PythonMath::exp2(test_vector_exp2);
    for (std::size_t i = 0; i < test_vector_exp2.size(); i++) {
        T exp2_answer = std::exp2(test_vector_exp2[i]);

        tester.expect_near(exp2_value_vector[i], exp2_answer, NEAR_LIMIT_STRICT,
            "check exp2 vector.");
    }

    std::array<T, 3> test_array_exp2 = {
        test_vector_exp2[0],
        test_vector_exp2[1],
        test_vector_exp2[2]
    };

    std::array<T, 3> exp2_value_array = PythonMath::exp2(test_array_exp2);
    for (std::size_t i = 0; i < test_vector_exp2.size(); i++) {
        T exp2_answer = std::exp2(test_vector_exp2[i]);

        tester.expect_near(exp2_value_array[i], exp2_answer, NEAR_LIMIT_STRICT,
            "check exp2 array.");
    }
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__

    std::vector<T> test_values_log({
        static_cast<T>(-1),
        static_cast<T>(0.001),
        static_cast<T>(0.1),
        static_cast<T>(0.5),
        static_cast<T>(0.77777),
        static_cast<T>(1),
        static_cast<T>(2),
        static_cast<T>(5),
        static_cast<T>(7.77777),
        static_cast<T>(10),
        static_cast<T>(50),
        static_cast<T>(100),
        static_cast<T>(500),
        static_cast<T>(1000),
        static_cast<T>(10000)
    });

    /* log */
    for (std::size_t i = 0; i < test_values_log.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T log_value = PythonMath::log(test_values_log[i]);
        T log_answer = static_cast<T>(0);
        if (test_values_log[i] <= static_cast<T>(0)) {
            log_value = static_cast<T>(Base::Math::LOG_OUTPUT_MIN);
            log_answer = static_cast<T>(Base::Math::LOG_OUTPUT_MIN);
        }
        else {
            log_answer = std::log(test_values_log[i]);
        }

        tester.expect_near(log_value, log_answer, NEAR_LIMIT_SOFT * std::abs(log_answer),
            "check log.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T log_value = PythonMath::log(test_values_log[i]);
        T log_answer = static_cast<T>(0);
        if (test_values_log[i] <= static_cast<T>(0)) {
            log_answer = static_cast<T>(Base::Math::LOG_OUTPUT_MIN);
        }
        else {
            log_answer = std::log(test_values_log[i]);
        }

        tester.expect_near(log_value, log_answer, NEAR_LIMIT_STRICT * std::abs(log_answer),
            "check log.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

#ifndef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    std::vector<T> test_vector_log({
        static_cast<T>(0.33),
        static_cast<T>(1.5),
        static_cast<T>(2.5)
        });

    std::vector<T> log_value_vector = PythonMath::log(test_vector_log);
    for (std::size_t i = 0; i < test_vector_log.size(); i++) {
        T log_answer = std::log(test_vector_log[i]);

        tester.expect_near(log_value_vector[i], log_answer, NEAR_LIMIT_STRICT,
            "check log vector.");
    }

    std::array<T, 3> test_array_log = {
        test_vector_log[0],
        test_vector_log[1],
        test_vector_log[2]
    };

    std::array<T, 3> log_value_array = PythonMath::log(test_array_log);
    for (std::size_t i = 0; i < test_vector_log.size(); i++) {
        T log_answer = std::log(test_vector_log[i]);

        tester.expect_near(log_value_array[i], log_answer, NEAR_LIMIT_STRICT,
            "check log array.");
    }
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__

    /* log2 */
    for (std::size_t i = 0; i < test_values_log.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T log2_value = PythonMath::log2(test_values_log[i]);
        T log2_answer = static_cast<T>(0);
        if (test_values_log[i] <= static_cast<T>(0)) {
            log2_value = static_cast<T>(Base::Math::LOG_OUTPUT_MIN)
                / static_cast<T>(Base::Math::LN_2);
            log2_answer = static_cast<T>(Base::Math::LOG_OUTPUT_MIN)
                / static_cast<T>(Base::Math::LN_2);
        }
        else {
            log2_answer = std::log2(test_values_log[i]);
        }

        tester.expect_near(log2_value, log2_answer, NEAR_LIMIT_SOFT * std::abs(log2_answer),
            "check log2.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T log2_value = PythonMath::log2(test_values_log[i]);
        T log2_answer = static_cast<T>(0);
        if (test_values_log[i] <= static_cast<T>(0)) {
            log2_answer = static_cast<T>(Base::Math::LOG_OUTPUT_MIN)
                / static_cast<T>(Base::Math::LN_2);
        }
        else {
            log2_answer = std::log2(test_values_log[i]);
        }

        tester.expect_near(log2_value, log2_answer, NEAR_LIMIT_STRICT * std::abs(log2_answer),
            "check log2.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

#ifndef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    std::vector<T> test_vector_log2({
        static_cast<T>(0.33),
        static_cast<T>(1.5),
        static_cast<T>(2.5)
        });

    std::vector<T> log2_value_vector = PythonMath::log2(test_vector_log2);
    for (std::size_t i = 0; i < test_vector_log2.size(); i++) {
        T log2_answer = std::log2(test_vector_log2[i]);

        tester.expect_near(log2_value_vector[i], log2_answer, NEAR_LIMIT_STRICT,
            "check log2 vector.");
    }

    std::array<T, 3> test_array_log2 = {
        test_vector_log2[0],
        test_vector_log2[1],
        test_vector_log2[2]
    };

    std::array<T, 3> log2_value_array = PythonMath::log2(test_array_log2);
    for (std::size_t i = 0; i < test_vector_log2.size(); i++) {
        T log2_answer = std::log2(test_vector_log2[i]);

        tester.expect_near(log2_value_array[i], log2_answer, NEAR_LIMIT_STRICT,
            "check log2 array.");
    }
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__

    /* log10 */
    for (std::size_t i = 0; i < test_values_log.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T log10_value = PythonMath::log10(test_values_log[i]);
        T log10_answer = static_cast<T>(0);
        if (test_values_log[i] <= static_cast<T>(0)) {
            log10_value = static_cast<T>(Base::Math::LOG_OUTPUT_MIN)
                / static_cast<T>(Base::Math::LN_10);
            log10_answer = static_cast<T>(Base::Math::LOG_OUTPUT_MIN)
                / static_cast<T>(Base::Math::LN_10);
        }
        else {
            log10_answer = std::log10(test_values_log[i]);
        }

        tester.expect_near(log10_value, log10_answer, NEAR_LIMIT_SOFT * std::abs(log10_answer),
            "check log10.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T log10_value = PythonMath::log10(test_values_log[i]);
        T log10_answer = static_cast<T>(0);
        if (test_values_log[i] <= static_cast<T>(0)) {
            log10_answer = static_cast<T>(Base::Math::LOG_OUTPUT_MIN)
                / static_cast<T>(Base::Math::LN_10);
        }
        else {
            log10_answer = std::log10(test_values_log[i]);
        }

        tester.expect_near(log10_value, log10_answer, NEAR_LIMIT_STRICT * std::abs(log10_answer),
            "check log10.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

#ifndef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    std::vector<T> test_vector_log10({
        static_cast<T>(0.33),
        static_cast<T>(1.5),
        static_cast<T>(2.5)
        });

    std::vector<T> log10_value_vector = PythonMath::log10(test_vector_log10);
    for (std::size_t i = 0; i < test_vector_log10.size(); i++) {
        T log10_answer = std::log10(test_vector_log10[i]);

        tester.expect_near(log10_value_vector[i], log10_answer, NEAR_LIMIT_STRICT,
            "check log10 vector.");
    }

    std::array<T, 3> test_array_log10 = {
        test_vector_log10[0],
        test_vector_log10[1],
        test_vector_log10[2]
    };

    std::array<T, 3> log10_value_array = PythonMath::log10(test_array_log10);
    for (std::size_t i = 0; i < test_vector_log10.size(); i++) {
        T log10_answer = std::log10(test_vector_log10[i]);

        tester.expect_near(log10_value_array[i], log10_answer, NEAR_LIMIT_STRICT,
            "check log10 array.");
    }
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__

    std::vector<std::pair<T, T>> test_values_pow({
        {static_cast<T>(0.1), static_cast<T>(0.1)},
        {static_cast<T>(0.1), static_cast<T>(0.5)},
        {static_cast<T>(0.1), static_cast<T>(1)},
        {static_cast<T>(0.1), static_cast<T>(2)},
        {static_cast<T>(0.5), static_cast<T>(0.1)},
        {static_cast<T>(0.5), static_cast<T>(0.5)},
        {static_cast<T>(0.5), static_cast<T>(1)},
        {static_cast<T>(0.5), static_cast<T>(2)},
        {static_cast<T>(1), static_cast<T>(0.1)},
        {static_cast<T>(1), static_cast<T>(0.5)},
        {static_cast<T>(1), static_cast<T>(1)},
        {static_cast<T>(1), static_cast<T>(2)},
        {static_cast<T>(2), static_cast<T>(0.1)},
        {static_cast<T>(2), static_cast<T>(0.5)},
        {static_cast<T>(2), static_cast<T>(1)},
        {static_cast<T>(2), static_cast<T>(2)},
        {static_cast<T>(5), static_cast<T>(0.1)},
        {static_cast<T>(5), static_cast<T>(0.5)},
        {static_cast<T>(5), static_cast<T>(1)},
        {static_cast<T>(5), static_cast<T>(2)},
        {static_cast<T>(10), static_cast<T>(0.1)},
        {static_cast<T>(10), static_cast<T>(0.5)},
        {static_cast<T>(10), static_cast<T>(1)},
        {static_cast<T>(10), static_cast<T>(2)},
        {static_cast<T>(50), static_cast<T>(0.1)},
        {static_cast<T>(50), static_cast<T>(0.5)},
        {static_cast<T>(50), static_cast<T>(1)},
        {static_cast<T>(50), static_cast<T>(2)},
        {static_cast<T>(100), static_cast<T>(0.1)},
        {static_cast<T>(100), static_cast<T>(0.5)},
        {static_cast<T>(100), static_cast<T>(1)},
        {static_cast<T>(100), static_cast<T>(2)},
        });

    /* pow */
    for (std::size_t i = 0; i < test_values_pow.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T pow_value = PythonMath::pow(test_values_pow[i].first, test_values_pow[i].second);
        T pow_answer = std::pow(test_values_pow[i].first, test_values_pow[i].second);

        tester.expect_near(pow_value, pow_answer, NEAR_LIMIT_SOFT * std::abs(pow_answer),
            "check pow.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T pow_value = PythonMath::pow(test_values_pow[i].first, test_values_pow[i].second);
        T pow_answer = std::pow(test_values_pow[i].first, test_values_pow[i].second);

        tester.expect_near(pow_value, pow_answer, NEAR_LIMIT_STRICT * std::abs(pow_answer),
            "check pow.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

#ifndef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    std::vector<T> test_vector_pow_x({
        static_cast<T>(0.5),
        static_cast<T>(0.5),
        static_cast<T>(1), 
        });
    std::vector<T> test_vector_pow_y({
        static_cast<T>(0.5),
        static_cast<T>(2),
        static_cast<T>(0.5),
        });

    std::vector<T> pow_value_vector = PythonMath::pow(test_vector_pow_x, test_vector_pow_y[0]);
    for (std::size_t i = 0; i < test_vector_pow_x.size(); i++) {
        T pow_answer = std::pow(test_vector_pow_x[i], test_vector_pow_y[0]);

        tester.expect_near(pow_value_vector[i], pow_answer, NEAR_LIMIT_STRICT,
            "check pow vector and element.");
    }

    pow_value_vector = PythonMath::pow(test_vector_pow_x[0], test_vector_pow_y);
    for (std::size_t i = 0; i < test_vector_pow_y.size(); i++) {
        T pow_answer = std::pow(test_vector_pow_x[0], test_vector_pow_y[i]);

        tester.expect_near(pow_value_vector[i], pow_answer, NEAR_LIMIT_STRICT,
            "check pow element and vector.");
    }

    pow_value_vector = PythonMath::pow(test_vector_pow_x, test_vector_pow_y);
    for (std::size_t i = 0; i < test_vector_pow_x.size(); i++) {
        T pow_answer = std::pow(test_vector_pow_x[i], test_vector_pow_y[i]);

        tester.expect_near(pow_value_vector[i], pow_answer, NEAR_LIMIT_STRICT,
            "check pow vector and vector.");
    }

    std::array<T, 3> test_array_pow_x = {
        test_vector_pow_x[0],
        test_vector_pow_x[1],
        test_vector_pow_x[2]
    };
    std::array<T, 3> test_array_pow_y = {
        test_vector_pow_y[0],
        test_vector_pow_y[1],
        test_vector_pow_y[2]
    };

    std::array<T, 3> pow_value_array = PythonMath::pow(test_array_pow_x, test_array_pow_y[0]);
    for (std::size_t i = 0; i < test_vector_pow_x.size(); i++) {
        T pow_answer = std::pow(test_vector_pow_x[i], test_array_pow_y[0]);

        tester.expect_near(pow_value_array[i], pow_answer, NEAR_LIMIT_STRICT,
            "check pow array and element.");
    }

    pow_value_array = PythonMath::pow(test_array_pow_x[0], test_array_pow_y);
    for (std::size_t i = 0; i < test_array_pow_y.size(); i++) {
        T pow_answer = std::pow(test_array_pow_x[0], test_array_pow_y[i]);

        tester.expect_near(pow_value_array[i], pow_answer, NEAR_LIMIT_STRICT,
            "check pow element and array.");
    }

    pow_value_array = PythonMath::pow(test_array_pow_x, test_array_pow_y);
    for (std::size_t i = 0; i < test_array_pow_x.size(); i++) {
        T pow_answer = std::pow(test_array_pow_x[i], test_array_pow_y[i]);

        tester.expect_near(pow_value_array[i], pow_answer, NEAR_LIMIT_STRICT,
            "check pow array and array.");
    }
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__


    tester.throw_error_if_test_failed();
}

template<typename T>
void check_python_math_trigonometric(void) {

    MCAPTester<T> tester;

#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    const T NEAR_LIMIT_SOFT = 1.0e-2F;
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    constexpr T NEAR_LIMIT_STRICT = std::is_same<T, double>::value ? T(1.0e-6) : T(1.0e-5);
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__


    std::vector<T> test_values_sin = {
        static_cast<T>(-10 * Base::Math::PI),
        static_cast<T>(-3 * Base::Math::PI),
        static_cast<T>(-2 * Base::Math::PI),
        static_cast<T>(-3 * Base::Math::PI / 2),
        static_cast<T>(-Base::Math::PI),
        static_cast<T>(-Base::Math::PI / 2),
        static_cast<T>(0),
        static_cast<T>(Base::Math::PI / 2),
        static_cast<T>(Base::Math::PI),
        static_cast<T>(3 * Base::Math::PI / 2),
        static_cast<T>(2 * Base::Math::PI),
        static_cast<T>(-5 * Base::Math::PI / 4),
        static_cast<T>(5 * Base::Math::PI / 4),
        static_cast<T>(-7 * Base::Math::PI / 4),
        static_cast<T>(7 * Base::Math::PI / 4),
        static_cast<T>(3 * Base::Math::PI),
        static_cast<T>(10 * Base::Math::PI)
    };

    /* sin */
    for (std::size_t i = 0; i < test_values_sin.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T sin_value = PythonMath::sin(test_values_sin[i]);
        T sin_answer = std::sin(test_values_sin[i]);

        tester.expect_near(sin_value, sin_answer, NEAR_LIMIT_SOFT,
            "check sin.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T sin_value = PythonMath::sin(test_values_sin[i]);
        T sin_answer = std::sin(test_values_sin[i]);

        tester.expect_near(sin_value, sin_answer, NEAR_LIMIT_STRICT,
            "check sin.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

    /* cos */
    for (std::size_t i = 0; i < test_values_sin.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T cos_value = PythonMath::cos(test_values_sin[i]);
        T cos_answer = std::cos(test_values_sin[i]);

        tester.expect_near(cos_value, cos_answer, NEAR_LIMIT_SOFT,
            "check cos.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T cos_value = PythonMath::cos(test_values_sin[i]);
        T cos_answer = std::cos(test_values_sin[i]);

        tester.expect_near(cos_value, cos_answer, NEAR_LIMIT_STRICT,
            "check cos.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

    /* tan */
    for (std::size_t i = 0; i < test_values_sin.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T tan_value = PythonMath::tan(test_values_sin[i]);
        T tan_answer = std::tan(test_values_sin[i]);

        if (std::abs(tan_value) > static_cast<T>(1 / NEAR_LIMIT_SOFT) &&
            std::abs(tan_answer) > static_cast<T>(1 / NEAR_LIMIT_SOFT)) {
            /* Do Nothing. */
        }
        else {
            tester.expect_near(tan_value, tan_answer, NEAR_LIMIT_SOFT,
                "check tan.");
        }
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T tan_value = PythonMath::tan(test_values_sin[i]);
        T tan_answer = std::tan(test_values_sin[i]);

        if (std::abs(tan_value) > static_cast<T>(1 / NEAR_LIMIT_STRICT) &&
            std::abs(tan_answer) > static_cast<T>(1 / NEAR_LIMIT_STRICT)) {
            /* Do Nothing. */
        }
        else {
            tester.expect_near(tan_value, tan_answer, NEAR_LIMIT_STRICT,
                "check tan.");
        }
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }


    std::vector<T> test_values_atan = {
    static_cast<T>(-10),
    static_cast<T>(-3),
    static_cast<T>(-1.1),
    static_cast<T>(-1),
    static_cast<T>(-0.9),
    static_cast<T>(-0.5),
    static_cast<T>(0),
    static_cast<T>(0.5),
    static_cast<T>(0.9),
    static_cast<T>(1),
    static_cast<T>(1.1),
    static_cast<T>(3),
    static_cast<T>(10)
    };

    /* atan */
    for (std::size_t i = 0; i < test_values_atan.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T atan_value = PythonMath::atan(test_values_atan[i]);
        T atan_answer = std::atan(test_values_atan[i]);

        tester.expect_near(atan_value, atan_answer, NEAR_LIMIT_SOFT,
            "check atan.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T atan_value = PythonMath::atan(test_values_atan[i]);
        T atan_answer = std::atan(test_values_atan[i]);

        tester.expect_near(atan_value, atan_answer, NEAR_LIMIT_STRICT,
            "check atan.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

    std::vector<std::pair<T, T>> test_values_atan2({
            {static_cast<T>(0), static_cast<T>(0)},    /* 原点 */
            {static_cast<T>(1), static_cast<T>(0)},    /* 正のx軸 */
            {static_cast<T>(0), static_cast<T>(1)},    /* 正のy軸 */
            {static_cast<T>(-1), static_cast<T>(0)},   /* 負のx軸 */
            {static_cast<T>(0), static_cast<T>(-1)},   /* 負のy軸 */
            {static_cast<T>(1), static_cast<T>(1)},    /* 第1象限 */
            {static_cast<T>(-1), static_cast<T>(1)},   /* 第2象限 */
            {static_cast<T>(-1), static_cast<T>(-1)},  /* 第3象限 */
            {static_cast<T>(1), static_cast<T>(-1)},   /* 第4象限 */
            {static_cast<T>(3), static_cast<T>(4)},    /* 一般的な正のケース */
            {static_cast<T>(-3), static_cast<T>(4)},   /* 第2象限 */
            {static_cast<T>(-3), static_cast<T>(-4)},  /* 第3象限 */
            {static_cast<T>(3), static_cast<T>(-4)},   /* 第4象限 */
            {static_cast<T>(0.5), static_cast<T>(0.5)}, /* 小さい正の値 */
            {static_cast<T>(-0.5), static_cast<T>(0.5)}, /* 小さい負の値 */
            {static_cast<T>(1e-6), static_cast<T>(1e-6)}, /* 非常に小さい正の値 */
            {static_cast<T>(-1e-6), static_cast<T>(1e-6)}, /* 非常に小さい負の値 */
            {static_cast<T>(1), static_cast<T>(1e6)},  /* 大きな差があるケース */
            {static_cast<T>(1e6), static_cast<T>(1)},  /* 大きな差が逆になるケース */
            {static_cast<T>(-1), static_cast<T>(1e6)}, /* 異なる象限のケース */
        });

    /* atan2 */
    for (std::size_t i = 0; i < test_values_atan2.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T atan2_value = PythonMath::atan2(test_values_atan2[i].second, test_values_atan2[i].first);
        T atan2_answer = std::atan2(test_values_atan2[i].second, test_values_atan2[i].first);

        tester.expect_near(atan2_value, atan2_answer, NEAR_LIMIT_SOFT,
            "check atan2.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T atan2_value = PythonMath::atan2(test_values_atan2[i].second, test_values_atan2[i].first);
        T atan2_answer = std::atan2(test_values_atan2[i].second, test_values_atan2[i].first);

        tester.expect_near(atan2_value, atan2_answer, NEAR_LIMIT_STRICT,
            "check atan2.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

#ifndef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__

    std::vector<T> test_values_atan2_y({
        static_cast<T>(-3),
        static_cast<T>(3),
        static_cast<T>(0.5)
        });
    std::vector<T> test_values_atan2_x({
        static_cast<T>(-4),
        static_cast<T>(-4),
        static_cast<T>(0.5)
        });

    std::vector<T> atan2_value_vector = PythonMath::atan2(test_values_atan2_y, test_values_atan2_x[0]);
    for (std::size_t i = 0; i < test_values_atan2_y.size(); i++) {
        T atan2_answer = std::atan2(test_values_atan2_y[i], test_values_atan2_x[0]);

        tester.expect_near(atan2_value_vector[i], atan2_answer, NEAR_LIMIT_STRICT,
            "check atan2 vector and element.");
    }

    atan2_value_vector = PythonMath::atan2(test_values_atan2_y[0], test_values_atan2_x);
    for (std::size_t i = 0; i < test_values_atan2_x.size(); i++) {
        T atan2_answer = std::atan2(test_values_atan2_y[0], test_values_atan2_x[i]);

        tester.expect_near(atan2_value_vector[i], atan2_answer, NEAR_LIMIT_STRICT,
            "check atan2 element and vector.");
    }

    atan2_value_vector = PythonMath::atan2(test_values_atan2_y, test_values_atan2_x);
    for (std::size_t i = 0; i < test_values_atan2_y.size(); i++) {
        T atan2_answer = std::atan2(test_values_atan2_y[i], test_values_atan2_x[i]);

        tester.expect_near(atan2_value_vector[i], atan2_answer, NEAR_LIMIT_STRICT,
            "check atan2 vector and vector.");
    }

    std::array<T, 3> test_array_atan2_y = {
        test_values_atan2_y[0],
        test_values_atan2_y[1],
        test_values_atan2_y[2]
    };
    std::array<T, 3> test_array_atan2_x = {
        test_values_atan2_x[0],
        test_values_atan2_x[1],
        test_values_atan2_x[2]
    };

    std::array<T, 3> atan2_value_array = PythonMath::atan2(test_array_atan2_y, test_array_atan2_x[0]);
    for (std::size_t i = 0; i < test_values_atan2_y.size(); i++) {
        T atan2_answer = std::atan2(test_values_atan2_y[i], test_array_atan2_x[0]);

        tester.expect_near(atan2_value_array[i], atan2_answer, NEAR_LIMIT_STRICT,
            "check atan2 array and element.");
    }

    atan2_value_array = PythonMath::atan2(test_array_atan2_y[0], test_array_atan2_x);
    for (std::size_t i = 0; i < test_array_atan2_x.size(); i++) {
        T atan2_answer = std::atan2(test_array_atan2_y[0], test_array_atan2_x[i]);

        tester.expect_near(atan2_value_array[i], atan2_answer, NEAR_LIMIT_STRICT,
            "check atan2 element and array.");
    }

    atan2_value_array = PythonMath::atan2(test_array_atan2_y, test_array_atan2_x);
    for (std::size_t i = 0; i < test_values_atan2_y.size(); i++) {
        T atan2_answer = std::atan2(test_values_atan2_y[i], test_values_atan2_x[i]);

        tester.expect_near(atan2_value_array[i], atan2_answer, NEAR_LIMIT_STRICT,
            "check atan2 array and array.");
    }
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__

    std::vector<T> test_values_asin_acos = {
        static_cast<T>(-2),
        static_cast<T>(-1),
        static_cast<T>(-0.86602540378), // -sqrt(3)/2
        static_cast<T>(-0.70710678118), // -sqrt(2)/2
        static_cast<T>(-0.5),
        static_cast<T>(-0.25881904510), // -1/2 * sqrt(3)/2
        static_cast<T>(0),
        static_cast<T>(0.25881904510),  // 1/2 * sqrt(3)/2
        static_cast<T>(0.5),
        static_cast<T>(0.70710678118),  // sqrt(2)/2
        static_cast<T>(0.86602540378),  // sqrt(3)/2
        static_cast<T>(1),
        static_cast<T>(2)
    };

    /* asin */
    for (std::size_t i = 0; i < test_values_asin_acos.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T asin_value = PythonMath::asin(test_values_asin_acos[i]);

        T asin_answer = static_cast<T>(0);
        if (test_values_asin_acos[i] < static_cast<T>(-1)) {
            asin_answer = std::asin(static_cast<T>(-1));
        }
        else if (test_values_asin_acos[i] > static_cast<T>(1)) {
            asin_answer = std::asin(static_cast<T>(1));
        }
        else
        {
            asin_answer = std::asin(test_values_asin_acos[i]);
        }

        tester.expect_near(asin_value, asin_answer, NEAR_LIMIT_SOFT,
            "check asin.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T asin_value = PythonMath::asin(test_values_asin_acos[i]);

        T asin_answer = static_cast<T>(0);
        if (test_values_asin_acos[i] < static_cast<T>(-1)) {
            asin_answer = std::asin(static_cast<T>(-1));
        }
        else if (test_values_asin_acos[i] > static_cast<T>(1)) {
            asin_answer = std::asin(static_cast<T>(1));
        }
        else
        {
            asin_answer = std::asin(test_values_asin_acos[i]);
        }

        tester.expect_near(asin_value, asin_answer, NEAR_LIMIT_STRICT,
            "check asin.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

    /* acos */
    for (std::size_t i = 0; i < test_values_asin_acos.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T acos_value = PythonMath::acos(test_values_asin_acos[i]);

        T acos_answer = static_cast<T>(0);
        if (test_values_asin_acos[i] < static_cast<T>(-1)) {
            acos_answer = std::acos(static_cast<T>(-1));
        }
        else if (test_values_asin_acos[i] > static_cast<T>(1)) {
            acos_answer = std::acos(static_cast<T>(1));
        }
        else
        {
            acos_answer = std::acos(test_values_asin_acos[i]);
        }

        tester.expect_near(acos_value, acos_answer, NEAR_LIMIT_SOFT,
            "check acos.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T acos_value = PythonMath::acos(test_values_asin_acos[i]);

        T acos_answer = static_cast<T>(0);
        if (test_values_asin_acos[i] < static_cast<T>(-1)) {
            acos_answer = std::acos(static_cast<T>(-1));
        }
        else if (test_values_asin_acos[i] > static_cast<T>(1)) {
            acos_answer = std::acos(static_cast<T>(1));
        }
        else
        {
            acos_answer = std::acos(test_values_asin_acos[i]);
        }

        tester.expect_near(acos_value, acos_answer, NEAR_LIMIT_STRICT,
            "check acos.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

    std::vector<T> test_values_hyperbolic({
        static_cast<T>(-10),
        static_cast<T>(-8.6),
        static_cast<T>(-5),
        static_cast<T>(-1),
        static_cast<T>(-0.1),
        static_cast<T>(0),
        static_cast<T>(0.1),
        static_cast<T>(1),
        static_cast<T>(5),
        static_cast<T>(8.7),
        static_cast<T>(10)
        });

    /* sinh */
    for (std::size_t i = 0; i < test_values_hyperbolic.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T sinh_value = PythonMath::sinh(test_values_hyperbolic[i]);
        T sinh_answer = std::sinh(test_values_hyperbolic[i]);

        tester.expect_near(sinh_value, sinh_answer, NEAR_LIMIT_SOFT * std::abs(sinh_answer),
            "check sinh.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T sinh_value = PythonMath::sinh(test_values_hyperbolic[i]);
        T sinh_answer = std::sinh(test_values_hyperbolic[i]);

        tester.expect_near(sinh_value, sinh_answer, NEAR_LIMIT_STRICT * std::abs(sinh_answer),
            "check sinh.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

    /* cosh */
    for (std::size_t i = 0; i < test_values_hyperbolic.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T cosh_value = PythonMath::cosh(test_values_hyperbolic[i]);
        T cosh_answer = std::cosh(test_values_hyperbolic[i]);

        tester.expect_near(cosh_value, cosh_answer, NEAR_LIMIT_SOFT * std::abs(cosh_answer),
            "check cosh.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T cosh_value = PythonMath::cosh(test_values_hyperbolic[i]);
        T cosh_answer = std::cosh(test_values_hyperbolic[i]);

        tester.expect_near(cosh_value, cosh_answer, NEAR_LIMIT_STRICT * std::abs(cosh_answer),
            "check cosh.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }

    /* tanh */
    for (std::size_t i = 0; i < test_values_hyperbolic.size(); i++) {
#ifdef __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T tanh_value = PythonMath::tanh(test_values_hyperbolic[i]);
        T tanh_answer = std::tanh(test_values_hyperbolic[i]);

        tester.expect_near(tanh_value, tanh_answer, NEAR_LIMIT_SOFT,
            "check tanh.");
#else // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
        T tanh_value = PythonMath::tanh(test_values_hyperbolic[i]);
        T tanh_answer = std::tanh(test_values_hyperbolic[i]);

        tester.expect_near(tanh_value, tanh_answer, NEAR_LIMIT_STRICT,
            "check tanh.");
#endif // __BASE_MATH_USE_ROUGH_BUT_FAST_APPROXIMATIONS__
    }


    tester.throw_error_if_test_failed();
}



template<typename T>
void check_python_math_calc(void) {

    check_python_math_arithmetic<T>();

    check_python_math_exponential_logarithmic<T>();

    check_python_math_trigonometric<T>();
}

int main() {

    check_base_math_calc<double>();

    check_base_math_calc<float>();

    check_python_math_calc<double>();

    check_python_math_calc<float>();
 

    return 0;
}
