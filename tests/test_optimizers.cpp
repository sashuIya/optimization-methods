#include <gtest/gtest.h>
#include "Optimizers.h"
#include <cmath>

TEST(OptimizerTest, GoldenSectionSearch) {
    auto f = [](double x) { return (x - 2.0) * (x - 2.0); };
    double res = LineSearch::gold_method(f, 0.0, 5.0);
    EXPECT_NEAR(res, 2.0, 1e-6);
}

TEST(OptimizerTest, LineSearchBestLambda) {
    auto f = [](double x) { return (x - 0.5) * (x - 0.5); };
    double res = LineSearch::get_best_lambda(f);
    EXPECT_NEAR(res, 0.5, 1e-4);
}
