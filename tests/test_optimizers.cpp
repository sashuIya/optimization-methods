#include <gtest/gtest.h>

#include <cmath>

#include "Optimizers.h"

TEST(OptimizerTest, GoldenSectionSearch) {
  auto f = [](double x) { return (x - 2.0) * (x - 2.0); };
  double res = LineSearch::GoldMethod(f, 0.0, 5.0);
  EXPECT_NEAR(res, 2.0, 1e-6);
}

TEST(OptimizerTest, LineSearchBestLambda) {
  auto f = [](double x) { return (x - 0.5) * (x - 0.5); };
  double res = LineSearch::GetBestLambda(f);
  EXPECT_NEAR(res, 0.5, 1e-4);
}
