#include <gtest/gtest.h>

#include <cmath>

#include "LennardJones.h"

TEST(PhysicsTest, BasicEnergy) {
  LennardJonesSystem system(2, 1);
  Vector x({0.0, 1.0});
  EXPECT_NEAR(system.ComputeTotalEnergy(x), 0.0, 1e-10);
}

TEST(PhysicsTest, Gradient) {
  LennardJonesSystem system(2, 1);
  Vector x({1.0, 2.0});
  Vector grad = system.ComputeGradient(x);
  EXPECT_NEAR(grad[0], 3.0, 1e-10);
  EXPECT_NEAR(grad[1], -3.0, 1e-10);
}
