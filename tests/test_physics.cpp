#include <gtest/gtest.h>
#include "LennardJones.hpp"
#include <cmath>

TEST(PhysicsTest, BasicEnergy) {
    LennardJonesSystem system(2, 1);
    // r = 1.0 => dist = 1.0
    // U = 1^-12 - 1^-6 = 0
    Vector x({0.0, 1.0});
    EXPECT_NEAR(system.compute_total_energy(x), 0.0, 1e-10);
}

TEST(PhysicsTest, Gradient) {
    LennardJonesSystem system(2, 1);
    Vector x({1.0, 2.0});
    Vector grad = system.compute_gradient(x);
    // dU/dx1 = -2(x1-x2) * (6 - 3d^3) d^-7 / n
    // d = 1.0, n = 2
    // grad1 = -2(1-2) * (6-3) * 1 / 2 = -2(-1) * 3 / 2 = 3.0
    EXPECT_NEAR(grad[0], 3.0, 1e-10);
    EXPECT_NEAR(grad[1], -3.0, 1e-10);
}
