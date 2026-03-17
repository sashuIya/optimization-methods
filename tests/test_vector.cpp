#include <gtest/gtest.h>

#include "Vector.h"

TEST(VectorTest, BasicOps) {
  Vector v1({1.0, 2.0, 3.0});
  Vector v2({4.0, 5.0, 6.0});

  Vector v3 = v1 + v2;
  EXPECT_DOUBLE_EQ(v3[0], 5.0);
  EXPECT_DOUBLE_EQ(v3[1], 7.0);
  EXPECT_DOUBLE_EQ(v3[2], 9.0);

  Vector v4 = v1 * 2.0;
  EXPECT_DOUBLE_EQ(v4[0], 2.0);
  EXPECT_DOUBLE_EQ(v4[1], 4.0);
  EXPECT_DOUBLE_EQ(v4[2], 6.0);
}

TEST(VectorTest, Norm) {
  Vector v({3.0, 4.0});
  EXPECT_DOUBLE_EQ(v.Norm(), 5.0);
  EXPECT_DOUBLE_EQ(v.NormSq(), 25.0);
}
