#include "storage/fieldvariable.h"
#include <gtest/gtest.h>
#include <math.h>

TEST(FieldVariableTest, ValueCheck) {
    std::array<int, 2> size = std::array < int,2 > {3, 3};
    std::array<double, 2> origin = std::array < double,2 > {0.0, 0.0};
    std::array<double, 2> meshWidth = std::array < double,2 > {1.0, 1.0};
    FieldVariable fv = FieldVariable(size, origin, meshWidth);

    fv(1, 1) = 1.0;
    fv(2, 1) = 2.0;
    ASSERT_EQ(0.0, fv.interpolateAt(0.0, 0.0));
    ASSERT_EQ(1.0, fv.interpolateAt(1.0, 1.0));
    ASSERT_EQ(0.25, fv.interpolateAt(0.5, 0.5));
    ASSERT_EQ(0.5, fv.interpolateAt(1.0, 0.5));
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}