#include <gtest/gtest.h>
#include "storage/fieldvariable.h"
#include "discretization/2_central_differences.h"
#include <cmath>

TEST(FieldVariableTest, ValueCheck) {
    auto size = std::array<int, 2>{3, 3};
    auto origin = std::array<double, 2>{0.0, 0.0};
    auto meshWidth = std::array<double, 2>{1.0, 1.0};
    FieldVariable fv = FieldVariable(size, origin, meshWidth);
    auto nCells = std::array<int, 2>{3, 3};
    auto discretization = new CentralDifferences(nCells, meshWidth);

    fv(1, 1) = 1.0;
    fv(2, 1) = 2.0;
    ASSERT_EQ(0.0, fv.interpolateAt(0.0, 0.0));
    ASSERT_EQ(1.0, fv.interpolateAt(1.0, 1.0));
    ASSERT_EQ(0.25, fv.interpolateAt(0.5, 0.5));
    ASSERT_EQ(0.5, fv.interpolateAt(1.0, 0.5));
}

TEST(SORTest, Test1) {
    auto nCells = std::array<int, 2>{3, 3};
    auto meshWidth = std::array<double, 2>{1.0, 1.0};
    auto d = new CentralDifferences(nCells, meshWidth);
    auto p = d->p();
    for (int i = d->pIBegin(); i <= d->pIEnd(); i++) {
        for (int j = d->pJBegin(); j <= d->pJEnd(); j++) {
            d->p(i, j) = 10 * i + j;
        }
    }
    for (int i = d->rhsIBegin(); i <= d->rhsIEnd(); i++) {
        for (int j = d->rhsJBegin(); j <= d->rhsJEnd(); j++) {
            std::cout << "iter(" << i << "," << j << ")" << std::endl;
            d->rhs(i, j) = (p(i+1,j)-2*p(i,j)+p(i-1,j))/pow(d->dx(),2); +
                            (p(i,j+1)-2*p(i,j)+p(i,j-1))/pow(d->dy(),2);
        }
    }

    /*ASSERT_EQ(0.0, fv.interpolateAt(0.0, 0.0));
    ASSERT_EQ(1.0, fv.interpolateAt(1.0, 1.0));
    ASSERT_EQ(0.25, fv.interpolateAt(0.5, 0.5));
    ASSERT_EQ(0.5, fv.interpolateAt(1.0, 0.5));*/
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}