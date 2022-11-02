#include <gtest/gtest.h>
#include "storage/fieldvariable.h"
#include "discretization/2_central_differences.h"
#include "pressure_solver/sor.h"
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
    auto origin = std::array<double, 2>{meshWidth[0] / 2.0, meshWidth[1] / 2.0};
    auto d = new CentralDifferences(nCells, meshWidth);
    auto pRef = FieldVariable({nCells[0] + 1, nCells[1] + 1}, origin, meshWidth);
    auto rhs = d->rhs();
    for (int i = d->pIBegin() + 1; i < d->pIEnd(); i++) {
        for (int j = d->pJBegin() + 1; j < d->pJEnd(); j++) {
            pRef(i, j) = 10 * i + j;
        }
    }
    std::cout << "p = ..." << std::endl;
    pRef.print();

    for (int i = d->rhsIBegin(); i <= d->rhsIEnd(); i++) {
        for (int j = d->rhsJBegin(); j <= d->rhsJEnd(); j++) {
            rhs(i, j) = (pRef(i + 1, j) - 2 * pRef(i, j) + pRef(i - 1, j)) / pow(d->dx(),2) +
                        (pRef(i, j + 1) - 2 * pRef(i, j) + pRef(i, j - 1)) / pow(d->dy(),2);
        }
    }
    std::cout << "rhs = ..." << std::endl;
    rhs.print();

    double epsilon = 0.001;
    int maximumNumberOfIterations = 100;
    double omega = 0.1;
    auto sor = SOR(static_cast<std::shared_ptr<Discretization>>(d), epsilon, maximumNumberOfIterations, omega);
    sor.solve();
    /*ASSERT_EQ(0.0, fv.interpolateAt(0.0, 0.0));
    ASSERT_EQ(1.0, fv.interpolateAt(1.0, 1.0));
    ASSERT_EQ(0.25, fv.interpolateAt(0.5, 0.5));
    ASSERT_EQ(0.5, fv.interpolateAt(1.0, 0.5));*/
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}