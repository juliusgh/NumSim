#include <gtest/gtest.h>
#include "storage/fieldvariable.h"
#include "discretization/2_central_differences.h"
#include "pressure_solver/sor.h"
#include "pressure_solver/gauss_seidel.h"
#include <cmath>

TEST(FieldVariableTest, ValueCheck) {
    auto size = std::array<int, 2>{3, 3};
    auto origin = std::array<double, 2>{0.0, 0.0};
    auto meshWidth = std::array<double, 2>{1.0, 1.0};
    FieldVariable fv = FieldVariable(size, origin, meshWidth);
    auto nCells = std::array<int, 2>{3, 3};

    fv(1, 1) = 1.0;
    fv(2, 1) = 2.0;
    ASSERT_EQ(0.0, fv.interpolateAt(0.0, 0.0));
    ASSERT_EQ(1.0, fv.interpolateAt(1.0, 1.0));
    ASSERT_EQ(0.25, fv.interpolateAt(0.5, 0.5));
    ASSERT_EQ(0.5, fv.interpolateAt(1.0, 0.5));
}

TEST(SORTest, Test1) {
    auto nCells = std::array<int, 2>{5, 10};
    auto meshWidth = std::array<double, 2>{0.1, 0.05};
    auto origin = std::array<double, 2>{meshWidth[0] / 2.0, meshWidth[1] / 2.0};
    auto d = new CentralDifferences(nCells, meshWidth);
    auto pRef = FieldVariable({nCells[0] + 2, nCells[1] + 2}, origin, meshWidth);
    for (int i = d->pIBegin() + 1; i < d->pIEnd(); i++) {
        for (int j = d->pJBegin() + 1; j < d->pJEnd(); j++) {
            pRef(i, j) = 10 * i + j;
        }
    }
    //std::cout << "p = ..." << std::endl;
    //pRef.print();

    for (int i = d->rhsIBegin(); i <= d->rhsIEnd(); i++) {
        for (int j = d->rhsJBegin(); j <= d->rhsJEnd(); j++) {
            d->rhs(i, j) = (pRef(i + 1, j) - 2 * pRef(i, j) + pRef(i - 1, j)) / pow(d->dx(),2) +
                           (pRef(i, j + 1) - 2 * pRef(i, j) + pRef(i, j - 1)) / pow(d->dy(),2);
        }
    }
    //std::cout << "rhs = ..." << std::endl;
    //d->rhs().print();

    double epsilon = 0.001;
    int maximumNumberOfIterations = 1000;
    double omega = 2 / (1 + sin(3.1 * meshWidth[0]));
    auto sor = SOR(static_cast<std::shared_ptr<Discretization>>(d), epsilon, maximumNumberOfIterations, omega);
    sor.solve();

    //std::cout << "p = ..." << std::endl;
    //d->p().print();

    // check interior of the domain
    for (int i = d->pIBegin() + 1; i < d->pIEnd(); i++) {
        for (int j = d->pJBegin() + 1; j < d->pJEnd(); j++) {
            ASSERT_LE(abs(d->p(i, j) - pRef(i, j)), epsilon);
        }
    }

    // check boundary of the domain
    for (int j = d->pJBegin() + 1; j < d->pJEnd(); j++) {

        // check left boundary
        ASSERT_EQ(d->p(d->pIBegin(), j), d->p(d->pIBegin() + 1, j));

        // check right boundary
        ASSERT_EQ(d->p(d->pIEnd(), j), d->p(d->pIEnd() - 1, j));
    }

    // check boundary of the domain
    for (int i = d->pIBegin() + 1; i < d->pIEnd(); i++) {

        // check left boundary
        ASSERT_EQ(d->p(i, d->pJBegin() ), d->p(i, d->pJBegin() + 1));

        // check right boundary
        ASSERT_EQ(d->p(i, d->pJEnd()), d->p(i, d->pJEnd() - 1));
    }
}


TEST(GaussSeidelTest, Test1) {
    auto nCells = std::array<int, 2>{5, 10};
    auto meshWidth = std::array<double, 2>{0.1, 0.05};
    auto origin = std::array<double, 2>{meshWidth[0] / 2.0, meshWidth[1] / 2.0};
    auto d = new CentralDifferences(nCells, meshWidth);
    auto pRef = FieldVariable({nCells[0] + 2, nCells[1] + 2}, origin, meshWidth);
    for (int i = d->pIBegin() + 1; i < d->pIEnd(); i++) {
        for (int j = d->pJBegin() + 1; j < d->pJEnd(); j++) {
            pRef(i, j) = 10 * i + j;
        }
    }
    //std::cout << "p = ..." << std::endl;
    //pRef.print();

    for (int i = d->rhsIBegin(); i <= d->rhsIEnd(); i++) {
        for (int j = d->rhsJBegin(); j <= d->rhsJEnd(); j++) {
            d->rhs(i, j) = (pRef(i + 1, j) - 2 * pRef(i, j) + pRef(i - 1, j)) / pow(d->dx(),2) +
                           (pRef(i, j + 1) - 2 * pRef(i, j) + pRef(i, j - 1)) / pow(d->dy(),2);
        }
    }
    //std::cout << "rhs = ..." << std::endl;
    //d->rhs().print();

    double epsilon = 0.001;
    int maximumNumberOfIterations = 1000;

    auto gs = GaussSeidel(static_cast<std::shared_ptr<Discretization>>(d), epsilon, maximumNumberOfIterations);
    gs.solve();

    //std::cout << "p = ..." << std::endl;
    //d->p().print();

    // check interior of the domain
    for (int i = d->pIBegin() + 1; i < d->pIEnd(); i++) {
        for (int j = d->pJBegin() + 1; j < d->pJEnd(); j++) {
            ASSERT_LE(abs(d->p(i, j) - pRef(i, j)), epsilon);
        }
    }

    // check boundary of the domain
    for (int j = d->pJBegin() + 1; j < d->pJEnd(); j++) {

        // check left boundary
        ASSERT_EQ(d->p(d->pIBegin(), j), d->p(d->pIBegin() + 1, j));

        // check right boundary
        ASSERT_EQ(d->p(d->pIEnd(), j), d->p(d->pIEnd() - 1, j));
    }

    // check boundary of the domain
    for (int i = d->pIBegin() + 1; i < d->pIEnd(); i++) {

        // check left boundary
        ASSERT_EQ(d->p(i, d->pJBegin() ), d->p(i, d->pJBegin() + 1));

        // check right boundary
        ASSERT_EQ(d->p(i, d->pJEnd()), d->p(i, d->pJEnd() - 1));
    }
}



int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}