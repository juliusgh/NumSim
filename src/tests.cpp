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

    for (int i = d->rhsIBegin(); i < d->rhsIEnd(); i++) {
        for (int j = d->rhsJBegin(); j < d->rhsJEnd(); j++) {
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
    for (int i = d->pInteriorIBegin(); i < d->pInteriorIEnd(); i++) {
        for (int j = d->pInteriorJBegin() + 1; j < d->pInteriorJEnd(); j++) {
            ASSERT_LE(abs(d->p(i, j) - pRef(i, j)), epsilon);
        }
    }

    // check boundary of the domain
    for (int j = d->pInteriorJBegin(); j < d->pInteriorJEnd(); j++) {

        // check left boundary
        ASSERT_EQ(d->p(d->pIBegin(), j), d->p(d->pInteriorIBegin(), j));

        // check right boundary
        ASSERT_EQ(d->p(d->pIEnd() - 1, j), d->p(d->pInteriorIEnd() - 1, j));
    }

    // check boundary of the domain
    for (int i = d->pInteriorIBegin(); i < d->pInteriorIEnd(); i++) {

        // check left boundary
        ASSERT_EQ(d->p(i, d->pJBegin() ), d->p(i, d->pInteriorJBegin()));

        // check right boundary
        ASSERT_EQ(d->p(i, d->pJEnd()), d->p(i, d->pInteriorJEnd() - 1));
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

    for (int i = d->rhsIBegin(); i < d->rhsIEnd(); i++) {
        for (int j = d->rhsJBegin(); j < d->rhsJEnd(); j++) {
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
    for (int i = d->pInteriorIBegin(); i < d->pInteriorIEnd(); i++) {
        for (int j = d->pInteriorJBegin(); j < d->pInteriorJEnd(); j++) {
            ASSERT_LE(abs(d->p(i, j) - pRef(i, j)), epsilon);
        }
    }

    // check boundary of the domain
    for (int j = d->pInteriorJBegin(); j < d->pInteriorJEnd(); j++) {

        // check left boundary
        ASSERT_EQ(d->p(d->pIBegin(), j), d->p(d->pInteriorIBegin(), j));

        // check right boundary
        ASSERT_EQ(d->p(d->pIEnd(), j), d->p(d->pInteriorIEnd() - 1, j));
    }

    // check boundary of the domain
    for (int i = d->pInteriorIBegin(); i < d->pInteriorIEnd(); i++) {

        // check left boundary
        ASSERT_EQ(d->p(i, d->pJBegin() ), d->p(i, d->pInteriorJBegin()));

        // check right boundary
        ASSERT_EQ(d->p(i, d->pJEnd()), d->p(i, d->pInteriorJEnd() - 1));
    }
}

TEST(DiscretizationTest, FirstOrder)
{
    auto nCells = std::array<int, 2>{4, 4};
    auto meshWidth = std::array<double, 2>{1.0, 1.0};

    CentralDifferences discr = CentralDifferences(nCells, meshWidth);

    // Assign values for p by p(x,y)=x*y
    for (int i = discr.pIBegin(); i <= discr.pIEnd(); i++)
    {
        const double x = i * meshWidth[0];
        for (int j = discr.pJBegin(); j <= discr.pJEnd(); j++)
        {
            const double y = j * meshWidth[1];
            discr.p(i, j) = x * y;
        }
    }

    // Check derivatives for p
    for (int i = discr.pInteriorIBegin(); i <= discr.pInteriorIEnd(); i++)
    {
        const double x = i * meshWidth[0];
        for (int j = discr.pInteriorJBegin(); j <= discr.pInteriorJEnd(); j++)
        {
            const double y = j * meshWidth[1];
            EXPECT_DOUBLE_EQ(discr.computeDpDx(i, j), y);
            EXPECT_DOUBLE_EQ(discr.computeDpDy(i, j), x);
        }
    }
}

TEST(DiscretizationTest, SecondOrder)
{
    auto nCells = std::array<int, 2>{3, 3};
    auto meshWidth = std::array<double, 2>{2.0, 1.0};

    CentralDifferences discr = CentralDifferences(nCells, meshWidth);

    // Assign values for v by v(x,y)=x*x*y
    for (int i = discr.vIBegin(); i <= discr.vIEnd(); i++)
    {
        const double x = i * meshWidth[0];
        for (int j = discr.vJBegin(); j <= discr.vJEnd(); j++)
        {
            const double y = j * meshWidth[1];
            discr.v(i, j) = x * x * y;
        }
    }
    // Assign values for u by u(x,y)=x*y*y
    for (int i = discr.uIBegin(); i <= discr.uIEnd(); i++)
    {
        const double x = i * meshWidth[0];
        for (int j = discr.uJBegin(); j <= discr.uJEnd(); j++)
        {
            const double y = j * meshWidth[1];
            discr.u(i, j) = x * y * y;
        }
    }

    // Check second order derivatives for u
    for (int i = discr.uInteriorIBegin(); i <= discr.uInteriorIEnd(); i++)
    {
        const double x = i * meshWidth[0];
        for (int j = discr.uInteriorJBegin(); j <= discr.uInteriorJEnd(); j++)
        {
            const double y = j * meshWidth[1];
            EXPECT_DOUBLE_EQ(discr.computeD2uDx2(i, j), 0);
            EXPECT_DOUBLE_EQ(discr.computeD2uDy2(i, j), 2 * x);
        }
    }
    // Check second order derivatives for v
    for (int i = discr.vInteriorIBegin(); i <= discr.vInteriorIEnd(); i++)
    {
        const double x = i * meshWidth[0];
        for (int j = discr.vInteriorJBegin(); j <= discr.vInteriorJEnd(); j++)
        {
            const double y = j * meshWidth[1];
            EXPECT_DOUBLE_EQ(discr.computeD2vDx2(i, j), 2 * y);
            EXPECT_DOUBLE_EQ(discr.computeD2vDy2(i, j), 0);
        }
    }
}

TEST(DiscretizationTest, FirstOrderSquared)
{
    auto nCells = std::array<int, 2>{3, 3};
    auto meshWidth = std::array<double, 2>{2.0, 1.0};

    CentralDifferences discr = CentralDifferences(nCells, meshWidth);

    // Assign values for v by v(x,y)=x*y
    for (int i = discr.vIBegin(); i <= discr.vIEnd(); i++)
    {
        const double x = i * meshWidth[0];
        for (int j = discr.vJBegin(); j <= discr.vJEnd(); j++)
        {
            const double y = j * meshWidth[1];
            discr.v(i, j) = x * y;
        }
    }
    // Assign values for u by u(x,y)=x*y
    for (int i = discr.uIBegin(); i <= discr.uIEnd(); i++)
    {
        const double x = i * meshWidth[0];
        for (int j = discr.uJBegin(); j <= discr.uJEnd(); j++)
        {
            const double y = j * meshWidth[1];
            discr.u(i, j) = x * y;
        }
    }

    // Check  derivatives for u^2
    for (int i = discr.uInteriorIBegin(); i <= discr.uInteriorIEnd(); i++)
    {
        const double x = i * meshWidth[0];
        for (int j = discr.uInteriorJBegin(); j <= discr.uInteriorJEnd(); j++)
        {
            const double y = j * meshWidth[1];
            EXPECT_DOUBLE_EQ(discr.computeDu2Dx(i, j), 2 * x * y * y);
        }
    }
    // Check derivatives for v^2
    for (int i = discr.vInteriorIBegin(); i <= discr.vInteriorIEnd(); i++)
    {
        const double x = i * meshWidth[0];
        for (int j = discr.vInteriorJBegin(); j <= discr.vInteriorJEnd(); j++)
        {
            const double y = j * meshWidth[1];
            EXPECT_DOUBLE_EQ(discr.computeDv2Dy(i, j), 2 * y * x * x);
        }
    }
}

TEST(DiscretizationTest, FirstOrderMixed)
{
    auto nCells = std::array<int, 2>{3, 3};
    auto meshWidth = std::array<double, 2>{2.0, 1.0};

    CentralDifferences discr = CentralDifferences(nCells, meshWidth);

    // Assign values for v by v(x,y)=x*y
    for (int i = discr.vIBegin(); i <= discr.vIEnd(); i++)
    {
        const double x = i * meshWidth[0];
        for (int j = discr.vJBegin(); j <= discr.vJEnd(); j++)
        {
            const double y = j * meshWidth[1];
            discr.v(i, j) = 1;
        }
    }
    // Assign values for u by u(x,y)=x*y
    for (int i = discr.uIBegin(); i <= discr.uIEnd(); i++)
    {
        const double x = i * meshWidth[0];
        for (int j = discr.uJBegin(); j <= discr.uJEnd(); j++)
        {
            const double y = j * meshWidth[1];
            discr.u(i, j) = 1;
        }
    }

    // Check  derivatives for u*v
    for (int i = discr.uInteriorIBegin(); i <= discr.uInteriorIEnd(); i++)
    {
        const double x = i * meshWidth[0];
        for (int j = discr.uInteriorJBegin(); j <= discr.uInteriorJEnd(); j++)
        {
            const double y = j * meshWidth[1];
            EXPECT_DOUBLE_EQ(discr.computeDuvDy(i, j), 0);
        }
    }
    // Check derivatives for u*v
    for (int i = discr.vInteriorIBegin(); i <= discr.vInteriorIEnd(); i++)
    {
        const double x = i * meshWidth[0];
        for (int j = discr.vInteriorJBegin(); j <= discr.vInteriorJEnd(); j++)
        {
            const double y = j * meshWidth[1];
            EXPECT_DOUBLE_EQ(discr.computeDuvDx(i, j), 0);
        }
    }
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
