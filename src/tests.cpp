#include <gtest/gtest.h>
#include "storage/fieldvariable.h"
#include "discretization/2_central_differences.h"
#include "pressure_solver/1_sor.h"
#include "pressure_solver/1_gauss_seidel.h"
#include <cmath>

TEST(FieldVariableTest, InterpolationCheck1) {
    auto size = std::array<int, 2>{3, 3};
    auto meshWidth = std::array<double, 2>{1.0, 1.0};
    auto origin = std::array<double, 2>{meshWidth[0], meshWidth[1]};
    FieldVariable fv = FieldVariable(size, origin, meshWidth);
    auto nCells = std::array<int, 2>{3, 3};

    fv(1, 1) = 1.0;
    fv(2, 1) = 2.0;
    ASSERT_EQ(0.0, fv.interpolateAt(0.0, 0.0));
    ASSERT_EQ(1.0, fv.interpolateAt(1.0, 1.0));
    ASSERT_EQ(0.25, fv.interpolateAt(0.5, 0.5));
    ASSERT_EQ(0.5, fv.interpolateAt(1.0, 0.5));
}

TEST(FieldVariableTest, InterpolationCheck2) {
    auto size = std::array<int, 2>{3, 3};
    auto meshWidth = std::array<double, 2>{1.0, 1.0};
    auto origin = std::array<double, 2>{meshWidth[0] / 2.0, meshWidth[1] / 2.0};
    FieldVariable fv = FieldVariable(size, origin, meshWidth);
    auto nCells = std::array<int, 2>{3, 3};

    fv(1, 1) = 1.0;
    fv(2, 1) = 2.0;
    ASSERT_EQ(0.25, fv.interpolateAt(0.0, 0.0));
    ASSERT_EQ(0.75, fv.interpolateAt(1.0, 1.0));
    ASSERT_EQ(1.0, fv.interpolateAt(0.5, 0.5));
    ASSERT_EQ(2.0, fv.interpolateAt(1.5, 0.5));
}

TEST(FieldVariableTest, InterpolationCheck) {
    auto size = std::array<int, 2>{6+2, 5+2};
    auto meshWidth = std::array<double, 2>{1.0, 1.0};
    auto origin = std::array<double, 2>{meshWidth[0]/2., meshWidth[1]};

    FieldVariable v = FieldVariable(size, origin, meshWidth);
    for (int i = 1; i < size[0] - 1; i++) {
        for (int j = 1; j < size[1] - 1; j++) {
            v(i, j) = i * 10 + j;
        }
    }
    //v.print();
    
    FieldVariable v_interp = FieldVariable({size[0] - 1, size[1] - 1}, origin, meshWidth);
    for (int i = 0; i < size[0] - 1; i++) {
        for (int j = 0; j < size[1] - 1; j++) {
            double x = meshWidth[0] * i;
            double y = meshWidth[1] * j;
            v_interp(i, j) = v.interpolateAt(x, y);
        }
    }
    //v_interp.print();
}

TEST(SORTest, Test1) {
    auto nCells = std::array<int, 2>{3, 3};
    auto meshWidth = std::array<double, 2>{1, 1};
    auto origin = std::array<double, 2>{meshWidth[0] / 2.0, meshWidth[1] / 2.0};
    auto d = new CentralDifferences(nCells, meshWidth);
    auto dRef = new CentralDifferences(nCells, meshWidth);
    for (int i = d->pInteriorIBegin(); i < d->pInteriorIEnd(); i++) {
        for (int j = d->pInteriorJBegin(); j < d->pInteriorJEnd(); j++) {
            dRef->p(i, j) = i + 10 * j;
        }
    }
#ifndef NDEBUG
    //std::cout << "p = ..." << std::endl;
    //pRef.print();

    // std::cout << "p = ..." << std::endl;
    // d->p().print();

    // std::cout << "pref = ..." << std::endl;
    // pRef.print();
#endif

    for (int i = d->rhsInteriorIBegin(); i < d->rhsInteriorIEnd(); i++) {
        for (int j = d->rhsInteriorJBegin(); j < d->rhsInteriorJEnd(); j++) {
            d->rhs(i, j) = (dRef->p(i + 1, j) - 2.0 * dRef->p(i, j) + dRef->p(i - 1, j)) / pow(d->dx(),2) +
                  (dRef->p(i, j + 1) - 2.0 * dRef->p(i, j) + dRef->p(i, j - 1)) / pow(d->dy(),2);
        }
    }
    // std::cout << "rhs = ..." << std::endl;
    dRef->p().print();
    d->rhs().print();

    double epsilon = 0.001;
    int maximumNumberOfIterations = 1000;
    double omega = 2.0 / (1.0 + sin(3.1 * meshWidth[0]));
    auto sor = SOR(static_cast<std::shared_ptr<Discretization>>(d), epsilon, maximumNumberOfIterations, omega);
    sor.solve();
    /*std::cout << "Iterations "  << sor.iterations() << std::endl;
    std::cout << "Norm "  << sor.residualNorm() << std::endl;*/

    //std::cout << "p = ..." << std::endl;
    d->p().print();

    // check interior of the domain
    for (int i = d->pInteriorIBegin(); i < d->pInteriorIEnd(); i++) {
        for (int j = d->pInteriorJBegin(); j < d->pInteriorJEnd(); j++) {
            EXPECT_NEAR(d->p(i, j), dRef->p(i, j), epsilon);
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
        ASSERT_EQ(d->p(i, d->pJBegin()), d->p(i, d->pInteriorJBegin()));

        // check right boundary
        ASSERT_EQ(d->p(i, d->pJEnd() - 1), d->p(i, d->pInteriorJEnd() - 1));
    }
}


/*TEST(GaussSeidelTest, Test1) {
    auto nCells = std::array<int, 2>{5, 10};
    auto meshWidth = std::array<double, 2>{0.1, 0.05};
    auto origin = std::array<double, 2>{meshWidth[0] / 2.0, meshWidth[1] / 2.0};
    auto d = new CentralDifferences(nCells, meshWidth);
    auto pRef = FieldVariable({nCells[0] + 2, nCells[1] + 2}, origin, meshWidth);
    for (int i = d->pInteriorIBegin(); i < d->pInteriorIEnd(); i++) {
        for (int j = d->pInteriorJBegin(); j < d->pInteriorJEnd(); j++) {
            pRef(i, j) = 10 * i + j ;
        }
    }
#ifndef NDEBUG
    //std::cout << "p = ..." << std::endl;
    //pRef.print();
    //std::cout << "p = ..." << std::endl;
    //d->p().print();

    //std::cout << "pref = ..." << std::endl;
    //pRef.print();
#endif
    for (int i = d->rhsInteriorIBegin(); i < d->rhsInteriorIEnd(); i++) {
        for (int j = d->rhsInteriorJBegin(); j < d->rhsInteriorJEnd(); j++) {
            double rhs = (pRef(i + 1, j) - 2 * pRef(i, j) + pRef(i - 1, j)) / pow(d->dx(),2) +
                  (pRef(i, j + 1) - 2 * pRef(i, j) + pRef(i, j - 1)) / pow(d->dy(),2);
            // std::cout << rhs << std::endl;
            d->rhs(i, j) = rhs;
        }
    }
#ifndef NDEBUG
    //std::cout << "rhs = ..." << std::endl;
    //d->rhs().print();
#endif

    double epsilon = 0.001;
    int maximumNumberOfIterations = 1000;
    auto gs = GaussSeidel(static_cast<std::shared_ptr<Discretization>>(d), epsilon, maximumNumberOfIterations);
    gs.solve();
    // check interior of the domain
    for (int i = d->pInteriorIBegin(); i < d->pInteriorIEnd(); i++) {
        for (int j = d->pInteriorJBegin(); j < d->pInteriorJEnd(); j++) {
            EXPECT_NEAR(d->p(i, j), pRef(i, j), epsilon);
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
        ASSERT_EQ(d->p(i, d->pJBegin()), d->p(i, d->pInteriorJBegin()));

        // check right boundary
        ASSERT_EQ(d->p(i, d->pJEnd() - 1), d->p(i, d->pInteriorJEnd() - 1));
    }
}*/

TEST(DiscretizationTest, FirstOrder)
{
    auto nCells = std::array<int, 2>{100, 100};
    auto meshWidth = std::array<double, 2>{0.2, 0.1};

    CentralDifferences discr = CentralDifferences(nCells, meshWidth);

    // Assign values for p by p(x,y)=x*y
    for (int i = discr.pIBegin(); i < discr.pIEnd(); i++)
    {
        const double x = i * meshWidth[0];
        for (int j = discr.pJBegin(); j < discr.pJEnd(); j++)
        {
            const double y = j * meshWidth[1];
            discr.p(i, j) = x * y;
        
        }
    }

    // Check derivatives for p
    for (int i = discr.pInteriorIBegin(); i < discr.pInteriorIEnd(); i++)
    {
        const double x = i * meshWidth[0];
        for (int j = discr.pInteriorJBegin(); j < discr.pInteriorJEnd(); j++)
        {
            const double y = j * meshWidth[1];
            EXPECT_NEAR(discr.computeDpDx(i, j), y, 0.0000001);
            EXPECT_NEAR(discr.computeDpDy(i, j), x, 0.0000001);
        }
    }
}

TEST(DiscretizationTest, SecondOrder)
{
    auto nCells = std::array<int, 2>{100, 100};
    auto meshWidth = std::array<double, 2>{0.2, 0.1};

    CentralDifferences discr = CentralDifferences(nCells, meshWidth);

    // Assign values for v by v(x,y)=x*x*y
    for (int i = discr.vIBegin(); i < discr.vIEnd(); i++)
    {
        const double x = i * meshWidth[0];
        for (int j = discr.vJBegin(); j < discr.vJEnd(); j++)
        {
            const double y = j * meshWidth[1];
            discr.v(i, j) = x * x * y;
        }
    }
    // Assign values for u by u(x,y)=x*y*y
    for (int i = discr.uIBegin(); i < discr.uIEnd(); i++)
    {
        const double x = i * meshWidth[0];
        for (int j = discr.uJBegin(); j < discr.uJEnd(); j++)
        {
            const double y = j * meshWidth[1];
            discr.u(i, j) = x * y * y;
        }
    }

    // Check second order derivatives for u
    for (int i = discr.uInteriorIBegin(); i < discr.uInteriorIEnd(); i++)
    {
        const double x = i * meshWidth[0];
        for (int j = discr.uInteriorJBegin(); j < discr.uInteriorJEnd(); j++)
        {
            const double y = j * meshWidth[1];
            EXPECT_NEAR(discr.computeD2uDx2(i, j), 0, 0.0000001);
            EXPECT_NEAR(discr.computeD2uDy2(i, j), 2 * x, 0.0000001);
        }
    }
    // Check second order derivatives for v
    for (int i = discr.vInteriorIBegin(); i < discr.vInteriorIEnd(); i++)
    {
        const double x = i * meshWidth[0];
        for (int j = discr.vInteriorJBegin(); j < discr.vInteriorJEnd(); j++)
        {
            const double y = j * meshWidth[1];
            EXPECT_NEAR(discr.computeD2vDx2(i, j), 2 * y, 0.0000001);
            EXPECT_NEAR(discr.computeD2vDy2(i, j), 0, 0.0000001);
        }
    }
}

TEST(DiscretizationTest, FirstOrderSquared)
{
    auto nCells = std::array<int, 2>{100, 100};
    auto meshWidth = std::array<double, 2>{0.2, 0.1};

    CentralDifferences discr = CentralDifferences(nCells, meshWidth);

    // Assign values for v by v(x,y)=x*y
    for (int i = discr.vIBegin(); i < discr.vIEnd(); i++)
    {
        const double x = i * meshWidth[0];
        for (int j = discr.vJBegin(); j < discr.vJEnd(); j++)
        {
            const double y = j * meshWidth[1];
            discr.v(i, j) = x * y;
        }
    }
    // Assign values for u by u(x,y)=x*y
    for (int i = discr.uIBegin(); i < discr.uIEnd(); i++)
    {
        const double x = i * meshWidth[0];
        for (int j = discr.uJBegin(); j < discr.uJEnd(); j++)
        {
            const double y = j * meshWidth[1];
            discr.u(i, j) = x * y;
        }
    }

    // Check  derivatives for u^2
    for (int i = discr.uInteriorIBegin(); i < discr.uInteriorIEnd(); i++)
    {
        const double x = i * meshWidth[0];
        for (int j = discr.uInteriorJBegin(); j < discr.uInteriorJEnd(); j++)
        {
            const double y = j * meshWidth[1];
            EXPECT_NEAR(discr.computeDu2Dx(i, j), 2 * x * y * y, 0.0000001);
        }
    }
    // Check derivatives for v^2
    for (int i = discr.vInteriorIBegin(); i < discr.vInteriorIEnd(); i++)
    {
        const double x = i * meshWidth[0];
        for (int j = discr.vInteriorJBegin(); j < discr.vInteriorJEnd(); j++)
        {
            const double y = j * meshWidth[1];
            EXPECT_NEAR(discr.computeDv2Dy(i, j), 2 * y * x * x, 0.0000001);
        }
    }
}

TEST(DiscretizationTest, FirstOrderMixed)
{
    auto nCells = std::array<int, 2>{100, 100};
    auto meshWidth = std::array<double, 2>{0.2, 0.1};

    CentralDifferences discr = CentralDifferences(nCells, meshWidth);

    // Assign values for v by v(x,y)=x*y
    for (int i = discr.vIBegin(); i < discr.vIEnd(); i++)
    {
        const double x = i * meshWidth[0];
        for (int j = discr.vJBegin(); j < discr.vJEnd(); j++)
        {
            const double y = j * meshWidth[1];
            discr.v(i, j) = 1;
        }
    }
    // Assign values for u by u(x,y)=x*y
    for (int i = discr.uIBegin(); i < discr.uIEnd(); i++)
    {
        const double x = i * meshWidth[0];
        for (int j = discr.uJBegin(); j < discr.uJEnd(); j++)
        {
            const double y = j * meshWidth[1];
            discr.u(i, j) = 1;
        }
    }

    // Check  derivatives for u*v
    for (int i = discr.uInteriorIBegin(); i < discr.uInteriorIEnd(); i++)
    {
        const double x = i * meshWidth[0];
        for (int j = discr.uInteriorJBegin(); j < discr.uInteriorJEnd(); j++)
        {
            const double y = j * meshWidth[1];
            EXPECT_NEAR(discr.computeDuvDy(i, j), 0, 0.0000001);
        }
    }
    // Check derivatives for u*v
    for (int i = discr.vInteriorIBegin(); i < discr.vInteriorIEnd(); i++)
    {
        const double x = i * meshWidth[0];
        for (int j = discr.vInteriorJBegin(); j < discr.vInteriorJEnd(); j++)
        {
            const double y = j * meshWidth[1];
            EXPECT_NEAR(discr.computeDuvDx(i, j), 0, 0.0000001);
        }
    }
}

//################################################ Consistency Order ####################################################

TEST(DiscretizationTest, ConsistencyOrder)
{
    double p_exact = 12;
    int n_exact = pow(2, p_exact);
    auto nCells_exact = std::array<int, 2>{n_exact, n_exact};
    auto meshWidth_exact = std::array<double, 2>{1.0/n_exact, 1.0/n_exact};
    double re = 1000.0;
    double dt = 0.1;
    double alpha = 0.0;

    auto discr_exact = DonorCell(nCells_exact, meshWidth_exact, alpha);
    //auto discr_exact = CentralDifferences(nCells_exact, meshWidth_exact);

    // Assign values for v by v(x,y)=x*y
    for (int i = discr_exact.vIBegin(); i < discr_exact.vIEnd(); i++)
    {
        const double x = i * meshWidth_exact[0];
        for (int j = discr_exact.vJBegin(); j < discr_exact.vJEnd(); j++)
        {
            const double y = j * meshWidth_exact[1];
            discr_exact.v(i, j) = sin(x * 2.0 * M_PI) * cos(y * 2.0 * M_PI);
        }
    }
    // Assign values for u by u(x,y)=x*y
    for (int i = discr_exact.uIBegin(); i < discr_exact.uIEnd(); i++)
    {
        const double x = i * meshWidth_exact[0];
        for (int j = discr_exact.uJBegin(); j < discr_exact.uJEnd(); j++)
        {
            const double y = j * meshWidth_exact[1];
            discr_exact.u(i, j) = sin(x * 2.0 * M_PI) * cos(y * 2.0 * M_PI);
        }
    }

    for (int i = discr_exact.uInteriorIBegin(); i < discr_exact.uInteriorIEnd(); i++) {
        const double x = i * meshWidth_exact[0];
        for (int j = discr_exact.uInteriorJBegin(); j < discr_exact.uInteriorJEnd(); j++) {
            const double y = j * meshWidth_exact[1];
            double lap_u = discr_exact.computeD2uDx2(i,j) + discr_exact.computeD2uDy2(i,j);
            double conv_u = discr_exact.computeDu2Dx(i,j) + discr_exact.computeDuvDy(i,j);
            // double f_num = discr_exact.u(i,j) + dt * (lap_u / re - conv_u);
            double f_num = conv_u;
            discr_exact.f(i,j) = f_num;
        }
    }

    // Compute G in the interior of the domain
    for (int i = discr_exact.vInteriorIBegin(); i < discr_exact.vInteriorIEnd(); i++) {
        const double x = i * meshWidth_exact[0];
        for (int j = discr_exact.vInteriorJBegin(); j < discr_exact.vInteriorJEnd(); j++) {
            const double y = j * meshWidth_exact[1];
            double lap_v = discr_exact.computeD2vDx2(i,j) + discr_exact.computeD2vDy2(i,j);
            double conv_v = discr_exact.computeDv2Dy(i,j) + discr_exact.computeDuvDx(i,j);
            double s1 = sin(2*M_PI*y);
            double s2 = sin(2*M_PI*x);
            double g_exact = s2*s1-pow(M_PI,2)*s2*s1/125
                             -2*M_PI*cos(2*M_PI*x)*s2*pow(s1,2)/5-2*M_PI*cos(2*M_PI*y)*pow(s2,2)*s1/5;
            double g_num = discr_exact.v(i,j) + dt * (lap_v / re - conv_v);
            discr_exact.g(i,j) = g_num;
        }
    }

    for (int p = 4; p < p_exact; p++) {
        int n = pow(2, p);
        auto nCells = std::array<int, 2>{n, n};
        auto meshWidth = std::array<double, 2>{1.0/n, 1.0/n};

        auto discr = DonorCell(nCells, meshWidth, alpha);
        // auto discr = CentralDifferences(nCells, meshWidth);

        // Assign values for v by v(x,y)=x*y
        for (int i = discr.vIBegin(); i < discr.vIEnd(); i++)
        {
            const double x = i * meshWidth[0];
            for (int j = discr.vJBegin(); j < discr.vJEnd(); j++)
            {
                const double y = j * meshWidth[1];
                discr.v(i, j) = sin(x * 2.0 * M_PI) * cos(y * 2.0 * M_PI);
            }
        }
        // Assign values for u by u(x,y)=x*y
        for (int i = discr.uIBegin(); i < discr.uIEnd(); i++)
        {
            const double x = i * meshWidth[0];
            for (int j = discr.uJBegin(); j < discr.uJEnd(); j++)
            {
                const double y = j * meshWidth[1];
                discr.u(i, j) = sin(x * 2.0 * M_PI) * cos(y * 2.0 * M_PI);
            }
        }

        // Compute F in the interior of the domain
        double f_error_norm2 = 0.0;
        double f_norm2 = 0.0;
        for (int i = discr.uInteriorIBegin(); i < discr.uInteriorIEnd(); i++) {
            const double x = i * meshWidth[0];
            for (int j = discr.uInteriorJBegin(); j < discr.uInteriorJEnd(); j++) {
                const double y = j * meshWidth[1];
                double lap_u = discr.computeD2uDx2(i,j) + discr.computeD2uDy2(i,j);
                double conv_u = discr.computeDu2Dx(i,j) + discr.computeDuvDy(i,j);
                double s1 = sin(2*M_PI*y);
                double s2 = sin(2*M_PI*x);
                double f_exact = s2*s1-pow(M_PI,2)*s2*s1/125
                        -2*M_PI*cos(2*M_PI*x)*s2*pow(s1,2)/5-2*M_PI*cos(2*M_PI*y)*pow(s2,2)*s1/5;
                f_exact = discr_exact.f().interpolateAt(x, y);
                double f_num = discr.u(i,j) + dt * (lap_u / re - conv_u);
                //double f_num = conv_u;
                double f_error = fabs(f_exact - f_num);
                f_error_norm2 += pow(f_error, 2);
                f_norm2 += pow(f_exact, 2);
                discr.f(i,j) = f_error;
            }
        }

        // Compute G in the interior of the domain
        double g_error_norm2 = 0.0;
        double g_norm2 = 0.0;
        for (int i = discr.vInteriorIBegin(); i < discr.vInteriorIEnd(); i++) {
            const double x = i * meshWidth[0];
            for (int j = discr.vInteriorJBegin(); j < discr.vInteriorJEnd(); j++) {
                const double y = j * meshWidth[1];
                double lap_v = discr.computeD2vDx2(i,j) + discr.computeD2vDy2(i,j);
                double conv_v = discr.computeDv2Dy(i,j) + discr.computeDuvDx(i,j);
                double s1 = sin(2*M_PI*y);
                double s2 = sin(2*M_PI*x);
                double g_exact = s2*s1-pow(M_PI,2)*s2*s1/125
                                 -2*M_PI*cos(2*M_PI*x)*s2*pow(s1,2)/5-2*M_PI*cos(2*M_PI*y)*pow(s2,2)*s1/5;
                g_exact = discr_exact.g().interpolateAt(x, y);
                double g_num = discr.v(i,j) + dt * (lap_v / re - conv_v);
                double g_error = fabs(g_exact - g_num);
                g_error_norm2 += pow(g_error, 2);
                g_norm2 += pow(g_exact, 2);
                discr.g(i,j) = g_error;
            }
        }
        double f_error_norm = sqrt(f_error_norm2) / n;
        double f_norm = sqrt(f_norm2) / n_exact;
        double g_error_norm = sqrt(g_error_norm2) / n;
        double g_norm = sqrt(g_norm2) / n_exact;
        double f_rel_error = f_error_norm / f_norm;
        double g_rel_error = g_error_norm / g_norm;
        std::cout << "n = " << n << ", f_error = " << f_rel_error << ", g_error = " << g_rel_error << std::endl;
    }

}

//############################################ SOR TEST ##################################################################

/*TEST(SORTest, Test1) {
    auto nCells = std::array<int, 2>{3, 3};
    auto meshWidth = std::array<double, 2>{1, 1};
    auto origin = std::array<double, 2>{meshWidth[0] / 2.0, meshWidth[1] / 2.0};
    auto d = new CentralDifferences(nCells, meshWidth);
    auto dRef = new CentralDifferences(nCells, meshWidth);
    for (int i = d->pInteriorIBegin(); i < d->pInteriorIEnd(); i++) {
        for (int j = d->pInteriorJBegin(); j < d->pInteriorJEnd(); j++) {
            dRef->p(i, j) = i + 10 * j;
        }
    }
#ifndef NDEBUG
    //std::cout << "p = ..." << std::endl;
    //pRef.print();

    // std::cout << "p = ..." << std::endl;
    // d->p().print();

    // std::cout << "pref = ..." << std::endl;
    // pRef.print();
#endif

    for (int i = d->rhsInteriorIBegin(); i < d->rhsInteriorIEnd(); i++) {
        for (int j = d->rhsInteriorJBegin(); j < d->rhsInteriorJEnd(); j++) {
            d->rhs(i, j) = (dRef->p(i + 1, j) - 2.0 * dRef->p(i, j) + dRef->p(i - 1, j)) / pow(d->dx(),2) +
                  (dRef->p(i, j + 1) - 2.0 * dRef->p(i, j) + dRef->p(i, j - 1)) / pow(d->dy(),2);
        }
    }
    // std::cout << "rhs = ..." << std::endl;
    dRef->p().print();
    d->rhs().print();

    double epsilon = 0.001;
    int maximumNumberOfIterations = 1000;
    double omega = 2.0 / (1.0 + sin(3.1 * meshWidth[0]));
    auto sor = SOR(static_cast<std::shared_ptr<Discretization>>(d), epsilon, maximumNumberOfIterations, omega);
    sor.solve();*/
    /*std::cout << "Iterations "  << sor.iterations() << std::endl;
    std::cout << "Norm "  << sor.residualNorm() << std::endl;*/
    /*
    //std::cout << "p = ..." << std::endl;
    d->p().print();

    // check interior of the domain
    for (int i = d->pInteriorIBegin(); i < d->pInteriorIEnd(); i++) {
        for (int j = d->pInteriorJBegin(); j < d->pInteriorJEnd(); j++) {
            EXPECT_NEAR(d->p(i, j), dRef->p(i, j), epsilon);
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
        ASSERT_EQ(d->p(i, d->pJBegin()), d->p(i, d->pInteriorJBegin()));

        // check right boundary
        ASSERT_EQ(d->p(i, d->pJEnd() - 1), d->p(i, d->pInteriorJEnd() - 1));
    }
}*/

//######################################## Gauss Seidel Test ############################################################

/*TEST(GaussSeidelTest, Test1) {
    auto nCells = std::array<int, 2>{5, 10};
    auto meshWidth = std::array<double, 2>{0.1, 0.05};
    auto origin = std::array<double, 2>{meshWidth[0] / 2.0, meshWidth[1] / 2.0};
    auto d = new CentralDifferences(nCells, meshWidth);
    auto pRef = FieldVariable({nCells[0] + 2, nCells[1] + 2}, origin, meshWidth);
    for (int i = d->pInteriorIBegin(); i < d->pInteriorIEnd(); i++) {
        for (int j = d->pInteriorJBegin(); j < d->pInteriorJEnd(); j++) {
            pRef(i, j) = 10 * i + j ;
        }
    }
#ifndef NDEBUG
    //std::cout << "p = ..." << std::endl;
    //pRef.print();
    //std::cout << "p = ..." << std::endl;
    //d->p().print();

    //std::cout << "pref = ..." << std::endl;
    //pRef.print();
#endif
    for (int i = d->rhsInteriorIBegin(); i < d->rhsInteriorIEnd(); i++) {
        for (int j = d->rhsInteriorJBegin(); j < d->rhsInteriorJEnd(); j++) {
            double rhs = (pRef(i + 1, j) - 2 * pRef(i, j) + pRef(i - 1, j)) / pow(d->dx(),2) +
                  (pRef(i, j + 1) - 2 * pRef(i, j) + pRef(i, j - 1)) / pow(d->dy(),2);
            // std::cout << rhs << std::endl;
            d->rhs(i, j) = rhs;
        }
    }
#ifndef NDEBUG
    //std::cout << "rhs = ..." << std::endl;
    //d->rhs().print();
#endif

    double epsilon = 0.001;
    int maximumNumberOfIterations = 1000;
    auto gs = GaussSeidel(static_cast<std::shared_ptr<Discretization>>(d), epsilon, maximumNumberOfIterations);
    gs.solve();
    // check interior of the domain
    for (int i = d->pInteriorIBegin(); i < d->pInteriorIEnd(); i++) {
        for (int j = d->pInteriorJBegin(); j < d->pInteriorJEnd(); j++) {
            EXPECT_NEAR(d->p(i, j), pRef(i, j), epsilon);
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
        ASSERT_EQ(d->p(i, d->pJBegin()), d->p(i, d->pInteriorJBegin()));

        // check right boundary
        ASSERT_EQ(d->p(i, d->pJEnd() - 1), d->p(i, d->pInteriorJEnd() - 1));
    }
}*/

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
