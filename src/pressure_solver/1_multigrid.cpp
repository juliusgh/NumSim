//
// Created by nhornischer on 28.01.23.
//

#include "1_multigrid.h"

Multigrid::Multigrid(std::shared_ptr<Discretization> discretization,
                     double epsilon,
                     int maximumNumberOfIterations,
                     int maximumNumberOfLayers,
                     int gamma,
                     int mu,
                     int nu,
                     double theta):
        PressureSolver(discretization,epsilon,maximumNumberOfIterations),
        maximumNumberOfLayers_(maximumNumberOfLayers),
        gamma_(gamma),
        mu_(mu),
        nu_(nu),
        theta_(theta)
        {

}

void Multigrid::solve(){
    const double eps2 = pow(epsilon_, 2);
    int iteration = 0;
    do {
        iteration++;

        // Execute the multigrid cycle
        Array2D v_L = MG_Cycle(maximumNumberOfLayers_, (Array2D &) discretization_->p(), discretization_->rhs());

        // Update the pressure field (Possible that indices are wrong)
        for (int i = discretization_->pInteriorIBegin() - discretization_->pIBegin(); i < discretization_->pInteriorIBegin() - discretization_->pIBegin(); i++) {
            for (int j = discretization_->pInteriorJBegin() - discretization_->pJBegin(); j < discretization_->pInteriorJEnd() - discretization_->pJBegin(); j++) {
                discretization_->p(i, j) = v_L(i, j);
            }
        }
    } while (residualNorm() > eps2 && iteration < maximumNumberOfIterations_);
    iterations_ = iteration;
};

void Multigrid::smoothing(Array2D& v, int smoothing_steps, bool convergenceCheck) {

    int steps = 0;
    double residual = 0;
    do{
        steps++;
        // Apply a Gauss-Seidel relaxation
        for (int i = 0; i < v.size()[0]; i++) {
            for (int j = 0; j < v.size()[0]; j++) {
                v(i, j) = (1 - theta_) * v(i, j) +
                          theta_ * (1 / (2 * (1 / discretization_->dx() + 1 / discretization_->dy()))) *
                          (discretization_->rhs(i, j) - (v(i - 1, j) + v(i + 1, j)) / discretization_->dx() -
                           (v(i, j - 1) + v(i, j + 1)) / discretization_->dy());
                if (convergenceCheck) {
                    residual += discretization_->rhs(i, j) - (v(i - 1, j) + v(i + 1, j)) / discretization_->dx() - (v(i, j - 1) + v(i, j + 1)) / discretization_->dy();
                }
            }

        }
    } while ((convergenceCheck && (residualNorm() > epsilon_)) || (steps < smoothing_steps));


};

Array2D Multigrid::defect_calculation(Array2D v_L , Array2D rhs) {
    Array2D d_L = Array2D({v_L.size()[0], v_L.size()[1]});
    for (int i = 0; i < d_L.size()[0]; i++) {
        for (int j = 0; j < d_L.size()[1]; j++) {
            // Calculate the Laplacian of the 2D v_L with the 5 point stencil
            d_L(i, j) = rhs(i, j) - ((v_L(i-1 ,j) + v_L(i+1 ,j) + 2 * v_L(i,j))/ discretization_->dx() + (v_L(i ,j-1) + v_L(i ,j+1) + 2 * v_L(i,j)) / discretization_->dy());
        }
    }
    return d_L;
};

Array2D Multigrid::restriction(Array2D d_L) {
    Array2D d_l = Array2D({d_L.size()[0] / 2, d_L.size()[1] / 2});
    for (int i = 0; i < d_l.size()[0]; i++) {
        for (int j = 0; j < d_l.size()[1]; j++) {
            // Interpolation of the 4 surrounding points
            d_l(i, j) = 0.25 * (d_L(2 * i, 2 * j) + d_L(2 * i + 1, 2 * j) + d_L(2 * i, 2 * j + 1) + d_L(2 * i + 1, 2 * j + 1));
        }
    }
    return d_l;
};

Array2D Multigrid::prolongation(Array2D c_l) {
    Array2D c_L = Array2D({c_l.size()[0] * 2, c_l.size()[1] * 2});
    for (int i = 0; i < c_l.size()[0]; i++) {
        for (int j = 0; j < c_l.size()[1]; j++) {
            c_L(2 * i, 2 * j) = c_l(i, j);
            c_L(2 * i + 1, 2 * j) = c_l(i, j);
            c_L(2 * i, 2 * j + 1) = c_l(i, j);
            c_L(2 * i + 1, 2 * j + 1) = c_l(i, j);
        }
    }
    return c_L;
};

void Multigrid::correction(Array2D c_L, Array2D& v_L) {
    for (int i = 0; i < c_L.size()[0]; i++) {
        for (int j = 0; j < c_L.size()[1]; j++) {
            v_L(i, j) = v_L(i, j) + theta_ * c_L(i, j);
        }
    }
};

Array2D Multigrid::MG_Cycle(int layer, Array2D& v_L, Array2D rhs) {
    //Solve exactly on the coarsest grid
    if (layer==0){
        smoothing(v_L, mu_ + nu_, true);
    }
    else {
        // Presmoothing
        smoothing(v_L, mu_);

        // defect calculation on the finer grid
        Array2D d_L = defect_calculation(v_L, rhs);

        // restriction to the coarser grid
        Array2D d_l = restriction(d_L);

        // Initialize the correction to zero
        Array2D c_l = Array2D(d_l.size());

        //residual equation(MGCycle)
        for (int i = 1; i < gamma_; i++) {
            MG_Cycle(layer - 1, c_l, d_l);
        }
        //Prolongation
        Array2D c_L = prolongation(c_l);

        //Correction
        correction(c_L, v_L);

        //Postsmoothing
        smoothing(v_L, nu_);
    }
    return v_L;
};