#include <iostream>
#include <fstream>
#include <utility>
#include "discretization/0_staggered_grid.h"

/**
 * Implement staggered grid, providing a variety of parameters
 * @param partitioning: encapsulate functionality corresponding to subdomain handling
 * @param meshWidth: cell width in all directions
 * @param settings: information about the settings received from parameter file
 */
StaggeredGrid::StaggeredGrid(const std::shared_ptr<Partitioning> &partitioning,
                             std::array<double, 2> meshWidth,
                             Settings settings) :
        partitioning_(partitioning),
        settings_(settings),
        nCells_(partitioning->nCellsLocal()),
        meshWidth_(meshWidth),
        marker_(pSize()),
        f_(uSize(), {meshWidth[0], meshWidth[1] / 2.}, meshWidth),
        g_(vSize(), {meshWidth[0] / 2., meshWidth[1]}, meshWidth),
        p_(pSize(), {meshWidth[0] / 2., meshWidth[1] / 2.}, meshWidth),
        rhs_(rhsSize(), {meshWidth[0] / 2., meshWidth[1] / 2.}, meshWidth),
        u_(uSize(), {meshWidth[0], meshWidth[1] / 2.}, meshWidth),
        v_(vSize(), {meshWidth[0] / 2., meshWidth[1]}, meshWidth),
        uLast_(uSize(), {meshWidth[0], meshWidth[1] / 2.}, meshWidth),
        vLast_(vSize(), {meshWidth[0] / 2., meshWidth[1]}, meshWidth),
        t_(tSize(), {meshWidth[0] / 2., meshWidth[1] / 2.}, meshWidth),
        tobs_(tSize(), {meshWidth[0] / 2., meshWidth[1] / 2.}, meshWidth),
        q_(tSize(), {meshWidth[0] / 2., meshWidth[1] / 2.}, meshWidth),
        particles_(0) {

    // set markers:
    marker_.setToFluid();
    q_.setToZero();
    tobs_.setToZero();

    // open input file for markers
    ifstream file(settings_.domainfile_path, ios::in);
    if (!file.is_open()) {
        cout << "Could not open domain file \"" << settings_.domainfile_path << "\"." << endl;
        cout << "Using lid driven cavity as default scenario on domain" << endl;
        // lid driven cavity scenario as default
        for (int i = pIBegin(); i < pIEnd(); i++) {
            // bottom
            marker(i, pJBegin()) = MARKER::NOSLIP;
            // top
            marker(i, pJEnd() - 1) = MARKER::INFLOW;
        }
        for (int j = pJBegin(); j < pJEnd(); j++) {
            // left
            marker(pIBegin(), j) = MARKER::NOSLIP;
            // right
            marker(pIEnd() - 1, j) = MARKER::NOSLIP;
        }
        marker_.print();
    } else {
        // read markers from file{
        cout << "Reading domain file \"" << settings_.domainfile_path << "\"." << endl;
        for (int j = 0;; j++) {
            string line;
            getline(file, line);

            if (file.eof())
                break;
            // set general marker
            for (int i = 0; i < line.size(); ++i) {
                int mi = i + pIBegin();
                int mj = pJEnd() - j - 1;
                switch (line[i]) {
                    case ' ':
                        marker(mi, mj) = MARKER::FREE;
                        break;
                    case '-':
                        marker(mi, mj) = MARKER::FLUID;
                        break;
                    case '*':
                        marker(mi, mj) = MARKER::FLUID;
                        q(mi, mj) = settings_.heatMagnitude;
                        break;
                    case 'i':
                        marker(mi, mj) = MARKER::INFLOW;
                        break;
                    case 'o':
                        marker(mi, mj) = MARKER::OUTFLOW;
                        break;
                    case 'n':
                        marker(mi, mj) = MARKER::NOSLIP;
                        break;
                    case 'x':
                        marker(mi, mj) = MARKER::OBSTACLE;
                        break;
                    case '#':
                        marker(mi, mj) = MARKER::OBSTACLE;
                        q(mi, mj) = settings_.heatMagnitude;
                        break;
                    case 'h':
                        marker(mi, mj) = MARKER::OBSTACLE;
                        tobs(mi, mj) = settings_.obstacleHot;
                        break;
                    case 'c':
                        marker(mi, mj) = MARKER::OBSTACLE;
                        tobs(mi, mj) = settings_.obstacleCold;
                        break;
                    default:
                        break;
                }
            }
        }
    }
    // check two-cell criterion. If not fulfilled end simulation
    for (int i = pInteriorIBegin(); i < pInteriorIEnd(); i++) {
        for (int j = pInteriorJBegin(); j < pInteriorJEnd(); j++) {
            if (marker(i,j) == MARKER::OBSTACLE){
                if ((j != pInteriorJBegin() and j != pInteriorJEnd()) and marker(i,j - 1) == MARKER::FLUID and marker(i,j + 1) == MARKER::FLUID){
                    throw std::domain_error("Given domain does not fulfill two cell criterion. Exiting programm:(");
                }
                if ((i != pInteriorIBegin() and i != pInteriorIEnd()) and marker(i - 1,j) == MARKER::FLUID and marker(i + 1,j) == MARKER::FLUID){
                    throw std::domain_error("Given domain does not fulfill two cell criterion. Exiting programm:(");
                }
            }
        }
    }               
    // TODO: activate commented out code when free surface is completed
    //setInitialParticles();
    setInitialTemperature();
    setObstacleMarkers();
    //updateCellTypes();
};

/**
 * miscellaneous
 */

/**
 * create and set particles uniformly distributed in fluid cells for free surface on a finer grid
 */
void StaggeredGrid::setInitialParticles() {
    // count fluid cells
    int fluidCells = 0;
    for (int i = pIBegin(); i < pIEnd(); i++) {
        for (int j = pJBegin(); j < pJEnd(); j++) {
            if (marker(i, j) == FLUID) {
                fluidCells++;
            }
        }
    }
    // place particle in fluid cell
    int particleNumber = settings_.particleRefinement * settings_.particleRefinement * fluidCells;
    particles_ = Particle2D(particleNumber);
    int k = 0;
    for (int i = pIBegin(); i < pIEnd(); i++) {
        for (int j = pJBegin(); j < pJEnd(); j++) {
            if (marker(i, j) != FLUID) {
                continue;
            }
            double originX = i * dx();
            double originY = j * dy();
            int partNumX = settings_.particleRefinement;
            int partNumY = settings_.particleRefinement;
            double hx = dx() / (partNumX + 1);
            double hy = dy() / (partNumY + 1);
            for (int px = 0; px < partNumX; px++) {
                for (int py = 0; py < partNumY; py++) {
                    particlePosX(k) = originX + (px + 1) * hx;
                    particlePosY(k) = originY + (py + 1) * hy;
                    k++;
                }
            }
        }
    }
};

/**
 * set temperature everywhere to an initial value (same value everywhere)
 */
void StaggeredGrid::setInitialTemperature() {
    for (int i = tIBegin(); i < tIEnd(); i++) {
        for (int j = tJBegin(); j < tJEnd(); j++) {
            t(i, j) = settings_.initialTemp;
        }
    }
};

/**
 * Obstacles and boundary conditions
 */

/**
 * specify obstacle states within an obstacle, because different positions require different calculations
 */
void StaggeredGrid::setObstacleMarkers() {
    // set obstacle markers
    for (int i = pIBegin(); i < pIEnd(); i++) {
        for (int j = pJBegin(); j < pJEnd(); j++) {
            if (marker(i, j) != MARKER::OBSTACLE)
                continue;
            // check right neighbor
            if (marker(i + 1, j) == MARKER::FLUID) {
                // check top neighbor
                if (marker(i, j + 1) == MARKER::FLUID) {
                    marker(i, j) = MARKER::OBSTACLE_RIGHT_TOP;
                }
                    // check bottom neighbor
                else if (marker(i, j - 1) == MARKER::FLUID) {
                    marker(i, j) = MARKER::OBSTACLE_RIGHT_BOTTOM;
                } else {
                    marker(i, j) = MARKER::OBSTACLE_RIGHT;
                }
            }
                // check left neighbor
            else if (marker(i - 1, j) == MARKER::FLUID) {
                // check top neighbor
                if (marker(i, j + 1) == MARKER::FLUID) {
                    marker(i, j) = MARKER::OBSTACLE_LEFT_TOP;
                }
                    // check bottom neighbor
                else if (marker(i, j - 1) == MARKER::FLUID) {
                    marker(i, j) = MARKER::OBSTACLE_LEFT_BOTTOM;
                } else {
                    marker(i, j) = MARKER::OBSTACLE_LEFT;
                }
            } else {
                // check top neighbor
                if (marker(i, j + 1) == MARKER::FLUID) {
                    marker(i, j) = MARKER::OBSTACLE_TOP;
                }
                    // check bottom neighbor
                else if (marker(i, j - 1) == MARKER::FLUID) {
                    marker(i, j) = MARKER::OBSTACLE_BOTTOM;
                }
            }
        }
    }
}

/**
 * set (preliminary) velocities v, u, F and G to zero in obstacles
 */
void StaggeredGrid::setObstacleValues() {
    for (int i = pInteriorIBegin(); i < pInteriorIEnd(); i++) {
        for (int j = pInteriorJBegin(); j < pInteriorJEnd(); j++) {
            switch (marker(i, j)) {
                case OBSTACLE:
                case OBSTACLE_LEFT:
                case OBSTACLE_RIGHT:
                case OBSTACLE_TOP:
                case OBSTACLE_BOTTOM:
                case OBSTACLE_LEFT_TOP:
                case OBSTACLE_RIGHT_TOP:
                case OBSTACLE_LEFT_BOTTOM:
                case OBSTACLE_RIGHT_BOTTOM:
                    f(i, j) = u(i, j) = 0.0;
                    g(i, j) = v(i, j) = 0.0;
                    break;
                default:
                    break;
            }
        }
    }
}

/**
 * set (preliminary) velocities v, u, F and G at obstacle boundaries according to derivation extended from the lecture
 * as well as conditions at the domains boundary values according to specified boundary conditions
 */
void StaggeredGrid::applyBoundaryVelocities() {
    for (int i = pInteriorIBegin(); i < pInteriorIEnd(); i++) {
        for (int j = pInteriorJBegin(); j < pInteriorJEnd(); j++) {
            switch (marker(i, j)) {
                case OBSTACLE:
                    f(i, j) = u(i, j) = 0.0;
                    g(i, j) = v(i, j) = 0.0;
                    break;
                case OBSTACLE_LEFT:
                    f(i, j) = u(i, j) = 0.0;
                    f(i - 1, j) = u(i - 1, j) = 0.0;
                    g(i, j) = v(i, j) = -v(i - 1, j);
                    break;
                case OBSTACLE_RIGHT:
                    f(i, j) = u(i, j) = 0.0;
                    g(i, j) = v(i, j) = -v(i + 1, j);
                    break;
                case OBSTACLE_TOP:
                    f(i, j) = u(i, j) = -u(i, j + 1);
                    g(i, j) = v(i, j) = 0.0;
                    break;
                case OBSTACLE_BOTTOM:
                    f(i, j) = u(i, j) = -u(i, j - 1);
                    g(i, j) = v(i, j) = 0.0;
                    g(i, j - 1) = v(i, j - 1) = 0.0;
                    break;
                case OBSTACLE_LEFT_TOP:
                    f(i, j) = u(i, j) = 0.0;
                    g(i, j) = v(i, j) = 0.0;
                    f(i - 1, j) = u(i - 1, j) = 0.0;
                    break;
                case OBSTACLE_RIGHT_TOP:
                    f(i, j) = u(i, j) = 0.0;
                    g(i, j) = v(i, j) = 0.0;
                    break;
                case OBSTACLE_LEFT_BOTTOM:
                    f(i, j) = u(i, j) = -u(i, j - 1);
                    f(i - 1, j) = u(i - 1, j) = 0.0;
                    g(i, j) = v(i, j) = -v(i - 1, j);
                    break;
                case OBSTACLE_RIGHT_BOTTOM:
                    f(i, j) = u(i, j) = 0.0;
                    g(i, j) = v(i, j) = -v(i + 1, j);
                    g(i, j - 1) = v(i, j - 1) = 0.0;
                    break;
                default:
                    break;
            }
        }
    }

    // set boundary values for u and v at bottom and top side (lower priority)
    for (int i = pIBegin(); i < pIEnd(); i++) {
        // set boundary values at bottom side
        switch (marker(i, pJBegin())) {
            case NOSLIP:
                if (i < pIEnd() - 1) {
                    f(i, uJBegin()) = u(i, uJBegin()) = -u(i, uInteriorJBegin());
                }
                g(i, vJBegin()) = v(i, vJBegin()) = 0.0;
                break;
            case INFLOW:
                if (i < pIEnd() - 1) {
                    f(i, uJBegin()) =
                    u(i, uJBegin()) = 2.0 * settings_.dirichletBcBottom[0]
                                              - u(i, uInteriorJBegin());
                }
                g(i, vJBegin()) =
                v(i, vJBegin()) = settings_.dirichletBcBottom[1];
                break;
            case OUTFLOW:
                if (i < pIEnd() - 1) {
                    f(i, uJBegin()) =
                    u(i, uJBegin()) = u(i, uInteriorJBegin());
                }
                g(i, vJBegin()) =
                v(i, vJBegin()) = v(i, vInteriorJBegin());
                break;
            default:
                break;
        }

        // set boundary values for u at top side
        switch (marker(i, pJEnd() - 1)) {
            case NOSLIP:
                if (i < pIEnd() - 1) {
                    f(i, uJEnd() - 1) = u(i, uJEnd() - 1) = -u(i, uInteriorJEnd() - 1);
                }
                g(i, vJEnd() - 1) = v(i, vJEnd() - 1) = 0.0;
                break;
            case INFLOW:
                if (i < pIEnd() - 1) {
                    f(i, uJEnd() - 1) = u(i, uJEnd() - 1) =
                            2.0 * settings_.dirichletBcTop[0] - u(i, uInteriorJEnd() - 1);
                }
                g(i, vJEnd() - 1) = v(i, vJEnd() - 1) = settings_.dirichletBcTop[1];
                break;
            case OUTFLOW:
                if (i < pIEnd() - 1) {
                    f(i, uJEnd() - 1) = u(i, uJEnd() - 1) = u(i, uInteriorJEnd() - 1);
                }
                g(i, vJEnd() - 1) = v(i, vJEnd() - 1) = v(i, vInteriorJEnd() - 1);
                break;
            default:
                break;
        }
    }

    // set boundary values for u and v at left and right side (higher priority)
    for (int j = pJBegin(); j < pJEnd(); j++) {
        // set boundary values for u at left side
        switch (marker(pIBegin(), j)) {
            case NOSLIP:
                f(uIBegin(), j) = u(uIBegin(), j) = 0.0;
                if (j < pJEnd() - 1) {
                    g(vIBegin(), j) = v(vIBegin(), j) = -v(vInteriorIBegin(), j);
                }
                break;
            case INFLOW:
                f(uIBegin(), j) = u(uIBegin(), j) = settings_.dirichletBcLeft[0];
                if (j < pJEnd() - 1) {
                    g(vIBegin(), j) = v(vIBegin(), j) = 2.0 * settings_.dirichletBcLeft[1]
                                                        - v(vInteriorIBegin(), j);
                }
                break;
            case OUTFLOW:
                f(uIBegin(), j) = u(uIBegin(), j) = u(uInteriorIBegin(), j);
                if (j < pJEnd() - 1) {
                    v(vIBegin(), j) = v(vInteriorIBegin(), j);
                }
                break;
            default:
                break;
        }
        // set boundary values for u at right side
        switch (marker(pIEnd() - 1, j)) {
            case NOSLIP:
                f(uIEnd() - 1, j) = u(uIEnd() - 1, j) = 0.0;
                if (j < pJEnd() - 1) {
                    g(vIEnd() - 1, j) = v(vIEnd() - 1, j) = -v(vInteriorIEnd() - 1, j);
                }
                break;
            case INFLOW:
                f(uIEnd() - 1, j) = u(uIEnd() - 1, j) = settings_.dirichletBcRight[0];
                if (j < pJEnd() - 1) {
                    g(vIEnd() - 1, j) = v(vIEnd() - 1, j) = settings_.dirichletBcRight[1]
                                                            - u(vInteriorIEnd() - 1, j);
                }
                break;
            case OUTFLOW:
                f(uIEnd() - 1, j) = u(uIEnd() - 1, j) = u(uInteriorIEnd() - 1, j);
                if (j < pJEnd() - 1) {
                    g(vIEnd() - 1, j) = v(vIEnd() - 1, j) = v(vInteriorIEnd() - 1, j);
                }
                break;
            default:
                break;
        }
    }
};

/**
 * Set pressure values at obstacle boundaries as well as at the domain boundaries according to specified boundary conditions
 */
void StaggeredGrid::applyBoundaryPressure() {
    for (int i = pIBegin(); i < pIEnd(); i++) {
        for (int j = pJBegin(); j < pJEnd(); j++) {
            switch (marker(i, j)) {
                case OBSTACLE_LEFT:
                    p(i, j) = p(i - 1, j);
                    break;
                case OBSTACLE_RIGHT:
                    p(i, j) = p(i + 1, j);
                    break;
                case OBSTACLE_TOP:
                    p(i, j) = p(i,j + 1);
                    break;
                case OBSTACLE_BOTTOM:
                    p(i, j) = p(i, j - 1);
                    break;
                case OBSTACLE_LEFT_TOP:
                    p(i, j) = (p(i - 1, j) + p(i,j + 1)) / 2.0;
                    break;
                case OBSTACLE_RIGHT_TOP:
                    p(i, j) = (p(i + 1, j) + p(i,j + 1)) / 2.0;
                    break;
                case OBSTACLE_LEFT_BOTTOM:
                    p(i, j) = (p(i - 1, j) + p(i,j - 1)) / 2.0;
                    break;
                case OBSTACLE_RIGHT_BOTTOM:
                    p(i, j) = (p(i + 1, j) + p(i,j - 1)) / 2.0;
                    break;
                default:
                    break;
            }
        }
    }

    // set boundary values for p at bottom and top side (lower priority)
    for (int i = pIBegin(); i < pIEnd(); i++) {
        // set boundary values at bottom side
        switch (marker(i, pJBegin())) {
            case INFLOW:
            case NOSLIP:
                p(i, pJBegin()) = p(i, pInteriorJBegin());
                break;
            case OUTFLOW:
                p(i, uJBegin()) = -p(i, uInteriorJBegin());
                break;
            default:
                break;
        }
        // set boundary values for p at top side
        switch (marker(i, pJEnd() - 1)) {
            case NOSLIP:
            case INFLOW:
                p(i, pJEnd() - 1) = p(i, pInteriorJEnd() - 1);
                break;
            case OUTFLOW:
                p(i, uJEnd() - 1) = -p(i, uInteriorJEnd() - 1);
                break;
            default:
                break;
        }
    }

    // set boundary values for p at left and right side (higher priority)
    for (int j = pJBegin(); j < pJEnd(); j++) {
        // set boundary values for p at left side
        switch (marker(pIBegin(), j)) {
            case NOSLIP:
            case INFLOW:
                p(pIBegin(), j) = p(pInteriorIBegin(), j);
                break;
            case OUTFLOW:
                p(pIBegin(), j) = -p(uInteriorIBegin(), j);
                break;
            default:
                break;
        }
        // set boundary values for p at right side
        switch (marker(pIEnd() - 1, j)) {
            case NOSLIP:
            case INFLOW:
                p(pIEnd() - 1, j) = p(pInteriorIEnd() - 1, j);
                break;
            case OUTFLOW:
                p(pIEnd() - 1, j) = -p(uInteriorIEnd() - 1, j);
                break;
            default:
                break;
        }
    }
};

/**
 * set temperature at obstacle boundaries as well as at domain boundaries
 */
void StaggeredGrid::applyBoundaryTemperature() {
    for (int i = pIBegin(); i < pIEnd(); i++) {
        for (int j = pJBegin(); j < pJEnd(); j++) {
            switch (marker(i, j)) {
                case OBSTACLE:
                    t(i, j) = tobs(i, j);
                    break;
                case OBSTACLE_LEFT:
                    t(i, j) = 2.0 * tobs(i, j) - t(i - 1, j);
                    break;
                case OBSTACLE_RIGHT:
                    t(i, j) = 2.0 * tobs(i, j) - t(i + 1, j);
                    break;
                case OBSTACLE_TOP:
                    t(i, j) = 2.0 * tobs(i, j) - t(i, j + 1);
                    break;
                case OBSTACLE_BOTTOM:
                    t(i, j) = 2.0 * tobs(i, j) - t(i, j - 1);
                    break;
                case OBSTACLE_LEFT_TOP:
                    t(i, j) = 2.0 * tobs(i, j) - (t(i - 1, j) + t(i, j + 1)) / 2.0;
                    break;
                case OBSTACLE_RIGHT_TOP:
                    t(i, j) = 2.0 * tobs(i, j) - (t(i + 1, j) + t(i, j + 1)) / 2.0;
                    break;
                case OBSTACLE_LEFT_BOTTOM:
                    t(i, j) = 2.0 * tobs(i, j) - (t(i - 1, j) + t(i, j - 1)) / 2.0;
                    break;
                case OBSTACLE_RIGHT_BOTTOM:
                    t(i, j) = 2.0 * tobs(i, j) - (t(i + 1, j) + t(i, j - 1)) / 2.0;
                    break;
                default:
                    break;
            }
        }
    }

    // set boundary values for t at bottom and top side (lower priotity)
    for (int i = tIBegin(); i < tIEnd(); i++) {
        if (settings_.setFixedTempBottom) {
            t(i, tJBegin()) = 2.0 * settings_.tempBcBottom - t(i, tInteriorJBegin());
        } else {
            t(i, tJBegin()) = t(i, tInteriorJBegin()) - dy() * settings_.tempBcBottom;
        }
        if (settings_.setFixedTempTop) {
            t(i, tJEnd() - 1) = 2.0 * settings_.tempBcTop - t(i, tInteriorJEnd() - 1);
        } else {
            t(i, tJEnd() - 1) = t(i, tInteriorJEnd() - 1) - dy() * settings_.tempBcTop;
        }
    }

    // set boundary values for t at left and right side (higher priority)
    for (int j = tJBegin(); j < tJEnd(); j++) {
        if (settings_.setFixedTempLeft) {
            t(tIBegin(), j) = 2.0 * settings_.tempBcLeft - t(tInteriorIBegin(), j);
        } else {
            t(tIBegin(), j) = t(tInteriorIBegin(), j) - dx() * settings_.tempBcLeft;
        }
        if (settings_.setFixedTempRight) {
            t(tIEnd() - 1, j) = 2.0 * settings_.tempBcRight - t(tInteriorIEnd() - 1, j);
        } else {
            t(tIEnd() - 1, j) = t(tInteriorIEnd() - 1, j) - dx() * settings_.tempBcRight;
        }
    }
};

/**
 * free surface
 */

/**
 * set velocites (u, v) and pressure p at surface between fluid and free cells
 * must differentiate between many cases which require different handling
 * @param dt: time step width
 */
void StaggeredGrid::setSurfaceValues(double dt) {
    // TODO: check order in which values are calculated! (e.g. gray depends on green)
    for (int i = pIBegin(); i < pIEnd(); i++) {
        for (int j = pJBegin(); j < pJEnd(); j++) {
            switch (marker(i, j)) {
                case SURFACE_RIGHT: // nx=1, ny=0, mx=0, my=-1
                    // tangential: du/dy + dv/dx = 0
                    v(i + 1, j - 1) = v(i, j - 1) - dx() / dy() * (u(i, j) - u(i, j - 1));
                    // normal: p - 2/Re * du/dx = 0
                    p(i, j) = 2.0 / settings_.re * (u(i, j) - u(i - 1, j)) / dx();
                    // continuity: du/dx + dv/dy = 0
                    u(i, j) = u(i - 1, j) - dx() / dy() * (v(i, j) - v(i, j - 1));
                    break;
                case SURFACE_LEFT: // nx=-1, ny=0, mx=0, my=1
                    // tangential: du/dy + dv/dx = 0
                    v(i - 1, j - 1) = v(i, j - 1) - dx() / dy() * (u(i - 1, j) - u(i - 1, j - 1));
                    // normal: p - 2/Re * du/dx = 0
                    p(i, j) = 2.0 / settings_.re * (u(i, j) - u(i - 1, j)) / dx();
                    // continuity: du/dx + dv/dy = 0
                    u(i - 1, j) = u(i, j) + dx() / dy() * (v(i, j) - v(i, j - 1));
                    break;
                case SURFACE_TOP: // nx=0, ny=1, mx=1, my=0
                    // tangential: du/dy + dv/dx = 0
                    u(i - 1, j) = u(i - 1, j - 1) - dy() / dx() * (v(i, j) - v(i - 1, j));
                    // normal: p - 2/Re * dv/dy = 0
                    p(i, j) = 2.0 / settings_.re * (v(i, j) - v(i, j - 1)) / dy();
                    // continuity: du/dx + dv/dy = 0
                    v(i, j) = v(i - 1, j) - dy() / dx() * (u(i, j) - u(i, j - 1));
                    break;
                case SURFACE_BOTTOM: // nx=0, ny=-1, mx=-1, my=0
                    // tangential: du/dy + dv/dx = 0
                    u(i - 1, j - 1) = u(i - 1, j) + dy() / dx() * (v(i, j - 1) - v(i - 1, j - 1));
                    // normal: p - 2/Re * dv/dy = 0
                    p(i, j) = 2.0 / settings_.re * (v(i, j) - v(i, j - 1)) / dy();
                    // continuity: du/dx + dv/dy = 0
                    v(i, j - 1) = v(i, j) + dy() / dx() * (u(i, j) - u(i - 1, j));
                    break;
                    // ------------------------
                case SURFACE_RIGHT_TOP: // nx=1/sqrt(2), ny=1/sqrt(2), mx=1/sqrt(2), my=-1/sqrt(2)
                    // tangential:
                    u(i - 1, j + 1) = u(i - 1, j) - dy() / dx() * (v(i, j) - v(i - 1, j));
                    v(i + 1, j - 1) = v(i, j - 1) - dx() / dy() * (u(i, j) - u(i, j - 1));
                    // continuity + tangential (first, red):
                    u(i, j) = u(i - 1, j);
                    v(i, j) = v(i, j - 1);
                    // continuity + tangential (second, gray):
                    u(i, j + 1) = u(i, j);
                    v(i + 1, j) = v(i, j);
                    // normal:
                    p(i, j) = 1.0 / settings_.re *
                              ((u(i, j) + u(i - 1, j) - u(i, j - 1) - u(i - 1, j - 1)) / (2.0 * dx()) +
                               (v(i, j) + v(i, j - 1) - v(i - 1, j) - v(i - 1, j - 1)) / (2.0 * dy()));
                    break;
                case SURFACE_RIGHT_BOTTOM: // nx=1/sqrt(2), ny=-1/sqrt(2), mx=-1/sqrt(2), my=-1/sqrt(2)
                    // tangential:
                    u(i - 1, j - 1) = u(i - 1, j) - dy() / dx() * (v(i, j - 1) - v(i - 1, j - 1));
                    // continuity + tangential (first, red):
                    u(i, j) = u(i - 1, j);
                    v(i, j - 1) = v(i, j);
                    // continuity + tangential (second, gray):
                    u(i, j - 1) = u(i, j);
                    v(i + 1, j - 1) = v(i, j - 1);
                    // normal:
                    p(i, j) = 1.0 / settings_.re *
                              ((u(i, j) + u(i, j + 1) - u(i - 1, j) - u(i - 1, j + 1)) / (2.0 * dx()) +
                               (v(i, j) + v(i, j - 1) - v(i - 1, j) - v(i - 1, j - 1)) / (2.0 * dy()));
                    break;
                case SURFACE_LEFT_BOTTOM: // nx=-1/sqrt(2), ny=-1/sqrt(2), mx=-1/sqrt(2), my=1/sqrt(2)
                    // tangential:
                    // nothing :)
                    // continuity + tangential (first, red):
                    u(i - 1, j) = u(i, j);
                    v(i, j - 1) = v(i, j);
                    // continuity + tangential (second, gray):
                    u(i - 1, j - 1) = u(i - 1, j);
                    v(i - 1, j - 1) = v(i, j - 1);
                    // normal:
                    p(i, j) = 1.0 / settings_.re *
                              ((u(i, j) + u(i, j + 1) - u(i - 1, j) - u(i - 1, j + 1)) / (2.0 * dx()) +
                               (v(i, j) + v(i + 1, j) - v(i, j - 1) - v(i + 1, j - 1)) / (2.0 * dy()));
                    break;
                case SURFACE_LEFT_TOP: // nx=-1/sqrt(2), ny=1/sqrt(2), mx=1/sqrt(2), my=1/sqrt(2)
                    // tangential:
                    v(i - 1, j - 1) = v(i, j - 1) + dx() / dy() * (u(i - 1, j) - u(i - 1, j - 1));
                    // continuity + tangential (first, red):
                    u(i - 1, j) = u(i, j);
                    v(i, j) = v(i, j - 1);
                    // continuity + tangential (second, gray):
                    u(i - 1, j + 1) = u(i - 1, j);
                    v(i - 1, j - 1) = v(i, j - 1);
                    // normal:
                    p(i, j) = 1.0 / settings_.re *
                              ((u(i, j) + u(i, j - 1) - u(i - 1, j) - u(i - 1, j - 1)) / (2.0 * dx()) +
                               (v(i, j) + v(i + 1, j) - v(i, j - 1) - v(i + 1, j - 1)) / (2.0 * dy()));
                    // ------------------------
                case SURFACE_LEFT_RIGHT:
                    // tangential:
                    v(i - 1, j - 1) = v(i, j - 1) + dx() / dy() * (u(i - 1, j) - u(i - 1, j - 1));
                    v(i + 1, j - 1) = v(i, j - 1) - dx() / dy() * (u(i, j) - u(i, j - 1));
                    // driven by gravity:
                    u(i, j) = u(i, j) + dt * settings_.g[0];
                    u(i - 1, j) = u(i - 1, j) + dt * settings_.g[0];
                    // normal:
                    p(i, j) = 2.0 / settings_.re * ((u(i, j) - u(i - 1, j)) / dx());
                case SURFACE_TOP_BOTTOM:
                    // tangential:
                    u(i - 1, j + 1) = u(i - 1, j) - dy() / dx() * (v(i, j) - v(i - 1, j));
                    u(i - 1, j - 1) = u(i - 1, j) + dy() / dx() * (v(i, j - 1) - v(i - 1, j - 1));
                    // driven by gravity:
                    v(i, j) = v(i, j) + dt * settings_.g[1];
                    v(i, j - 1) = v(i, j - 1) + dt * settings_.g[1];
                    // normal:
                    p(i, j) = 2.0 / settings_.re * ((v(i, j) - v(i, j - 1)) / dy());
                    break;
                    // ------------------------
                case SURFACE_LEFT_TOP_RIGHT:
                    // tangential:
                    v(i - 1, j - 1) = v(i, j - 1) + dx() / dy() * (u(i - 1, j) - u(i - 1, j - 1));
                    v(i + 1, j - 1) = v(i, j - 1) - dx() / dy() * (u(i, j) - u(i, j - 1));
                    // driven by gravity:
                    u(i, j) = u(i, j) + dt * settings_.g[0];
                    u(i - 1, j) = u(i - 1, j) + dt * settings_.g[0];
                    // continuity + tangential (second, gray):
                    u(i, j + 1) = u(i, j);
                    v(i + 1, j) = v(i, j);
                    // continuity + tangential (second, gray):
                    u(i - 1, j + 1) = u(i - 1, j);
                    v(i - 1, j - 1) = v(i, j - 1);
                    // continuity (blue)
                    v(i, j) = v(i, j - 1) - dy() / dx() * (u(i, j) - u(i - 1, j));
                    // normal:
                    p(i, j) = 2.0 / settings_.re * ((u(i, j) - u(i - 1, j)) / dx());
                    break;
                case SURFACE_TOP_RIGHT_BOTTOM:
                    // tangential:
                    u(i - 1, j + 1) = u(i - 1, j) - dy() / dx() * (v(i, j) - v(i - 1, j));
                    u(i - 1, j - 1) = u(i - 1, j) + dy() / dx() * (v(i, j - 1) - v(i - 1, j - 1));
                    // driven by gravity:
                    v(i, j) = v(i, j) + dt * settings_.g[1];
                    v(i, j - 1) = v(i, j - 1) + dt * settings_.g[1];
                    // continuity + tangential (second, gray):
                    u(i, j + 1) = u(i, j);
                    v(i + 1, j) = v(i, j);
                    // continuity + tangential (second, gray):
                    u(i, j - 1) = u(i, j);
                    v(i + 1, j - 1) = v(i, j - 1);
                    // continuity (blue)
                    u(i, j) = u(i - 1, j) - dx() / dy() + (v(i, j) - v(i, j - 1));
                    // normal:
                    p(i, j) = 2.0 / settings_.re * ((v(i, j) - v(i, j - 1)) / dy());
                    break;
                case SURFACE_RIGHT_BOTTOM_LEFT:
                    // driven by gravity:
                    u(i, j) = u(i, j) + dt * settings_.g[0];
                    u(i - 1, j) = u(i - 1, j) + dt * settings_.g[0];
                    // continuity + tangential (second, gray):
                    u(i - 1, j - 1) = u(i - 1, j);
                    v(i - 1, j - 1) = v(i, j - 1);
                    // continuity + tangential (second, gray):
                    u(i, j - 1) = u(i, j);
                    v(i + 1, j - 1) = v(i, j - 1);
                    // continuity (blue)
                    v(i, j - 1) = v(i, j) + dy() / dx() * (u(i, j) - u(i - 1, j));
                    // normal:
                    p(i, j) = 2.0 / settings_.re * ((u(i, j) - u(i - 1, j)) / dx());
                    break;
                case SURFACE_BOTTOM_LEFT_TOP:
                    // driven by gravity:
                    v(i, j) = v(i, j) + dt * settings_.g[1];
                    v(i, j - 1) = v(i, j - 1) + dt * settings_.g[1];
                    // continuity + tangential (second, gray):
                    u(i - 1, j + 1) = u(i - 1, j);
                    v(i - 1, j - 1) = v(i, j - 1);
                    // continuity + tangential (second, gray):
                    u(i - 1, j - 1) = u(i - 1, j);
                    v(i - 1, j - 1) = v(i, j - 1);
                    // continuity (blue)
                    u(i - 1, j) = u(i, j) + dx() / dy() + (v(i, j) - v(i, j - 1));
                    // normal:
                    p(i, j) = 2.0 / settings_.re * ((v(i, j) - v(i, j - 1)) / dy());
                    break;
                    // ------------------------
                case SURFACE_TOP_RIGHT_BOTTOM_LEFT:
                    // driven by gravity:
                    u(i, j) = u(i, j) + dt * settings_.g[0];
                    u(i - 1, j) = u(i - 1, j) + dt * settings_.g[0];
                    // driven by gravity:
                    v(i, j) = v(i, j) + dt * settings_.g[1];
                    v(i, j - 1) = v(i, j - 1) + dt * settings_.g[1];
                    // continuity + tangential (second, gray):
                    u(i - 1, j + 1) = u(i - 1, j);
                    v(i - 1, j - 1) = v(i, j - 1);
                    // continuity + tangential (second, gray):
                    u(i - 1, j - 1) = u(i - 1, j);
                    v(i - 1, j - 1) = v(i, j - 1);
                    // continuity + tangential (second, gray):
                    u(i, j + 1) = u(i, j);
                    v(i + 1, j) = v(i, j);
                    // continuity + tangential (second, gray):
                    u(i, j - 1) = u(i, j);
                    v(i + 1, j - 1) = v(i, j - 1);
                    // normal:
                    p(i, j) = 2.0 / settings_.re * ((v(i, j) - v(i, j - 1)) / dy());
                    break;
                default:
                    break;
            }
        }
    }
}

/**
 * mesh information
 */

/**
 * get the mesh width
 * @return mesh width
 */
const std::array<double, 2> StaggeredGrid::meshWidth() const {
    return meshWidth_;
};

/**
 * get number of cells
 * @return number of cells
 */
const std::array<int, 2> StaggeredGrid::nCells() const {
    return nCells_;
};

/**
 * get mesh width in x-direction
 * @return mesh width in x-direction
 */
double StaggeredGrid::dx() const {
    return meshWidth_[0];
};

/**
 * get mesh width in y-direction
 * @return mesh width in y-direction
 */
double StaggeredGrid::dy() const {
    return meshWidth_[1];
};

/**
 * pressure variable p
 */

/**
 * get first valid index for p in x direction
 * @return first valid index for p in x direction
 */
int StaggeredGrid::pIBegin() const {
    return -1;
};

/**
 * get one after last valid index for p in x direction
 * @return one after last valid index for p in x direction
 */
int StaggeredGrid::pIEnd() const {
    return nCells_[0] + 1;
};

/**
 * get first valid index for p in y direction
 * @return first valid index for p in y direction
 */
int StaggeredGrid::pJBegin() const {
    return -1;
};

/**
 * get one after last valid index for p in y direction
 * @return one after last valid index for p in y direction
 */
int StaggeredGrid::pJEnd() const {
    return nCells_[1] + 1;
};

/**
 * get size of FieldVariable p
 * @return size of FieldVariable p
 */
std::array<int, 2> StaggeredGrid::pSize() const {
    return {pIEnd() - pIBegin(), pJEnd() - pJBegin()};
}

/**
 * get first valid Interior index for p in x direction
 * @return first valid Interior index for p in x direction
 */
int StaggeredGrid::pInteriorIBegin() const {
    return pIBegin() + 1;
};

/**
 * get one after last valid Interior index for p in x direction
 * @return one after last valid Interior index for p in x direction
 */
int StaggeredGrid::pInteriorIEnd() const {
    return pIEnd() - 1;
};

/**
 * get first valid Interior index for p in y direction
 * @return first valid Interior index for p in y direction
 */
int StaggeredGrid::pInteriorJBegin() const {
    return pJBegin() + 1;
};

/**
 * get first valid Interior index for p in y direction
 * @return first valid Interior index for p in y direction
 */
int StaggeredGrid::pInteriorJEnd() const {
    return pJEnd() - 1;
};

/**
 * evaluate field variable p in an element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return field variable p in an element (i,j)
 */
MARKER StaggeredGrid::marker(int i, int j) const {
#ifndef NDEBUG
    assert((pIBegin() <= i) && (i <= pIEnd()));
    assert((pJBegin() <= j) && (j <= pJEnd()));
#endif
    return (MARKER) marker_(i - pIBegin(), j - pJBegin());
};

/**
 * evaluate field variable p in an element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return field variable p in an element (i,j)
 */
MARKER &StaggeredGrid::marker(int i, int j) {
#ifndef NDEBUG
    assert((pIBegin() <= i) && (i <= pIEnd()));
    assert((pJBegin() <= j) && (j <= pJEnd()));
#endif
    return (MARKER &) marker_(i - pIBegin(), j - pJBegin());
};

/**
 * get reference to field variable p
 * @return reference to field variable p
 */
const FieldVariable &StaggeredGrid::p() const {
    return p_;
};

/**
 * evaluate field variable p in an element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return field variable p in an element (i,j)
 */
double StaggeredGrid::p(int i, int j) const {
#ifndef NDEBUG
    assert((pIBegin() <= i) && (i <= pIEnd()));
    assert((pJBegin() <= j) && (j <= pJEnd()));
#endif
    return p_(i - pIBegin(), j - pJBegin());
};

/**
 * evaluate field variable p in an element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return field variable p in an element (i,j)
 */
double &StaggeredGrid::p(int i, int j) {
#ifndef NDEBUG
    assert((pIBegin() <= i) && (i <= pIEnd()));
    assert((pJBegin() <= j) && (j <= pJEnd()));
#endif
    return p_(i - pIBegin(), j - pJBegin());
};

/**
 * velocity in x-direction u
 */

/**
 * get first valid index for u in x direction
 * @return first valid index for u in x direction
 */
int StaggeredGrid::uIBegin() const {
    if (partitioning_->ownPartitionContainsLeftBoundary()) {
        return -1;
    }
    return -2;
};

/**
 * get one after last valid index for u in x direction
 * @return one after last valid index for u in x direction
 */
int StaggeredGrid::uIEnd() const {
    if (partitioning_->ownPartitionContainsRightBoundary()) {
        return nCells_[0];
    }
    return nCells_[0] + 1;
};

/**
 * get first valid index for u in y direction
 * @return first valid index for u in y direction
 */
int StaggeredGrid::uJBegin() const {
    return -1;
};

/**
 * get one after last valid index for u in y direction
 * @return one after last valid index for u in y direction
 */
int StaggeredGrid::uJEnd() const {
    return nCells_[1] + 1;
};

/**
 * get size of FieldVariable u
 * @return size of FieldVariable u
 */
std::array<int, 2> StaggeredGrid::uSize() const {
    return {uIEnd() - uIBegin(), uJEnd() - uJBegin()};
}

/**
 * get first valid field index for u in x direction
 * @return first valid field index for u in x direction
 */
int StaggeredGrid::uInteriorIBegin() const {

    return uIBegin() + 1;
};

/**
 * get first valid field index for u in x direction
 * @return first valid field index for u in x direction
 */
int StaggeredGrid::uInteriorIEnd() const {

    return uIEnd() - 1;

};

/**
 * first valid Interior index for u in y direction
 * @return first inner grid value of u in j direction
 */
int StaggeredGrid::uInteriorJBegin() const {

    return uJBegin() + 1;

};

/**
 * get one after last valid Interior index for u in y direction
 * @return one after last valid Interior index for u in y direction
 */
int StaggeredGrid::uInteriorJEnd() const {

    return uJEnd() - 1;

};

/**
 * get a reference to field variable u
 * @return reference to field variable u
 */
const FieldVariable &StaggeredGrid::u() const {
    return u_;
};

/**
 * access value of u in element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return value of u in element (i,j)
 */
double StaggeredGrid::u(int i, int j) const {
#ifndef NDEBUG
    assert((uIBegin() <= i) && (i <= uIEnd()));
    assert((uJBegin() <= j) && (j <= uJEnd()));
#endif
    return u_(i - uIBegin(), j - uJBegin());
};

/**
 * access value of u in element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return value of u in element (i,j)
 */
double &StaggeredGrid::u(int i, int j) {
#ifndef NDEBUG
    assert((uIBegin() <= i) && (i <= uIEnd()));
    assert((uJBegin() <= j) && (j <= uJEnd()));
#endif
    return u_(i - uIBegin(), j - uJBegin());
};

/**
 * get a reference to field variable u
 * @return reference to field variable u
 */
const FieldVariable &StaggeredGrid::uLast() const {
    return uLast_;
};

/**
 * access value of u in element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return value of u in element (i,j)
 */
double StaggeredGrid::uLast(int i, int j) const {
#ifndef NDEBUG
    assert((uIBegin() <= i) && (i <= uIEnd()));
    assert((uJBegin() <= j) && (j <= uJEnd()));
#endif
    return uLast_(i - uIBegin(), j - uJBegin());
};

/**
 * access value of u in element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return value of u in element (i,j)
 */
double &StaggeredGrid::uLast(int i, int j) {
#ifndef NDEBUG
    assert((uIBegin() <= i) && (i <= uIEnd()));
    assert((uJBegin() <= j) && (j <= uJEnd()));
#endif
    return uLast_(i - uIBegin(), j - uJBegin());
};

/**
 * velocity in y-direction v
 */

/**
 * get first valid index for v in x direction
 * @return first valid index for v in x direction
 */
int StaggeredGrid::vIBegin() const {
    return -1;
};

/**
 * get one after last valid index for v in x direction
 * @return one after last valid index for v in x direction
 */
int StaggeredGrid::vIEnd() const {
    return nCells_[0] + 1;
};

/**
 * get first valid index for v in y direction
 * @return first valid index for v in y direction
 */
int StaggeredGrid::vJBegin() const {
    if (partitioning_->ownPartitionContainsBottomBoundary()) {
        return -1;
    }
    return -2;
};

/**
 * get one after last valid index for v in y direction
 * @return one after last valid index for v in y direction
 */
int StaggeredGrid::vJEnd() const {
    if (partitioning_->ownPartitionContainsTopBoundary()) {
        return nCells_[1];
    }
    return nCells_[1] + 1;
};

/**
 * get size of FieldVariable v
 * @return size of FieldVariable v
 */
std::array<int, 2> StaggeredGrid::vSize() const {
    return {vIEnd() - vIBegin(), vJEnd() - vJBegin()};
}

/**
 * get first valid Interior index for v in x direction
 * @return first valid Interior index for v in x direction
 */
int StaggeredGrid::vInteriorIBegin() const {
    return vIBegin() + 1;

};


/**
 * get one after last valid Interior index for v in x direction
 * @return one after last valid Interior index for v in x direction
 */
int StaggeredGrid::vInteriorIEnd() const {
    return vIEnd() - 1;

};

/**
 * get first valid Interior index for v in y direction
 * @return first valid Interior index for v in y direction
 */
int StaggeredGrid::vInteriorJBegin() const {
    return vJBegin() + 1;

};

/**
 * get one after last valid Interior index for v in y direction
 * @return one after last valid Interior index for v in y direction
 */
int StaggeredGrid::vInteriorJEnd() const {
    return vJEnd() - 1;

};

/**
 * get a reference to field variable v
 * @return a reference to field variable v
 */
const FieldVariable &StaggeredGrid::v() const {
    return v_;
};

/**
 * access value of v in element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return value of v in element (i,j)
 */
double StaggeredGrid::v(int i, int j) const {
#ifndef NDEBUG
    assert((vIBegin() <= i) && (i <= vIEnd()));
    assert((vJBegin() <= j) && (j <= vJEnd()));
#endif
    return v_(i - vIBegin(), j - vJBegin());
};

/**
 * access value of v in element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return value of v in element (i,j)
 */
double &StaggeredGrid::v(int i, int j) {
#ifndef NDEBUG
    assert((vIBegin() <= i) && (i <= vIEnd()));
    assert((vJBegin() <= j) && (j <= vJEnd()));
#endif
    return v_(i - vIBegin(), j - vJBegin());
};

/**
 * get a reference to field variable v
 * @return a reference to field variable v
 */
const FieldVariable &StaggeredGrid::vLast() const {
    return vLast_;
};

/**
 * access value of v in element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return value of v in element (i,j)
 */
double StaggeredGrid::vLast(int i, int j) const {
#ifndef NDEBUG
    assert((vIBegin() <= i) && (i <= vIEnd()));
    assert((vJBegin() <= j) && (j <= vJEnd()));
#endif
    return vLast_(i - vIBegin(), j - vJBegin());
};

/**
 * access value of v in element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return value of v in element (i,j)
 */
double &StaggeredGrid::vLast(int i, int j) {
#ifndef NDEBUG
    assert((vIBegin() <= i) && (i <= vIEnd()));
    assert((vJBegin() <= j) && (j <= vJEnd()));
#endif
    return vLast_(i - vIBegin(), j - vJBegin());
};

/**
 * right hand side rhs
 */

/**
 * get first valid index for rhs in x direction
 * @return first valid index for rhs in x direction
 */
int StaggeredGrid::rhsIBegin() const {
    return -1;
};

/**
 * get one after last valid index for rhs in x direction
 * @return one after last valid index for rhs in x direction
 */
int StaggeredGrid::rhsIEnd() const {
    return nCells_[0] + 1;
};


/**
 * get first valid index for rhs in y direction
 * @return first valid index for rhs in y direction
 */
int StaggeredGrid::rhsJBegin() const {
    return -1;
};

/**
 * get one after last valid index for rhs in y direction
 * @return one after last valid index for rhs in y direction
 */
int StaggeredGrid::rhsJEnd() const {
    return nCells_[1] + 1;
};

/**
 * get size of FieldVariable rhs
 * @return size of FieldVariable rhs
 */
std::array<int, 2> StaggeredGrid::rhsSize() const {
    return {rhsIEnd() - rhsIBegin(), rhsJEnd() - rhsJBegin()};
}

/**
 * get first valid field index for u in x direction
 * @return first valid field index for u in x direction
 */
int StaggeredGrid::rhsInteriorIBegin() const {

    return rhsIBegin() + 1;
};

/**
 * get one after last valid Interior index for u in x direction
 * @return one after last valid Interior index for u in x direction
 */
int StaggeredGrid::rhsInteriorIEnd() const {

    return rhsIEnd() - 1;

};

/**
 * get first valid Interior index for u in y direction
 * @return first valid Interior index for u in y direction
 */
int StaggeredGrid::rhsInteriorJBegin() const {

    return rhsJBegin() + 1;

};

/**
 * get one after last valid Interior index for u in y direction
 * @return one after last valid Interior index for u in y direction
 */
int StaggeredGrid::rhsInteriorJEnd() const {

    return rhsJEnd() - 1;

};

/**
 * access value of rhs in element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return value of rhs in element (i,j)
 */
double StaggeredGrid::rhs(int i, int j) const {
#ifndef NDEBUG
    assert((rhsIBegin() <= i) && (i <= rhsIEnd()));
    assert((rhsJBegin() <= j) && (j <= rhsJEnd()));
#endif
    return rhs_(i - rhsIBegin(), j - rhsJBegin());
};

/**
 * access value of rhs in element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return value of rhs in element (i,j)
 */
double &StaggeredGrid::rhs(int i, int j) {
#ifndef NDEBUG
    assert((rhsIBegin() <= i) && (i <= rhsIEnd()));
    assert((rhsJBegin() <= j) && (j <= rhsJEnd()));
#endif
    return rhs_(i - rhsIBegin(), j - rhsJBegin());
};

/**
 * get reference to field variable rhs
 * @return reference to field variable rhs
 */
const FieldVariable &StaggeredGrid::rhs() const {
    return rhs_;
};

/**
 * preliminary velocity F
 */

/**
 * access value of F in element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return value of F in element (i,j)
 */
double StaggeredGrid::f(int i, int j) const {
#ifndef NDEBUG
    assert((uIBegin() <= i) && (i <= uIEnd()));
    assert((uJBegin() <= j) && (j <= uJEnd()));
#endif
    return f_(i - uIBegin(), j - uJBegin());
};

/**
 * access value of F in element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return value of F in element (i,j)
 */
double &StaggeredGrid::f(int i, int j) {
#ifndef NDEBUG
    assert((uIBegin() <= i) && (i <= uIEnd()));
    assert((uJBegin() <= j) && (j <= uJEnd()));
#endif
    return f_(i - uIBegin(), j - uJBegin());
};


/**
 * get reference to field variable F
 * @return reference to field variable F
 */
const FieldVariable &StaggeredGrid::f() const {
    return f_;
};

/**
 * preliminary velocity G
 */


/**
 * access value of G in element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return value of G in element (i,j)
 */
double StaggeredGrid::g(int i, int j) const {
#ifndef NDEBUG
    assert((vIBegin() <= i) && (i <= vIEnd()));
    assert((vJBegin() <= j) && (j <= vJEnd()));
#endif
    return g_(i - vIBegin(), j - vJBegin());
};

/**
 * access value of G in element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return value of G in element (i,j)
 */
double &StaggeredGrid::g(int i, int j) {
#ifndef NDEBUG
    assert((vIBegin() <= i) && (i <= vIEnd()));
    assert((vJBegin() <= j) && (j <= vJEnd()));
#endif
    return g_(i - vIBegin(), j - vJBegin());
};

/**
 * get reference to field variable G
 * @return reference to field variable G
 */
const FieldVariable &StaggeredGrid::g() const {
    return g_;
};

/**
 * temperature variable t
 */

/**
 * get first valid index for t in x direction
 * @return first valid index for t in x direction
 */
int StaggeredGrid::tIBegin() const {
    return -1;
};

/**
 * get one after last valid index for t in x direction
 * @return one after last valid index for t in x direction
 */
int StaggeredGrid::tIEnd() const {
    return nCells_[0] + 1;
};

/**
 * get first valid index for t in y direction
 * @return first valid index for t in y direction
 */
int StaggeredGrid::tJBegin() const {
    return -1;
};

/**
 * get one after last valid index for t in y direction
 * @return one after last valid index for t in y direction
 */
int StaggeredGrid::tJEnd() const {
    return nCells_[1] + 1;
};

/**
 * get size of FieldVariable t
 * @return size of FieldVariable t
 */
std::array<int, 2> StaggeredGrid::tSize() const {
    return {tIEnd() - tIBegin(), tJEnd() - tJBegin()};
}

/**
 * get first valid Interior index for t in x direction
 * @return first valid Interior index for t in x direction
 */
int StaggeredGrid::tInteriorIBegin() const {
    return tIBegin() + 1;
};

/**
 * get one after last valid Interior index for t in x direction
 * @return one after last valid Interior index for t in x direction
 */
int StaggeredGrid::tInteriorIEnd() const {
    return tIEnd() - 1;
};

/**
 * get first valid Interior index for t in y direction
 * @return first valid Interior index for t in y direction
 */
int StaggeredGrid::tInteriorJBegin() const {
    return tJBegin() + 1;
};

/**
 * get first valid Interior index for t in y direction
 * @return first valid Interior index for t in y direction
 */
int StaggeredGrid::tInteriorJEnd() const {
    return tJEnd() - 1;
};

/**
 * get reference to field variable t
 * @return reference to field variable t
 */
const FieldVariable &StaggeredGrid::t() const {
    return t_;
};

/**
 * evaluate field variable t in an element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return field variable t in an element (i,j)
 */
double StaggeredGrid::t(int i, int j) const {
#ifndef NDEBUG
    assert((tIBegin() <= i) && (i <= tIEnd()) && "temperature i failed in const");
    assert((tJBegin() <= j) && (j <= tJEnd()) && "temperature j failed in const");
#endif
    return t_(i - tIBegin(), j - tJBegin());
};

/**
 * evaluate field variable t in an element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return field variable t in an element (i,j)
 */
double &StaggeredGrid::t(int i, int j) {
#ifndef NDEBUG
    assert((tIBegin() <= i) && (i <= tIEnd()) && "temperature i failed");
    assert((tJBegin() <= j) && (j <= tJEnd()) && "temperature j failed");
#endif
    return t_(i - tIBegin(), j - tJBegin());
};

/**
 * get reference to field variable t
 * @return reference to field variable t
 */
const FieldVariable &StaggeredGrid::tobs() const {
    return tobs_;
};

/**
 * evaluate field variable t in an element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return field variable t in an element (i,j)
 */
double StaggeredGrid::tobs(int i, int j) const {
#ifndef NDEBUG
    assert((tIBegin() <= i) && (i <= tIEnd()) && "temperature i failed in const");
    assert((tJBegin() <= j) && (j <= tJEnd()) && "temperature j failed in const");
#endif
    return tobs_(i - tIBegin(), j - tJBegin());
};

/**
 * evaluate field variable t in an element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return field variable t in an element (i,j)
 */
double &StaggeredGrid::tobs(int i, int j) {
#ifndef NDEBUG
    assert((tIBegin() <= i) && (i <= tIEnd()) && "temperature i failed");
    assert((tJBegin() <= j) && (j <= tJEnd()) && "temperature j failed");
#endif
    return tobs_(i - tIBegin(), j - tJBegin());
};

/**
 * sources and sinks for the heat equation q
 */

/**
 * get first valid index for q in x direction
 * @return first valid index for q in x direction
 */
int StaggeredGrid::qIBegin() const {
    return -1;
};

/**
 * get one after last valid index for q in x direction
 * @return one after last valid index for q in x direction
 */
int StaggeredGrid::qIEnd() const {
    return nCells_[0] + 1;
};

/**
 * get first valid index for q in y direction
 * @return first valid index for q in y direction
 */
int StaggeredGrid::qJBegin() const {
    return -1;
};

/**
 * get one after last valid index for q in y direction
 * @return one after last valid index for q in y direction
 */
int StaggeredGrid::qJEnd() const {
    return nCells_[1] + 1;
};

/**
 * get size of FieldVariable t
 * @return size of FieldVariable t
 */
std::array<int, 2> StaggeredGrid::qSize() const {
    return {qIEnd() - qIBegin(), qJEnd() - qJBegin()};
}

/**
 * get first valid Interior index for q in x direction
 * @return first valid Interior index for q in x direction
 */
int StaggeredGrid::qInteriorIBegin() const {
    return qIBegin() + 1;
};

/**
 * get one after last valid Interior index for q in x direction
 * @return one after last valid Interior index for q in x direction
 */
int StaggeredGrid::qInteriorIEnd() const {
    return qIEnd() - 1;
};

/**
 * get first valid Interior index for q in y direction
 * @return first valid Interior index for q in y direction
 */
int StaggeredGrid::qInteriorJBegin() const {
    return qJBegin() + 1;
};

/**
 * get first valid Interior index for q in y direction
 * @return first valid Interior index for q in y direction
 */
int StaggeredGrid::qInteriorJEnd() const {
    return qJEnd() - 1;
};

/**
 * get reference to field variable q
 * @return reference to field variable q
 */
const FieldVariable &StaggeredGrid::q() const {
    return q_;
};

/**
 * evaluate field variable temperature souce term q in an element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return field variable q in an element (i,j)
 */
double StaggeredGrid::q(int i, int j) const {
#ifndef NDEBUG
    assert((qIBegin() <= i) && (i <= qIEnd()) && "Q i failed in const");
    assert((qJBegin() <= j) && (j <= qJEnd()) && "Q j failed in const");
#endif
    return q_(i - qIBegin(), j - qJBegin());
};

/**
 * evaluate field variable temperature souce term q in an element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return field variable q in an element (i,j)
 */
double &StaggeredGrid::q(int i, int j) {
#ifndef NDEBUG
    assert((qIBegin() <= i) && (i <= qIEnd()) && "Q i failed");
    assert((qJBegin() <= j) && (j <= qJEnd()) && "Q j failed");
#endif
    return q_(i - qIBegin(), j - qJBegin());
};

/**
 * free surface
 */

/**
 * get number of particles
 * @return number of particles
 */
int StaggeredGrid::particleNumber() const {
    return particles_.number();
};
/**
 * get list of particles
 * @return particles
 */
const Particle2D &StaggeredGrid::particle() const {

    return particles_;
}

/**
 * get particle position in x direction
 * @param k: index of the particle
 * @return particle position in x direction
 */

double StaggeredGrid::particlePosX(int k) const{

    return particles_(k,0);
};

/**
 * get particle position in x direction
 * @param k: index of the particle
 * @return particle position in x direction
 */
double &StaggeredGrid::particlePosX(int k){
    return particles_(k,0);
};

/**
 * get particle position in y direction
 * @param k: index of the particle
 * @return particle position in y direction
 */
double StaggeredGrid::particlePosY(int k) const{

    return particles_(k,1);
};

/**
 * get particle position in y direction
 * @param k: index of the particle
 * @return particle position in y direction
 */
double &StaggeredGrid::particlePosY(int k) {

    return particles_(k, 1);
};

/**
 * split particle index k into indices i and j for x and y direction, respectively.
 * @param k: index of the particle
 * @return tuple of i and j for the particle position on a two dimensional grid
 */
std::array<int, 2> StaggeredGrid::particleCell(int k) const {
#ifndef NDEBUG
    assert((0 <= k) && (k <= particleNumber()) && "particle i failed in const");
#endif
    int i = (int) (particlePosX(k) / dx());
    int j = (int) (particlePosY(k) / dy());

    return {i, j};
}

/**
 * set fluid, free and surface cells accordingly. See lecture notes for a detailed description of the surface cases.
 */
void StaggeredGrid::updateCellTypes() {
    // Reset all fluid cells to free
    for (int i = pIBegin(); i < pIEnd(); i++) {
        for (int j = pJBegin(); j < pJEnd(); j++) {
            if (marker(i, j) == FLUID){
                marker(i, j) = FREE;
            }
        }
    }
    // Detect all cells containing any particles and define them as fluid
    for (int k = 0; k < particleNumber(); k++){
        std::array<int, 2> indices = particleCell(k);
        marker(indices[0], indices[1]) = MARKER::FLUID;
    }

    // Detect and mark the surface cells
    for (int i = pIBegin(); i < pIEnd(); i++) {
        for (int j = pJBegin(); j < pJEnd(); j++) {
            if (marker(i, j) != FLUID){
                continue;
            }
            int emptyCells = countFreeNeighbors(i, j);
            switch (emptyCells) {
                case 0:
                    marker(i, j) = FLUID;
                    break;
                case 1:
                    if (marker(i, j + 1) == FREE)
                        marker(i, j) = SURFACE_TOP;
                    else if (marker(i, j - 1) == FREE)
                        marker(i, j) = SURFACE_BOTTOM;
                    else if (marker(i + 1, j) == FREE)
                        marker(i, j) = SURFACE_RIGHT;
                    else if (marker(i - 1, j) == FREE)
                        marker(i, j) = SURFACE_LEFT;
                    break;
                case 2:
                    if (marker(i, j + 1) == FREE) {
                        if (marker(i + 1, j) == FREE) {
                            marker(i, j) = SURFACE_RIGHT_TOP;
                        }
                        else if (marker(i - 1, j) == FREE) {
                            marker(i, j) = SURFACE_LEFT_TOP;
                        }
                        else
                            marker(i, j) = SURFACE_TOP_BOTTOM;
                    }
                    else if (marker(i, j - 1) == FREE) {
                        if (marker(i + 1, j) == FREE) {
                            marker(i, j) = SURFACE_RIGHT_BOTTOM;
                        }
                        else if (marker(i - 1, j) == FREE) {
                            marker(i, j) = SURFACE_LEFT_BOTTOM;
                        }
                    }
                    else {
                        marker(i, j) = SURFACE_LEFT_RIGHT;
                    }
                    break;
                case 3:
                    if (marker(i, j + 1) == FLUID)
                        marker(i, j) = SURFACE_RIGHT_BOTTOM_LEFT;
                    else if (marker(i, j - 1) == FLUID)
                        marker(i, j) = SURFACE_LEFT_TOP_RIGHT;
                    else if (marker(i + 1, j) == FLUID)
                        marker(i, j) = SURFACE_TOP_RIGHT_BOTTOM;
                    else
                        marker(i, j) = SURFACE_BOTTOM_LEFT_TOP;
                    break;
                default: // 4
                    marker(i, j) = SURFACE_TOP_RIGHT_BOTTOM_LEFT;
                    break;
            }
        }
    }
};

/**
 * gets number of free neighbors a cell has. Diagonal neighbours are irrelevant.
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return number of free neighbor cells
 */
int StaggeredGrid::countFreeNeighbors(int i, int j) {
    int count = 0;
    if (marker(i, j + 1) == FREE)
        count++;
    if (marker(i, j - 1) == FREE)
        count++;
    if (marker(i + 1, j) == FREE)
        count++;
    if (marker(i - 1, j) == FREE)
        count++;
    return count;
};

/**
 * check if cell is inner fluid cell (has no free neighbor cell)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return boolean: true - is inner fluid cell or false (is not inner fluid cell)
 */
bool StaggeredGrid::isInnerFluid(int i, int j) {
    return (marker(i, j) == FLUID && countFreeNeighbors(i, j) == 0);
};