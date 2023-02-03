#pragma once

#include <vector>
#include <array>


/**
 * This class represents a 2D array of the positions of the particles.
 * Internally they are stored consecutively in memory.
 * The entries can be accessed by one index i.
 */
class Particle2D {
public:
    /**
    * constructor
    * @param size: number of cells
    */
    Particle2D(int number);

    /**
    * get number of particles
    * @return number of particles
    */
    int number() const;

    /**
    * access the position of the i-th particle, declared const, i.e. it is not possible to change the value
    * @param i: index of the particle
    * @return position in the grid
    */
    double &operator()(int k, int axis);

    /**
    * get the position of the i-th particle, declared const, i.e. it is not possible to change the value
    * @param i: index of the particle
    * @return position in the grid
    */
    double operator()(int k, int axis) const;

    void print() const;

    //! return data vector
    void *data();

    void setToZero();

protected:
    std::vector<double> data_;      //< storage array values, in row-major order
    const int number_; //< width, height of the domain
};
