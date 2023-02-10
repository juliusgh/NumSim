#pragma once

#include <vector>
#include <array>


/**
 * This class represents a 2D array of the positions of the particles.
 * Internally they are stored consecutively in memory.
 * The entries can be accessed by one index k.
 */
class Particle2D {
public:
    /**
    * constructor
    * @param number: number of particles
    */
    Particle2D(int number);

    /**
    * get number of particles
    * @return number of particles
    */
    int number() const;

    /**
    * access the position of the i-th particle, declared const, i.e. it is not possible to change the value
    * @param k: index of the particle
    * @return coordinate in the grid
    */
    double &operator()(int k, int axis);

    /**
    * get the position of the i-th particle, declared const, i.e. it is not possible to change the value
    * @param k: index of the particle
    * @return coordinate in the grid
    */
    double operator()(int k, int axis) const;

    void print() const;

    /**
     * get vector with particle information stored in the data vector
     */
    void *data();

    /**
     * set particle vector values to zero
     */
    void setToZero();

protected:
    std::vector<double> data_;      //< storage array values, in row-major order
    int number_; //< width, height of the domain
};
