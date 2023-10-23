#include "storage/particle2d.h"
#include <iostream>
#include <cassert>

/**
 * This class represents a 2D array of double values.
 *  Internally they are stored consecutively in memory.
 *  The entries can be accessed by two indices i,j.
 * @param number: number of particles
 */
Particle2D::Particle2D(int number) :
        number_(number) {
    // allocate data, initialize to FLUID
    setToZero();
}

/**
 * get number of particles
 * @return number of cells
 */
int Particle2D::number() const {
    return number_;
}

/**
* access the position of the k-th particle as reference
* @param k: index of the particle
* @param axis: axis (0 for x, 1 for y)
* @return coordinate in the grid along given axis (0 for x, 1 for y)
*/
double &Particle2D::operator()(int k, int axis) {
#ifndef NDEBUG
    assert((0 <= k) && (k < number_));
    assert((0 <= axis) && (axis < 2));
#endif
    const int index = k * 2 + axis;

    // assert that indices are in range
#ifndef NDEBUG
    assert((0 <= index) && (index < number_ * 2));
#endif

    return data_[index];
}

/**
* get the position of the k-th particle, declared const, i.e. it is not possible to change the value
* @param k: index of the particle
* @param axis: axis (0 for x, 1 for y)
* @return coordinate in the grid along given axis (0 for x, 1 for y)
*/
double Particle2D::operator()(int k, int axis) const {
#ifndef NDEBUG
    assert((0 <= k) && (k < number_));
    assert((0 <= axis) && (axis < 2));
#endif
    const int index = k * 2 + axis;

    // assert that indices are in range
#ifndef NDEBUG
    assert((0 <= index) && (index < number_ * 2));
#endif

    return data_[index];
}

/**
 * print out grid in two dimensions
 */
void Particle2D::print() const {
    std::cout << std::endl;
    for (int k = 0; k < number_; k++) {
        std::cout << (*this)(k, 0) << "|"<< (*this)(k,1) << std::endl;
    }
}

/**
 * get vector with particle information stored in the data vector
 */
void Particle2D::setToZero() {
    data_.resize(number_ * 2, 0.0);
}

/**
 * set particle vector values to zero
 */
void *Particle2D::data() {
    return data_.data();
}
