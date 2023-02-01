#include "storage/array2d.h"
#include <iostream>
#include <iomanip>
#include <cassert>

/**
 * This class represents a 2D array of double values.
 *  Internally they are stored consecutively in memory.
 *  The entries can be accessed by two indices i,j.
 * @param size: number of cells
 */

Array2D::Array2D(std::array<int, 2> size) :
        size_(size) {
    // allocate data, initialize to 0
    data_.resize(size_[0] * size_[1], 0.0);
}

/**
 * get number of cells
 * @return number of cells
 */
std::array<int, 2> Array2D::size() const {
    return size_;
}

/**
 * access the value at coordinate (i,j), declared not const, i.e. the value can be changed
 * @param i: discretized position in x direction
 * @param j: discretized position in y direction
 * @return reference to value at the grid cell (i,j)
 */
double &Array2D::operator()(int i, int j) {
    const int index = j * size_[0] + i;

    // assert that indices are in range
#ifndef NDEBUG
    assert((0 <= i) && (i < size_[0]));
    assert((0 <= j) && (j < size_[1]));
    assert(j * size_[0] + i < (int) data_.size());
#endif

    return data_[index];
}

/**
 * get the value at coordinate (i,j), declared const, i.e. it is not possible to change the value
 * @param i: discretized position in x direction
 * @param j: discretized position in y direction
 * @return value at the grid cell (i,j)
 */
double Array2D::operator()(int i, int j) const {
    const int index = j * size_[0] + i;

    // assert that indices are in range
#ifndef NDEBUG
    assert((0 <= i) && (i < size_[0]));
    assert((0 <= j) && (j < size_[1]));
    assert(j * size_[0] + i < (int) data_.size());
#endif

    return data_[index];
}

/**
 * print out grid in two dimensions
 */

void Array2D::print() const {
    std::cout << std::endl << "----------" << std::endl;
    for (int j = size_[1] - 1; j >= 0; j--) {
        for (int i = 0; i < size_[0]; i++) {
            std::cout << std::setw(6) << std::setprecision(3) <<  (*this)(i, j) << " | ";
        }
        std::cout << std::endl;
    }
    std::cout << "----------" << std::endl;
}

void Array2D::setToZero() {
    data_.resize(size_[0] * size_[1], 0.0);
}

void *Array2D::data() {
    return data_.data();
}