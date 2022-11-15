#include "storage/array2d.h"
#include <iostream>
#ifndef NDEBUG
#include <cassert>
#endif

/**
 * This class represents a 2D array of double values.
 *  Internally they are stored consecutively in memory.
 *  The entries can be accessed by two indices i,j.
 * @param size
 */

Array2D::Array2D(std::array<int, 2> size) :
        size_(size) {
    // allocate data, initialize to 0
    data_.resize(size_[0] * size_[1], 0.0);
}

//! get the size
std::array<int, 2> Array2D::size() const {
    return size_;
}

double &Array2D::operator()(int i, int j) {
    const int index = j * size_[0] + i;

    // assert that indices are in range
    #ifndef NDEBUG
    assert((0 <= i) && (i < size_[0]));
    assert((0 <= j) && (j < size_[1]));
    #endif
    assert(j * size_[0] + i < (int) data_.size());

    return data_[index];
}

double Array2D::operator()(int i, int j) const {
    const int index = j * size_[0] + i;

    // assert that indices are in range
    #ifndef NDEBUG
    assert((0 <= i) && (i < size_[0]));
    assert((0 <= j) && (j < size_[1]));
    #endif
    assert(j * size_[0] + i < (int) data_.size());

    return data_[index];
}

void Array2D::print() const {
    std::cout << std::endl << "----------" << std::endl;
    for (int j = size_[1] - 1; j >= 0; j--) {
        for (int i = 0; i < size_[0]; i++) {
            std::cout << (*this)(i, j) << " | ";
        }
        std::cout << std::endl << "----------" << std::endl;
    }
}