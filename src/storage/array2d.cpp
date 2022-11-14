#include "storage/array2d.h"
#include <iostream>
#include <cassert>

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
    assert((0 <= i) && (i < size_[0]));
    assert((0 <= j) && (j < size_[1]));
    assert(j * size_[0] + i < (int) data_.size());

    return data_[index];
}

double Array2D::operator()(int i, int j) const {
    const int index = j * size_[0] + i;

    // assert that indices are in range
    assert((0 <= i) && (i < size_[0]));
    assert((0 <= j) && (j < size_[1]));
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