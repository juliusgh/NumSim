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
    //std::cout << "&Array2d(" << i << ", " << j << ")" << std::endl;
    const int index = j * size_[0] + i;

    // assert that indices are in range
    assert((0 <= i) && (i < size_[0]));
    assert((0 <= j) && (j < size_[1]));
    assert(j * size_[0] + i < (int) data_.size());

    return data_[index];
}

double Array2D::operator()(int i, int j) const {
    //std::cout << "Array2d(" << i << ", " << j << ")" << std::endl;
    const int index = j * size_[0] + i;

    // assert that indices are in range
    assert((0 <= i) && (i < size_[0]));
    assert((0 <= j) && (j < size_[1]));
    assert(j * size_[0] + i < (int) data_.size());

    return data_[index];
}

void Array2D::print() const {
    std::cout << std::endl << "----------" << std::endl;
    for (int i = size_[0] - 1; i >= 0; i--) {
        for (int j = 0; j < size_[1]; j++) {
            std::cout << (*this)(i, j) << " | ";
        }
        std::cout << std::endl << "----------" << std::endl;
    }
}