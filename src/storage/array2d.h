#pragma once

#include <vector>
#include <array>

/** This class represents a 2D array of double values.
 *  Internally they are stored consecutively in memory.
 *  The entries can be accessed by two indices i,j.
 */
class Array2D
{
public:
  //! constructor
  Array2D(std::array<int,2> size);

  //! get the size
  std::array<int,2> size() const;

  //! access the value at coordinate (i,j), declared not const, i.e. the value can be changed
  double &operator()(int i, int j);

  //! get the value at coordinate (i,j), declared const, i.e. it is not possible to change the value
  double operator()(int i, int j) const;

protected:

  std::vector<double> data_;  //< storage array values, in row-major order
  const std::array<int,2> size_;    //< width, height of the domain
};
