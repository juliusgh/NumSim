#pragma once

#include <vector>
#include <array>

enum MARKER {
    FLUID, // 0
    FREE, // 1
    NOSLIP, // 2
    INFLOW, // 3
    OUTFLOW, // 4
    OBSTACLE, // 5
    OBSTACLE_LEFT, // 6
    OBSTACLE_RIGHT, // 7
    OBSTACLE_TOP, // 8
    OBSTACLE_BOTTOM, // 9
    OBSTACLE_LEFT_TOP, // 10
    OBSTACLE_RIGHT_TOP, // 11
    OBSTACLE_LEFT_BOTTOM, // 12
    OBSTACLE_RIGHT_BOTTOM, // 13
    SURFACE, // 14

};

/**
 * This class represents a 2D array of MARKER values.
 * Internally they are stored consecutively in memory.
 * The entries can be accessed by two indices i,j.
 */
class Marker2D {
public:
    /**
    * constructor
    * @param size: number of cells
    */
    Marker2D(std::array<int, 2> size);

    /**
    * get number of cells
    * @return number of cells
    */
    std::array<int, 2> size() const;

    /**
    * access the value at coordinate (i,j), declared not const, i.e. the value can be changed
    * @param i: discretized position in x direction
    * @param j: discretized position in y direction
    * @return reference to value at the grid cell (i,j)
    */
    MARKER &operator()(int i, int j);

    /**
    * get the value at coordinate (i,j), declared const, i.e. it is not possible to change the value
    * @param i: discretized position in x direction
    * @param j: discretized position in y direction
    * @return value at the grid cell (i,j)
    */
    MARKER operator()(int i, int j) const;

    void print() const;

    //! set all data of 2D array to FLUID cells
    void setToFluid();

    //! return data vector
    void *data();

protected:
    std::vector<MARKER> data_;      //< storage array values, in row-major order
    const std::array<int, 2> size_; //< width, height of the domain
};
