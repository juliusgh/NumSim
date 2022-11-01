cmake_minimum_required(VERSION 3.8)

# Define the project name.
project(numsim)

# Specify the name of the executable (${PROJECT_NAME} which is equal to what was set in the project() command).
# Also specify the source files.
add_executable(${PROJECT_NAME}
        main.cpp
        output_writer/output_writer.cpp output_writer/output_writer.h
        output_writer/output_writer_paraview.cpp output_writer/output_writer_paraview.h
        output_writer/output_writer_text.cpp output_writer/output_writer_text.h
        settings.cpp settings.h
        storage/array2d.cpp storage/array2d.h
        storage/fieldvariable.cpp storage/fieldvariable.h
        discretization/0_staggered_grid.cpp discretization/0_staggered_grid.h
        discretization/1_discretization.cpp discretization/1_discretization.h
        discretization/2_central_differences.cpp discretization/2_central_differences.h
        discretization/2_donor_cell.cpp discretization/2_donor_cell.h
        pressure_solver/pressure_solver.cpp pressure_solver/pressure_solver.h
        pressure_solver/gauss_seidel.cpp pressure_solver/gauss_seidel.h
        pressure_solver/sor.cpp pressure_solver/sor.h
)

# Add the project directory to include directories, to be able to include all project header files from anywhere
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})

# Locate GTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# Link runTests with what we want to test and the GTest and pthread library
add_executable(runTests
        tests.cpp
        output_writer/output_writer.cpp output_writer/output_writer.h
        output_writer/output_writer_text.cpp output_writer/output_writer_text.h
        settings.cpp settings.h
        storage/array2d.cpp storage/array2d.h
        storage/fieldvariable.cpp storage/fieldvariable.h
        discretization/0_staggered_grid.cpp discretization/0_staggered_grid.h
        discretization/1_discretization.cpp discretization/1_discretization.h
        discretization/2_central_differences.cpp discretization/2_central_differences.h
        discretization/2_donor_cell.cpp discretization/2_donor_cell.h
        pressure_solver/pressure_solver.cpp pressure_solver/pressure_solver.h
        pressure_solver/gauss_seidel.cpp pressure_solver/gauss_seidel.h
        pressure_solver/sor.cpp pressure_solver/sor.h
        )

target_link_libraries(runTests ${GTEST_LIBRARIES} pthread)
target_include_directories(runTests PUBLIC ${PROJECT_SOURCE_DIR})

# Search for the external package "VTK"
find_package(VTK)

# Output various cmake variables for demonstration purpose
message("If VTK was found on the system: VTK_FOUND: ${VTK_FOUND}")
message("The directory of VTK: VTK_DIR: ${VTK_DIR}")
message("The include directory of VTK: VTK_INCLUDE_DIRS: ${VTK_INCLUDE_DIRS}")

# If an installation of vtk was found on the system
if (VTK_FOUND)
  include_directories(${VTK_INCLUDE_DIRS})               # add the include directory where the header files are for the compiler
  target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES}) # add the libraries for the linker
endif(VTK_FOUND)

# install numsim executable in build directory
install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/../build)

# Add additonial compile options to enable more warnings
add_compile_options(-Wall -Wextra)

# Set the version of the C++ standard to use, we use C++14, published in 2014
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

message("CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
