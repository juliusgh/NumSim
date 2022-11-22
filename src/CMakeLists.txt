cmake_minimum_required(VERSION 3.8)

# Define the project name.
project(numsim)

# Specify the name of the executable (${PROJECT_NAME} which is equal to what was set in the project() command).
# Also specify the source files.
add_executable(${PROJECT_NAME}
    main.cpp
    settings.cpp

    computation/0_computation.cpp
    computation/1_computation_parallel.cpp

    storage/array2d.cpp
    storage/fieldvariable.cpp

    discretization/0_staggered_grid.cpp
    discretization/1_discretization.cpp
    discretization/2_central_differences.cpp
    discretization/2_donor_cell.cpp

    pressure_solver/0_pressure_solver.cpp
    pressure_solver/1_gauss_seidel.cpp
    pressure_solver/1_sor.cpp

    output_writer/output_writer.cpp
    output_writer/output_writer_text.cpp
    output_writer/output_writer_paraview.cpp
)

# Locate GTest
find_package(GTest)

# If an installation of vtk was found on the system
if(GTest_FOUND)
    # Link runTests with what we want to test and the GTest and pthread library
    add_executable(runTests
        tests.cpp
        settings.cpp

        computation/0_computation.cpp
        computation/1_computation_parallel.cpp

        storage/array2d.cpp
        storage/fieldvariable.cpp

        discretization/0_staggered_grid.cpp
        discretization/1_discretization.cpp
        discretization/2_central_differences.cpp
        discretization/2_donor_cell.cpp

        pressure_solver/0_pressure_solver.cpp
        pressure_solver/1_gauss_seidel.cpp
        pressure_solver/1_sor.cpp

        output_writer/output_writer.cpp
        output_writer/output_writer_text.cpp
        output_writer/output_writer_paraview.cpp
    )

    target_link_libraries(runTests ${GTEST_LIBRARIES} pthread)
    target_include_directories(runTests PUBLIC ${PROJECT_SOURCE_DIR})
    include_directories(${GTEST_INCLUDE_DIRS})              # add the include directory where the header files are for the compiler
    target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES}) # add the libraries for the linker
    target_link_libraries(runTests ${VTK_LIBRARIES}) # add the libraries for the linker
endif(GTest_FOUND)

# Add the project directory to include directories, to be able to include all project header files from anywhere
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})

# Search for the external package "VTK"
find_package(VTK)

# If an installation of vtk was found on the system
if(VTK_FOUND)
    include_directories(${VTK_INCLUDE_DIRS})               # add the include directory where the header files are for the compiler
    target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES}) # add the libraries for the linker
    if(GTest_FOUND)
        target_link_libraries(runTests ${VTK_LIBRARIES}) # add the libraries for the linker
    endif(GTest_FOUND)
endif(VTK_FOUND)

# install numsim executable in build directory
install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/../build)

# Add additonial compile options to enable more warnings
add_compile_options(-Wall -Wextra)
add_compile_options(-O2)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

message("CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
