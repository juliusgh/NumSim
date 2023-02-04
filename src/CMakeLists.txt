cmake_minimum_required(VERSION 3.8)

# Define the project name.
project(numsim)

# Files for the project
set(PROJECT_FILES
        settings.cpp

        computation/0_computation.cpp
        computation/1_computation_parallel.cpp
        partitioning/partitioning.cpp

        storage/array2d.cpp
        storage/field_variable.cpp

        discretization/0_staggered_grid.cpp
        discretization/1_discretization.cpp
        discretization/2_central_differences.cpp
        discretization/2_donor_cell.cpp

        pressure_solver/0_pressure_solver.cpp
        pressure_solver/1_pressure_solver_parallel.cpp
        pressure_solver/1_gauss_seidel.cpp
        pressure_solver/1_sor.cpp
        pressure_solver/1_conjugate_gradient.cpp
        pressure_solver/2_red_black.cpp
        pressure_solver/2_red_black_sor.cpp
        pressure_solver/2_conjugate_gradient_parallel.cpp

        output_writer/output_writer.cpp
        output_writer/output_writer_text.cpp
        output_writer/output_writer_paraview.cpp
        output_writer/output_writer_text_parallel.cpp
        output_writer/output_writer_paraview_parallel.cpp
        storage/marker2d.cpp
        storage/particle2d.cpp
        )

# Executable for serial version
set(TARGET_SERIAL ${PROJECT_NAME})
add_executable(${TARGET_SERIAL} main.cpp ${PROJECT_FILES})
target_compile_definitions(${TARGET_SERIAL} PRIVATE -DNPARALLEL)
# Add the project directory to include directories, to be able to include all project header files from anywhere
target_include_directories(${TARGET_SERIAL} PUBLIC ${PROJECT_SOURCE_DIR})

# Executable for parallel version
find_package(MPI)
if (MPI_FOUND)
    set(TARGET_PARALLEL ${PROJECT_NAME}_parallel)
    add_executable(${TARGET_PARALLEL} main.cpp ${PROJECT_FILES})
    include_directories(${MPI_INCLUDE_PATH})
    target_link_libraries(${TARGET_PARALLEL} ${MPI_LIBRARIES})
    if (MPI_COMPILE_FLAGS)
        set_target_properties(${TARGET_PARALLEL} PROPERTIES
                COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
    endif (MPI_COMPILE_FLAGS)
    if (MPI_LINK_FLAGS)
        set_target_properties(${TARGET_PARALLEL} PROPERTIES
                LINK_FLAGS "${MPI_LINK_FLAGS}")
    endif (MPI_LINK_FLAGS)
    # Add the project directory to include directories, to be able to include all project header files from anywhere
    target_include_directories(${TARGET_PARALLEL} PUBLIC ${PROJECT_SOURCE_DIR})
    # install executable in build directory
    install(TARGETS ${TARGET_PARALLEL} RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/../build)
endif (MPI_FOUND)

# Locate GTest
find_package(GTest)

# If an installation of vtk was found on the system
if (GTest_FOUND)
    # Link runTests with what we want to test and the GTest and pthread library
    set(TARGET_TESTS ${PROJECT_NAME}_tests)
    add_executable(${TARGET_TESTS} tests.cpp ${PROJECT_FILES})

    target_link_libraries(${TARGET_TESTS} ${GTEST_LIBRARIES} pthread)
    target_include_directories(${TARGET_TESTS} PUBLIC ${PROJECT_SOURCE_DIR})
    include_directories(${GTEST_INCLUDE_DIRS})              # add the include directory where the header files are for the compiler
    target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES}) # add the libraries for the linker
    target_link_libraries(${TARGET_TESTS} ${VTK_LIBRARIES}) # add the libraries for the linker
    # install executable in build directory
    install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/../build)
endif (GTest_FOUND)

# Search for the external package "VTK"
find_package(VTK REQUIRED)

# If an installation of vtk was found on the system
if (VTK_FOUND)
    include_directories(${VTK_INCLUDE_DIRS})               # add the include directory where the header files are for the compiler
    target_link_libraries(${TARGET_SERIAL} ${VTK_LIBRARIES}) # add the libraries for the linker
    target_link_libraries(${TARGET_PARALLEL} ${VTK_LIBRARIES}) # add the libraries for the linker
    if (GTest_FOUND)
        target_link_libraries(${TARGET_TESTS} ${VTK_LIBRARIES}) # add the libraries for the linker
    endif (GTest_FOUND)
endif (VTK_FOUND)

# Add additonial compile options to enable more warnings
add_compile_options(-Wall -Wextra)
add_compile_options(-O2)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

message("CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

add_executable(compare_output comparison/main.cpp)

# VTK
find_package(VTK)
include_directories(${VTK_INCLUDE_DIRS})
target_link_libraries(compare_output ${VTK_LIBRARIES})

# Boost
find_package(Boost COMPONENTS filesystem system)
target_include_directories(compare_output PUBLIC ${BOOST_INCLUDE_DIRS})
target_link_libraries(compare_output ${Boost_LIBRARIES})


option(PROFILE "Add flags to profile the program with gprof." OFF)
if (PROFILE)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg -g")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg -g")
endif ()
