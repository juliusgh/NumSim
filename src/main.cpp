#include "output_writer/output_writer_paraview.h"
#include "settings.h"
#include "discretization/2_central_differences.h"
#include "pressure_solver/sor.h"
#include "pressure_solver/gauss_seidel.h"
#include "computation/computation.h"
#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace std;

int main(int argc, char *argv[]) {
    // if the number of given command line arguments is only 1 (= the program name), print out usage information and exit
    if (argc == 1) {
        cout << "usage: " << argv[0] << " <filename>" << endl;

        return EXIT_FAILURE;
    }

    // read in the first argument
    string program = argv[0];
    string filename = argv[1];

    cout << "Program: \"" << program << "\"" << endl;

    // print message
    cout << "Filename: \"" << filename << "\"" << endl;

    auto computation = Computation();
    computation.initialize(argc, argv);

    // write 5 output files
    /*for (int i = 0; i < 5; i++)
    {
        writeParaviewOutput(i);
    }*/

    return EXIT_SUCCESS;
}

