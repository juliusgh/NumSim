//#include "output_writer/output_writer_paraview.h"
#include "settings.h"
#include "computation/computation.h"
#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace std;

int main(int argc, char *argv[]) {
    // if the number of given command line arguments is only 1 (= the program name), print out usage information and exit
    string program = argv[0];
    string filename = "";
    if (argc == 1) {
        //cout << "usage: " << argv[0] << " <filename>" << endl;
        cout << "enter parameter file path: " << endl;
        cin >> filename;

        //return EXIT_FAILURE;
    }
    else {
        // read in the first argument
        filename = argv[1];
    }


    cout << "Program: \"" << program << "\"" << endl;

    // print message
    cout << "Filename: \"" << filename << "\"" << endl;

    auto computation = Computation();
    computation.initialize(filename);
    computation.runSimulation();

    // TODO: write actual output

    // write 5 output files
    /*for (int i = 0; i < 5; i++)
    {
        writeParaviewOutput(i);
    }*/

    return EXIT_SUCCESS;
}

