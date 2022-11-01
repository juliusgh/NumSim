#include "output_writer/output_writer_paraview.h"
#include "settings.h"
#include <iostream>
#include <cstdlib>
using namespace std;

int main(int argc, char *argv[])
{
    // if the number of given command line arguments is only 1 (= the program name), print out usage information and exit
    if (argc == 1)
    {
        cout << "usage: " << argv[0] << " <filename>" << endl;

        return EXIT_FAILURE;
    }

    // read in the first argument
    string program = argv[0];
    string filename = argv[1];

    cout << "Program: \"" << program << "\"" << endl;

    // print message
    cout << "Filename: \"" << filename << "\"" << endl;


    Settings settings;
    // load settings from file
    settings.loadFromFile(filename);

    // display all settings on console
    settings.printSettings();

    // write 5 output files
    for (int i = 0; i < 5; i++)
    {
        writeParaviewOutput(i);
    }

    return EXIT_SUCCESS;
}

