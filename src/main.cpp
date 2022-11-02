#include "output_writer/output_writer_paraview.h"
#include "settings.h"
#include "discretization/2_central_differences.h"
#include "pressure_solver/sor.h"
#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace std;

int main(int argc, char *argv[]) {
    auto nCells = std::array<int, 2>{3, 3};
    auto meshWidth = std::array<double, 2>{1.0, 1.0};
    auto origin = std::array<double, 2>{meshWidth[0] / 2.0, meshWidth[1] / 2.0};
    auto d = new CentralDifferences(nCells, meshWidth);
    auto pRef = FieldVariable({nCells[0] + 2, nCells[1] + 2}, origin, meshWidth);
    auto rhs = d->rhs();
    for (int i = d->pIBegin() + 1; i < d->pIEnd(); i++) {
        for (int j = d->pJBegin() + 1; j < d->pJEnd(); j++) {
            pRef(i, j) = 10 * i + j;
        }
    }
    std::cout << "p = ..." << std::endl;
    pRef.print();

    for (int i = d->rhsIBegin(); i <= d->rhsIEnd(); i++) {
        for (int j = d->rhsJBegin(); j <= d->rhsJEnd(); j++) {
            rhs(i, j) = (pRef(i + 1, j) - 2 * pRef(i, j) + pRef(i - 1, j)) / pow(d->dx(),2) +
                        (pRef(i, j + 1) - 2 * pRef(i, j) + pRef(i, j - 1)) / pow(d->dy(),2);
        }
    }
    std::cout << "rhs = ..." << std::endl;
    rhs.print();

    double epsilon = 0.001;
    int maximumNumberOfIterations = 100;
    double omega = 0.1;
    auto sor = SOR(static_cast<std::shared_ptr<Discretization>>(d), epsilon, maximumNumberOfIterations, omega);
    sor.solve();
    return EXIT_SUCCESS;
    /*
    std::array < int, 2 > size = std::array < int, 2 > {3, 3};
    std::array < double, 2 > origin = std::array < double, 2 > {0.0, 0.0};
    std::array < double, 2 > meshWidth = std::array < double, 2 > {1.0, 1.0};
    FieldVariable fv = FieldVariable(size, origin, meshWidth);
    for (int i = 0; i < size[0]; i++) {
        for (int j = 0; j < size[1]; j++) {
            cout << "fv(" << i << "," << j << ") = " << fv(i, j) << endl;
        }
    }
    fv(1, 1) = 1.0;
    fv(2, 1) = 20;
    double x = 1.0;
    double y = 1.0;
    cout << "Enter x value: " << endl;
    cin >> x;
    cout << "Enter y value: " << endl;
    cin >> y;
    cout << "interpolation value: " << fv.interpolateAt(x, y) << endl;

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


    Settings settings;
    // load settings from file
    settings.loadFromFile(filename);

    // display all settings on console
    settings.printSettings();

    // write 5 output files
    for (int i = 0; i < 5; i++)
    {
        writeParaviewOutput(i);
    }*/

    return EXIT_SUCCESS;
}

