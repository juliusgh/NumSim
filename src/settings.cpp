#include <fstream>   // for file operations
#include <iostream>  // for cout
#include <iomanip>
#include "settings.h"
using namespace std;

void Settings::loadFromFile(string filename) {
    // open file
    ifstream file(filename, ios::in);

    // check if file is open
    if (!file.is_open())
    {
        cout << "Could not open parameter file \"" << filename << "\"." << endl;
        return;
    }

    // loop over lines of file
    for (int lineNo = 0; ; lineNo++)
    {
        // read line
        string line;
        getline(file, line);

        // at the end of the file break for loop
        if (file.eof())
            break;

        // remove whitespace at beginning of line (if there is any)
        if (line.find_first_of(" \t") != string::npos)
        {
            line = line.substr(line.find_first_not_of(" \t"));
        }
        // if first character is a '#', skip line (line[0] == '#')
        if (line[0] == '#')
            continue;
        // if line does not contain a '=' sign, skip line
        if (line.find('=') == string::npos)
            continue;
        // parse parameter name
        string parameterName = line.substr(0, line.find('='));
        // remove trailing spaces from parameterName
        if (parameterName.find_first_of(" \t") != string::npos)
        {
            parameterName.erase(parameterName.find_first_of(" \t"));
        }
        // parse value
        string value = line.substr(line.find('=') + 1);
        // remove whitespace at beginning of value
        if (value.find_first_of(" \t") != string::npos)
        {
            value = value.substr(value.find_first_not_of(" \t"));
        }
        // remove comments at end of value
        if (value.find_first_of('#') != string::npos)
        {
            value = value.substr(0, value.find_first_of('#'));
        }
        // remove whitespace at end of value
        if (value.find_first_of(" \t") != string::npos)
        {
            value = value.substr(0, value.find_first_of(" \t"));
        }

        // parse actual value and set corresponding parameter
        if (parameterName == "nCellsX") {
            nCells[0] = atoi(value.c_str());
            continue;
        }
        if (parameterName == "nCellsY")
        {
            nCells[1] = atoi(value.c_str());
            continue;
        }
        if (parameterName == "physicalSizeX") {
            physicalSize[0] = atof(value.c_str());
            continue;
        }
        if (parameterName == "physicalSizeY") {
            physicalSize[1] = atof(value.c_str());
            continue;
        }
        if (parameterName == "re") {
            re = atof(value.c_str());
            continue;
        }
        if (parameterName == "beta") {
            beta = atof(value.c_str());
            continue;
        }
        if (parameterName == "endTime") {
            endTime = atof(value.c_str());
            continue;
        }
        if (parameterName == "tau") {
            tau = atof(value.c_str());
            continue;
        }
        if (parameterName == "maximumDt") {
            maximumDt = atof(value.c_str());
            continue;
        }
        if (parameterName == "gX") {
            g[0] = atof(value.c_str());
            continue;
        }
        if (parameterName == "gY") {
            g[1] = atof(value.c_str());
            continue;
        }
        if (parameterName == "useDonorCell") {
            istringstream(value.c_str()) >> boolalpha >> useDonorCell;
            continue;
        }
        if (parameterName == "alpha") {
            alpha = atof(value.c_str());
            continue;
        }
        if (parameterName == "gamma") {
            gamma = atof(value.c_str());
            continue;
        }
        if (parameterName == "dirichletBottomX") {
            dirichletBcBottom[0] = atof(value.c_str());
            continue;
        }
        if (parameterName == "dirichletBottomY") {
            dirichletBcBottom[1] = atof(value.c_str());
            continue;
        }
        if (parameterName == "dirichletTopX") {
            dirichletBcTop[0] = atof(value.c_str());
            continue;
        }
        if (parameterName == "dirichletTopY") {
            dirichletBcTop[1] = atof(value.c_str());
            continue;
        }
        if (parameterName == "dirichletLeftX") {
            dirichletBcLeft[0] = atof(value.c_str());
            continue;
        }
        if (parameterName == "dirichletLeftY") {
            dirichletBcLeft[1] = atof(value.c_str());
            continue;
        }
        if (parameterName == "dirichletRightX") {
            dirichletBcRight[0] = atof(value.c_str());
            continue;
        }
        if (parameterName == "dirichletRightY") {
            dirichletBcRight[1] = atof(value.c_str());
            continue;
        }
        if (parameterName == "pressureSolver") {
            pressureSolver = value;
            continue;
        }
        if (parameterName == "omega") {
            omega = atof(value.c_str());
            continue;
        }
        if (parameterName == "epsilon") {
            epsilon = atof(value.c_str());
            continue;
        }
        if (parameterName == "maximumNumberOfIterations") {
            maximumNumberOfIterations = static_cast<int>(atof(value.c_str()));
            continue;
        }
    }
}

void Settings::printSettings()
{
    cout << "Settings: " << endl
              << "  physicalSize: " << physicalSize[0] << " x " << physicalSize[1] << ", nCells: " << nCells[0] << " x " << nCells[1] << endl
              << "  endTime: " << endTime << " s, re: " << re << " beta: " << beta << ", g: (" << g[0] << "," << g[1] << "), tau: " << tau << ", maximum dt: " << maximumDt << endl
              << "  dirichletBC: bottom: (" << dirichletBcBottom[0] << "," << dirichletBcBottom[1]  << ")"
              << ", top: ("  << dirichletBcTop[0] << "," << dirichletBcTop[1]  << ")"
              << ", left: ("  << dirichletBcLeft[0] << "," << dirichletBcLeft[1] << ")"
              << ", right: ("  << dirichletBcRight[0] << "," << dirichletBcRight[1] << ")" << endl
              << "  useDonorCell: " << boolalpha << useDonorCell << ", alpha: " << alpha << ", gamma: " << gamma << endl
              << "  pressureSolver: " << pressureSolver << ", omega: " << omega << ", epsilon: " << epsilon << ", maximumNumberOfIterations: " << maximumNumberOfIterations << endl;
}
