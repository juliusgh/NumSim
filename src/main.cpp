#include "output_writer/write_paraview_output.h"

#include <iostream>
#include <cstdlib>

int main(int argc, char *argv[])
{
  // write 5 output files
  for (int i = 0; i < 5; i++)
  {
    writeParaviewOutput(i);
  }

  return EXIT_SUCCESS;
}

