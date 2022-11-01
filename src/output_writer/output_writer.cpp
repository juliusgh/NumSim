#include "output_writer/output_writer.h"
#include <cassert>

//! constructor
OutputWriter::OutputWriter(std::shared_ptr<Discretization> discretization) :
discretization_(discretization)
{

};