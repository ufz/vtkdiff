/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <cstdlib>

#include <tclap/CmdLine.h>

#include <vtkXMLUnstructuredGridReader.h>


int
main(int argc, char* argv[])
{
    // Parse CLI arguments.
    TCLAP::CmdLine cmd("VtkDiff software.\n"
            "Copyright (c) 2015, OpenGeoSys Community "
            "(http://www.opengeosys.org) "
            "Distributed under a Modified BSD License. "
            "See accompanying file LICENSE.txt or "
            "http://www.opengeosys.org/project/license",
        ' ',
        "0.1");

    TCLAP::UnlabeledValueArg<std::string> vtk_input_arg(
        "input-file",
        "Path to the VTK unstructured grid input file.",
        true,
        "",
        "VTK FILE");
    cmd.add(vtk_input_arg);

    TCLAP::ValueArg<std::string> data_array_a_arg(
        "a",
        "first",
        "First data array name for comparison",
        true,
        "",
        "NAME");
    cmd.add(data_array_a_arg);

    TCLAP::ValueArg<std::string> data_array_b_arg(
        "b",
        "second",
        "Second data array name for comparison",
        true,
        "",
        "NAME");
    cmd.add(data_array_b_arg);

    cmd.parse(argc, argv);

    return EXIT_SUCCESS;
}
