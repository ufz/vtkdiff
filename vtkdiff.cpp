/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <ios>

#include <tclap/CmdLine.h>

#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
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

    TCLAP::SwitchArg quite_arg(
        "q",
        "quite",
        "Suppress all but error output.");
    cmd.add(quite_arg);

    cmd.parse(argc, argv);

    bool const quite = quite_arg.getValue();

    // Setup the stdandard output and error stream numerical formats.
    std::cout << std::scientific << std::setprecision(16);
    std::cerr << std::scientific << std::setprecision(16);

    // Read input file.
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader
        = vtkXMLUnstructuredGridReader::New();

    reader->SetFileName(vtk_input_arg.getValue().c_str());
    reader->Update();

    vtkDataArray* a = reader->GetOutput()->GetPointData()->GetScalars(
        data_array_a_arg.getValue().c_str());
    if (a == nullptr)
    {
        std::cerr << "Error: Scalars data array "
            << "\'" << data_array_a_arg.getValue().c_str() << "\'"
            << "not found in point data.\n";
        return EXIT_FAILURE;
    }

    vtkDataArray* b = reader->GetOutput()->GetPointData()->GetScalars(
        data_array_b_arg.getValue().c_str());
    if (b == nullptr)
    {
        std::cerr << "Error: Scalars data array "
            << "\'" << data_array_b_arg.getValue().c_str() << "\'"
            << "not found in point data.\n";
        return EXIT_FAILURE;
    }

    // Check similarity of the data arrays.

    // Is numeric
    if (!a->IsNumeric())
    {
        std::cerr << "Data in data array a is not numeric:\n"
            << "data type is " << a->GetDataTypeAsString() << "\n";

        return EXIT_FAILURE;
    }
    if (!b->IsNumeric())
    {
        std::cerr << "Data in data array b is not numeric.\n"
            << "data type is " << b->GetDataTypeAsString() << "\n";
        return EXIT_FAILURE;
    }

    // Number of components
    if (a->GetNumberOfComponents() != b->GetNumberOfComponents())
    {
        std::cerr << "Number of components differ:\n"
            << a->GetNumberOfComponents() << " in data array a and "
            << b->GetNumberOfComponents() << " in data array b\n";
        return EXIT_FAILURE;
    }

    // For now only scalar type arrays are allowed.
    if (a->GetNumberOfComponents() != 1)
    {
        std::cerr << "Only scalar data arrays are supported.\n";
        return EXIT_FAILURE;
    }

    // Calculate difference of the data arrays.

    // Absolute error and norms.
    double abs_err_norm_l1 = 0;
    double abs_err_norm_2_2 = 0;
    double abs_err_norm_max = 0;
    vtkSmartPointer<vtkDoubleArray> abs_err = vtkDoubleArray::New();
    abs_err->SetNumberOfComponents(a->GetNumberOfComponents());
    abs_err->SetNumberOfTuples(a->GetNumberOfTuples());

    // Relative error and norms.
    double rel_err_norm_l1 = 0;
    double rel_err_norm_2_2 = 0;
    double rel_err_norm_max = 0;
    vtkSmartPointer<vtkDoubleArray> rel_err = vtkDoubleArray::New();
    rel_err->SetNumberOfComponents(a->GetNumberOfComponents());
    rel_err->SetNumberOfTuples(a->GetNumberOfTuples());

    for (vtkIdType i = 0; i < a->GetNumberOfTuples(); ++i)
    {
        // absolute error and its norms:
        abs_err->SetTuple1(i, a->GetTuple1(i) - b->GetTuple1(i));
        double const abs_err_i = abs_err->GetTuple1(i);

        abs_err_norm_l1 += abs_err_i;
        abs_err_norm_2_2 += abs_err_i*abs_err_i;
        abs_err_norm_max = std::max(abs_err_norm_max, abs_err_i);

        // relative error (to the data array a) and its norms:
        double const abs_a_i = std::abs(a->GetTuple1(i));
        if (abs_err_i == 0)
        {
            rel_err->SetTuple1(i, 0);
        }
        else if (abs_a_i == 0)
        {
            rel_err->SetTuple1(i, std::numeric_limits<double>::infinity());
        }
        else
        {
            rel_err->SetTuple1(i, abs_err_i / abs_a_i);
        }
        double const rel_err_i = rel_err->GetTuple1(i);

        rel_err_norm_l1 += rel_err_i;
        rel_err_norm_2_2 += rel_err_i*rel_err_i;
        rel_err_norm_max = std::max(rel_err_norm_max, rel_err_i);
    }

    // Error information
    if (!quite)
        std::cout << "Computed difference between data arrays:\n"
            << "abs l1 norm = " << abs_err_norm_l1 << "\n"
            << "abs 2-norm^2 = " << abs_err_norm_2_2 << "\n"
            << "abs 2-norm = " << std::sqrt(abs_err_norm_2_2) << "\n"
            << "abs maximum norm = " << abs_err_norm_max << "\n"
            << "\n"
            << "rel l1 norm = " << rel_err_norm_l1 << "\n"
            << "rel 2-norm^2 = " << rel_err_norm_2_2 << "\n"
            << "rel 2-norm = " << std::sqrt(rel_err_norm_2_2) << "\n"
            << "rel maximum norm = " << rel_err_norm_max << "\n";


    return EXIT_SUCCESS;
}
