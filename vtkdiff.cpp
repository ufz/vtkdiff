/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <ios>
#include <sstream>
#include <tuple>
#include <type_traits>

#include <tclap/CmdLine.h>

#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGridReader.h>

template <typename T>
auto float_to_string(T const& v) -> std::string
{
    static_assert(std::is_floating_point<T>::value,
        "float_to_string requires a floating point input type.");

    std::stringstream double_eps_sstream;
    double_eps_sstream << std::scientific << std::setprecision(16) << v;
    return double_eps_sstream.str();
}

bool stringEndsWith(std::string const& str, std::string const& ending)
{
    if (str.length() < ending.length())
        return false;

    // now the difference is non-negative, no underflow possible.
    auto const string_end_length = str.length() - ending.length();
    return str.compare(string_end_length, ending.length(), ending) == 0;
}

struct Args
{
    bool const quite;
    double const abs_err_thr;
    double const rel_err_thr;
    std::string const vtk_input_a;
    std::string const vtk_input_b;
    std::string const data_array_a;
    std::string const data_array_b;
};

auto parseCommandLine(int argc, char* argv[]) -> Args
{
    TCLAP::CmdLine cmd("VtkDiff software.\n"
            "Copyright (c) 2015, OpenGeoSys Community "
            "(http://www.opengeosys.org) "
            "Distributed under a Modified BSD License. "
            "See accompanying file LICENSE.txt or "
            "http://www.opengeosys.org/project/license",
        ' ',
        "0.1");

    TCLAP::UnlabeledValueArg<std::string> vtk_input_a_arg(
        "input-file-a",
        "Path to the VTK unstructured grid input file.",
        true,
        "",
        "VTK FILE");
    cmd.add(vtk_input_a_arg);

    TCLAP::UnlabeledValueArg<std::string> vtk_input_b_arg(
        "input-file-b",
        "Path to the second VTK unstructured grid input file.",
        false,
        "",
        "VTK FILE");
    cmd.add(vtk_input_b_arg);

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

    auto const double_eps_string = float_to_string(
        std::numeric_limits<double>::epsilon());

    TCLAP::ValueArg<double> abs_err_thr_arg(
        "",
        "abs",
        "Threshold of the absolute error (" + double_eps_string + ")",
        false,
        std::numeric_limits<double>::epsilon(),
        "FLOAT");
    cmd.add(abs_err_thr_arg);

    TCLAP::ValueArg<double> rel_err_thr_arg(
        "",
        "rel",
        "Threshold of the relative error (" + double_eps_string + ")",
        false,
        std::numeric_limits<double>::epsilon(),
        "FLOAT");
    cmd.add(rel_err_thr_arg);

    cmd.parse(argc, argv);

    return Args
        { quite_arg.getValue()
        , abs_err_thr_arg.getValue()
        , abs_err_thr_arg.getValue()
        , vtk_input_a_arg.getValue()
        , vtk_input_b_arg.getValue()
        , data_array_a_arg.getValue()
        , data_array_b_arg.getValue()
        };
}

template<typename T>
auto readDataArraysFromFile(
        std::string const& file_a_name,
        std::string const& file_b_name,
        std::string const& data_array_a_name,
        std::string const& data_array_b_name)
    -> std::tuple<bool, vtkSmartPointer<vtkDataArray>, vtkSmartPointer<vtkDataArray>>
{
    // Read input file.
    auto reader_a = T::New();

    reader_a->SetFileName(file_a_name.c_str());
    reader_a->Update();

    // Get arrays
    auto a = vtkSmartPointer<vtkDataArray> {
        reader_a->GetOutput()->GetPointData()->GetScalars(
        data_array_a_name.c_str()) };
    vtkSmartPointer<vtkDataArray> b;
    if(file_b_name.size() == 0)
    {
        b = vtkSmartPointer<vtkDataArray> {
            reader_a->GetOutput()->GetPointData()->GetScalars(
            data_array_b_name.c_str()) };
    }
    else
    {
        auto reader_b = T::New();
        reader_b->SetFileName(file_b_name.c_str());
        reader_b->Update();
        b = vtkSmartPointer<vtkDataArray> {
            reader_b->GetOutput()->GetPointData()->GetScalars(
            data_array_b_name.c_str()) };
        reader_b->Delete();
    }

    reader_a->Delete();

    // Check arrays' validity
    if (!a)
    {
        std::cerr << "Error: Scalars data array "
            << "\'" << data_array_a_name.c_str() << "\'"
            << "not found in point data.\n";
        return std::make_tuple(false, nullptr, nullptr);
    }

    if (!b)
    {
        std::cerr << "Error: Scalars data array "
            << "\'" << data_array_b_name.c_str() << "\'"
            << "not found in point data.\n";
        return std::make_tuple(false, nullptr, nullptr);
    }

    return std::make_tuple(true, a, b);
}

int
main(int argc, char* argv[])
{
    auto const args = parseCommandLine(argc, argv);

    // Setup the stdandard output and error stream numerical formats.
    std::cout << std::scientific << std::setprecision(16);
    std::cerr << std::scientific << std::setprecision(16);

    // Read arrays from input file.
    bool read_successful;
    vtkSmartPointer<vtkDataArray> a;
    vtkSmartPointer<vtkDataArray> b;

    if (stringEndsWith(args.vtk_input_a, ".vtu"))
        std::tie(read_successful, a, b) =
            readDataArraysFromFile<vtkXMLUnstructuredGridReader>(
                args.vtk_input_a,
                args.vtk_input_b,
                args.data_array_a,
                args.data_array_b);
    else if (stringEndsWith(args.vtk_input_a, ".vtk"))
        std::tie(read_successful, a, b) =
            readDataArraysFromFile<vtkUnstructuredGridReader>(
                args.vtk_input_a,
                args.vtk_input_b,
                args.data_array_a,
                args.data_array_b);
    else
    {
        std::cerr << "Invalid file type! "
                     "Only .vtu and .vtk files are supported.\n";
        return EXIT_FAILURE;
    }

    if (!read_successful)
        return EXIT_FAILURE;

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
    auto abs_err = vtkSmartPointer<vtkDoubleArray>::New();
    abs_err->SetNumberOfComponents(a->GetNumberOfComponents());
    abs_err->SetNumberOfTuples(a->GetNumberOfTuples());

    // Relative error and norms.
    double rel_err_norm_l1 = 0;
    double rel_err_norm_2_2 = 0;
    double rel_err_norm_max = 0;
    auto rel_err = vtkSmartPointer<vtkDoubleArray>::New();
    rel_err->SetNumberOfComponents(a->GetNumberOfComponents());
    rel_err->SetNumberOfTuples(a->GetNumberOfTuples());

    for (auto i = 0; i < a->GetNumberOfTuples(); ++i)
    {
        // absolute error and its norms:
        abs_err->SetTuple1(i, a->GetTuple1(i) - b->GetTuple1(i));
        auto const abs_err_i = abs_err->GetTuple1(i);

        abs_err_norm_l1 += abs_err_i;
        abs_err_norm_2_2 += abs_err_i*abs_err_i;
        abs_err_norm_max = std::max(abs_err_norm_max, abs_err_i);

        // relative error (to the data array a) and its norms:
        auto const abs_a_i = std::abs(a->GetTuple1(i));
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
        auto const rel_err_i = rel_err->GetTuple1(i);

        rel_err_norm_l1 += rel_err_i;
        rel_err_norm_2_2 += rel_err_i*rel_err_i;
        rel_err_norm_max = std::max(rel_err_norm_max, rel_err_i);

        if (abs_err_i > args.abs_err_thr && rel_err_i > args.rel_err_thr && !args.quite)
        {
            std::cout << i << ": abs err = " << abs_err_i << "\n";
            std::cout << i << ": rel err = " << rel_err_i << "\n";
        }
    }

    // Error information
    if (!args.quite)
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

    if (abs_err_norm_max > args.abs_err_thr && rel_err_norm_max > args.rel_err_thr)
    {
        if (!args.quite)
            std::cout << "Absolute or relative error maximum is larger than "
                "the corresponding threshold.\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
