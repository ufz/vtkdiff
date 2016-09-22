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

#include <vtkCommand.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>

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
    bool const quiet;
    bool const verbose;
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

    TCLAP::SwitchArg quiet_arg(
        "q",
        "quiet",
        "Suppress all but error output.");
    cmd.add(quiet_arg);

    TCLAP::SwitchArg verbose_arg(
        "v",
        "verbose",
        "Also print which values differ.");
    cmd.add(verbose_arg);

    auto const double_eps_string = float_to_string(
        std::numeric_limits<double>::epsilon());

    TCLAP::ValueArg<double> abs_err_thr_arg(
        "",
        "abs",
        "Tolerance for the absolute error in the maximum norm (" + double_eps_string + ")",
        false,
        std::numeric_limits<double>::epsilon(),
        "FLOAT");
    cmd.add(abs_err_thr_arg);

    TCLAP::ValueArg<double> rel_err_thr_arg(
        "",
        "rel",
        "Tolerance for the componentwise relative error (" + double_eps_string + ")",
        false,
        std::numeric_limits<double>::epsilon(),
        "FLOAT");
    cmd.add(rel_err_thr_arg);

    cmd.parse(argc, argv);

    return Args
        { quiet_arg.getValue()
        , verbose_arg.getValue()
        , abs_err_thr_arg.getValue()
        , rel_err_thr_arg.getValue()
        , vtk_input_a_arg.getValue()
        , vtk_input_b_arg.getValue()
        , data_array_a_arg.getValue()
        , data_array_b_arg.getValue()
        };
}

template <typename T>
class ErrorCallback : public vtkCommand
{
public:
    vtkTypeMacro(ErrorCallback, vtkCommand);

    static ErrorCallback<T>* New() { return new ErrorCallback<T>; }

    void Execute(vtkObject* caller, unsigned long vtkNotUsed(eventId),
                 void* vtkNotUsed(callData))
    {
        auto* reader = static_cast<T*>(caller);
        std::cerr << "Error reading file `" << reader->GetFileName()
                  << "'. Aborting." << std::endl;
        std::abort();
    }
};

template<typename T>
std::tuple<bool, vtkSmartPointer<vtkDataArray>, vtkSmartPointer<vtkDataArray>>
readDataArraysFromFile(
        std::string const& file_a_name,
        std::string const& file_b_name,
        std::string const& data_array_a_name,
        std::string const& data_array_b_name)
{
    vtkSmartPointer<ErrorCallback<T>> errorCallback =
        vtkSmartPointer<ErrorCallback<T>>::New();

    // Read input file.
    vtkSmartPointer<T> reader_a = vtkSmartPointer<T>::New();
    reader_a->AddObserver(vtkCommand::ErrorEvent, errorCallback);
    reader_a->SetFileName(file_a_name.c_str());
    reader_a->Update();

    bool point_data(false);
    if (reader_a->GetOutput()->GetPointData()->HasArray(
            data_array_a_name.c_str()))
        point_data = true;
    else if (reader_a->GetOutput()->GetCellData()->HasArray(
                 data_array_a_name.c_str()))
        point_data = false;
    else
    {
        std::cerr << "Error: Scalars data array "
            << "\'" << data_array_a_name.c_str() << "\'"
            << " neither found in point data nor in cell data.\n";
        return std::make_tuple(false, nullptr, nullptr);
    }

    // Get arrays
    vtkSmartPointer<vtkDataArray> a;
    if (point_data)
        a = vtkSmartPointer<vtkDataArray>{
            reader_a->GetOutput()->GetPointData()->GetScalars(
                data_array_a_name.c_str())};
    else
        a = vtkSmartPointer<vtkDataArray>{
            reader_a->GetOutput()->GetCellData()->GetScalars(
                data_array_a_name.c_str())};

    // Check arrays' validity
    if (!a)
    {
        std::cerr << "Error: Scalars data array "
            << "\'" << data_array_a_name.c_str() << "\'"
            << " could not be read.\n";
        return std::make_tuple(false, nullptr, nullptr);
    }

    vtkSmartPointer<vtkDataArray> b;
    if (file_b_name.empty()) {
        if (data_array_a_name == data_array_b_name) {
            std::cerr << "Error: You are trying to compare data array `"
                      << data_array_a_name << "' from file `" << file_a_name
                      << "' to itself. Aborting.\n";
            std::abort();
        }
        if (point_data)
            b = vtkSmartPointer<vtkDataArray>{
                reader_a->GetOutput()->GetPointData()->GetScalars(
                    data_array_b_name.c_str())};
        else
            b = vtkSmartPointer<vtkDataArray>{
                reader_a->GetOutput()->GetCellData()->GetScalars(
                    data_array_b_name.c_str())};
    } else {
        vtkSmartPointer<T> reader_b = vtkSmartPointer<T>::New();
        reader_b->AddObserver(vtkCommand::ErrorEvent, errorCallback);
        reader_b->SetFileName(file_b_name.c_str());
        reader_b->Update();
        if (point_data)
            b = vtkSmartPointer<vtkDataArray>{
                reader_b->GetOutput()->GetPointData()->GetScalars(
                    data_array_b_name.c_str())};
        else
            b = vtkSmartPointer<vtkDataArray>{
                reader_b->GetOutput()->GetCellData()->GetScalars(
                    data_array_b_name.c_str())};
    }

    if (!b)
    {
        std::cerr << "Error: Scalars data array "
            << "\'" << data_array_b_name.c_str() << "\'"
            << " not found.\n";
        return std::make_tuple(false, nullptr, nullptr);
    }

    return std::make_tuple(true, a, b);
}

int main(int argc, char* argv[])
{
    auto const digits10 = std::numeric_limits<double>::digits10;
    auto const args = parseCommandLine(argc, argv);

    // Setup the standard output and error stream numerical formats.
    std::cout << std::scientific << std::setprecision(digits10);
    std::cerr << std::scientific << std::setprecision(digits10);

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
    else
    {
        std::cerr << "Invalid file type! "
                     "Only .vtu files are supported.\n";
        return EXIT_FAILURE;
    }

    if (!read_successful)
        return EXIT_FAILURE;

    if (!args.quiet)
        std::cout << "Comparing data array `" << args.data_array_a
            << "' from file `" << args.vtk_input_a
            << "' to data array `" << args.data_array_b
            << "' from file `" << args.vtk_input_b << "'.\n";

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

    auto const num_tuples = a->GetNumberOfTuples();
    // Number of components
    if (num_tuples != b->GetNumberOfTuples())
    {
        std::cerr << "Number of tuples differ:\n"
            << num_tuples << " in data array a and "
            << b->GetNumberOfTuples() << " in data array b\n";
        return EXIT_FAILURE;
    }

    auto const num_components = a->GetNumberOfComponents();
    // Number of components
    if (num_components != b->GetNumberOfComponents())
    {
        std::cerr << "Number of components differ:\n"
            << num_components << " in data array a and "
            << b->GetNumberOfComponents() << " in data array b\n";
        return EXIT_FAILURE;
    }

    // Calculate difference of the data arrays.

    // Absolute error and norms.
    double abs_err_norm_l1 = 0;
    double abs_err_norm_2_2 = 0;
    double abs_err_norm_max = 0;

    // Relative error and norms.
    double rel_err_norm_l1 = 0;
    double rel_err_norm_2_2 = 0;
    double rel_err_norm_max = 0;

    for (auto tuple_idx = 0; tuple_idx < num_tuples; ++tuple_idx)
    {
        for (auto component_idx = 0; component_idx < num_components; ++component_idx)
        {
            auto const a_comp = a->GetComponent(tuple_idx, component_idx);
            auto const b_comp = b->GetComponent(tuple_idx, component_idx);
            auto const abs_err = std::abs(a_comp - b_comp);

            abs_err_norm_l1  += abs_err;
            abs_err_norm_2_2 += abs_err*abs_err;
            abs_err_norm_max  = std::max(abs_err_norm_max, abs_err);

            // relative error and its norms:
            double rel_err;

            if (abs_err == 0.0) {
                rel_err = 0.0;
            } else if (a_comp == 0.0 || b_comp == 0.0) {
                rel_err = std::numeric_limits<double>::infinity();
            } else {
                rel_err = abs_err / std::min(std::abs(a_comp), std::abs(b_comp));
            }

            rel_err_norm_l1  += rel_err;
            rel_err_norm_2_2 += rel_err*rel_err;
            rel_err_norm_max  = std::max(rel_err_norm_max, rel_err);

            if (abs_err > args.abs_err_thr && rel_err > args.rel_err_thr && args.verbose)
            {
                std::cout << std::setw(4) << tuple_idx
                          << ": abs err = " << std::setw(digits10+7) << abs_err
                          << ", rel err = " << std::setw(digits10+7) << rel_err << "\n";
            }
        }
    }

    // Error information
    if (!args.quiet)
        std::cout << "Computed difference between data arrays:\n"
            << "abs l1 norm      = " << abs_err_norm_l1 << "\n"
            << "abs l2-norm^2    = " << abs_err_norm_2_2 << "\n"
            << "abs l2-norm      = " << std::sqrt(abs_err_norm_2_2) << "\n"
            << "abs maximum norm = " << abs_err_norm_max << "\n"
            << "\n"
            << "rel l1 norm      = " << rel_err_norm_l1 << "\n"
            << "rel l2-norm^2    = " << rel_err_norm_2_2 << "\n"
            << "rel l2-norm      = " << std::sqrt(rel_err_norm_2_2) << "\n"
            << "rel maximum norm = " << rel_err_norm_max << "\n";

    if (abs_err_norm_max > args.abs_err_thr && rel_err_norm_max > args.rel_err_thr)
    {
        if (!args.quiet)
            std::cout << "Absolute and relative error (maximum norm) are larger"
                         " than the corresponding thresholds.\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
