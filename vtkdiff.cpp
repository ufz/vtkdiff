/**
 * \copyright
 * Copyright (c) 2015-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <ios>
#include <iterator>
#include <sstream>
#include <tuple>
#include <type_traits>

#include <tclap/CmdLine.h>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCommand.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVersion.h>
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

template <typename T>
std::ostream& operator<<(std::ostream& os, std::vector<T> const& vector)
{
    if (vector.empty())
    {
        return os << "[]";
    }

    // print first n-1 elements
    os << "[";
    std::size_t const size = vector.size();
    for (std::size_t i = 0; i < size - 1; ++i)
    {
        os << vector[i] << ", ";
    }
    return os << vector.back() << "]";
}

struct Args
{
    bool const quiet;
    bool const verbose;
    bool const meshcheck;
    double const abs_err_thr;
    double const rel_err_thr;
    std::string const vtk_input_a;
    std::string const vtk_input_b;
    std::string const data_array_a;
    std::string const data_array_b;
};

auto parseCommandLine(int argc, char* argv[]) -> Args
{
    TCLAP::CmdLine cmd(
        "VtkDiff software.\n"
        "Copyright (c) 2015-2020, OpenGeoSys Community "
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
        "a", "first_data_array", "First data array name for comparison", true,
        "", "NAME");

    TCLAP::ValueArg<std::string> data_array_b_arg(
        "b",
        "second_data_array",
        "Second data array name for comparison",
        false,
        "",
        "NAME");
    cmd.add(data_array_b_arg);

    TCLAP::SwitchArg meshcheck_arg(
        "m", "mesh_check", "Compare mesh geometries using absolute tolerance.");
    cmd.xorAdd(data_array_a_arg, meshcheck_arg);

    TCLAP::SwitchArg quiet_arg("q", "quiet", "Suppress all but error output.");
    cmd.add(quiet_arg);

    TCLAP::SwitchArg verbose_arg("v", "verbose",
                                 "Also print which values differ.");
    cmd.add(verbose_arg);

    auto const double_eps_string =
        float_to_string(std::numeric_limits<double>::epsilon());

    TCLAP::ValueArg<double> abs_err_thr_arg(
        "",
        "abs",
        "Tolerance for the absolute error in the maximum norm (" +
            double_eps_string + ")",
        false,
        std::numeric_limits<double>::epsilon(),
        "FLOAT");
    cmd.add(abs_err_thr_arg);

    TCLAP::ValueArg<double> rel_err_thr_arg(
        "",
        "rel",
        "Tolerance for the componentwise relative error (" + double_eps_string +
            ")",
        false,
        std::numeric_limits<double>::epsilon(),
        "FLOAT");
    cmd.add(rel_err_thr_arg);

    cmd.parse(argc, argv);

    return Args{quiet_arg.getValue(),       verbose_arg.getValue(),
                meshcheck_arg.getValue(),   abs_err_thr_arg.getValue(),
                rel_err_thr_arg.getValue(), vtk_input_a_arg.getValue(),
                vtk_input_b_arg.getValue(), data_array_a_arg.getValue(),
                data_array_b_arg.getValue()};
}

template <typename T>
class ErrorCallback : public vtkCommand
{
public:
    vtkTypeMacro(ErrorCallback, vtkCommand);

    static ErrorCallback<T>* New() { return new ErrorCallback<T>; }

    void Execute(vtkObject* caller, unsigned long vtkNotUsed(eventId),
                 void* callData) override
    {
        auto* reader = static_cast<T*>(caller);
        std::cerr << "Error reading file `" << reader->GetFileName() << "'\n"
                  << static_cast<char*>(callData) << "\nAborting." << std::endl;
        std::exit(2);
    }
};

vtkSmartPointer<vtkUnstructuredGrid> readMesh(std::string const& filename)
{
    if (filename.empty())
    {
        return nullptr;
    }

    if (!stringEndsWith(filename, ".vtu"))
    {
        std::cerr << "Error: Expected a file with .vtu extension."
                  << "File '" << filename << "' not read.";
        return nullptr;
    }

    vtkSmartPointer<ErrorCallback<vtkXMLUnstructuredGridReader>> errorCallback =
        vtkSmartPointer<ErrorCallback<vtkXMLUnstructuredGridReader>>::New();

    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
        vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->AddObserver(vtkCommand::ErrorEvent, errorCallback);
    reader->SetFileName(filename.c_str());
    reader->Update();
    return reader->GetOutput();
}

std::tuple<vtkSmartPointer<vtkUnstructuredGrid>,
           vtkSmartPointer<vtkUnstructuredGrid>>
readMeshes(std::string const& file_a_name, std::string const& file_b_name)
{
    return {readMesh(file_a_name), readMesh(file_b_name)};
}

std::tuple<bool, vtkSmartPointer<vtkDataArray>, vtkSmartPointer<vtkDataArray>>
readDataArraysFromMeshes(
    std::tuple<vtkSmartPointer<vtkUnstructuredGrid>,
               vtkSmartPointer<vtkUnstructuredGrid>> const& meshes,
    std::string const& data_array_a_name,
    std::string const& data_array_b_name)
{
    if (std::get<0>(meshes) == nullptr)
    {
        std::cerr << "First mesh was not read correctly and is a nullptr.\n";
        return {false, nullptr, nullptr};
    }

    bool point_data(false);
    if (std::get<0>(meshes)->GetPointData()->HasArray(
            data_array_a_name.c_str()))
    {
        point_data = true;
    }
    else if (std::get<0>(meshes)->GetCellData()->HasArray(
                 data_array_a_name.c_str()))
    {
        point_data = false;
    }
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
    {
        a = vtkSmartPointer<vtkDataArray>{
            std::get<0>(meshes)->GetPointData()->GetScalars(
                data_array_a_name.c_str())};
    }
    else
    {
        a = vtkSmartPointer<vtkDataArray>{
            std::get<0>(meshes)->GetCellData()->GetScalars(
                data_array_a_name.c_str())};
    }

    // Check arrays' validity
    if (!a)
    {
        std::cerr << "Error: Scalars data array "
                  << "\'" << data_array_a_name.c_str() << "\'"
                  << " could not be read.\n";
        return std::make_tuple(false, nullptr, nullptr);
    }

    vtkSmartPointer<vtkDataArray> b;
    if (std::get<1>(meshes) == nullptr)
    {
        if (data_array_a_name == data_array_b_name)
        {
            std::cerr << "Error: You are trying to compare data array `"
                      << data_array_a_name
                      << "' from first file to itself. Aborting.\n";
            std::exit(3);
        }
        if (point_data)
        {
            b = vtkSmartPointer<vtkDataArray>{
                std::get<0>(meshes)->GetPointData()->GetScalars(
                    data_array_b_name.c_str())};
        }
        else
        {
            b = vtkSmartPointer<vtkDataArray>{
                std::get<0>(meshes)->GetCellData()->GetScalars(
                    data_array_b_name.c_str())};
        }
    }
    else
    {
        if (point_data)
        {
            b = vtkSmartPointer<vtkDataArray>{
                std::get<1>(meshes)->GetPointData()->GetScalars(
                    data_array_b_name.c_str())};
        }
        else
        {
            b = vtkSmartPointer<vtkDataArray>{
                std::get<1>(meshes)->GetCellData()->GetScalars(
                    data_array_b_name.c_str())};
        }
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

bool compareCellTopology(vtkCellArray* const cells_a,
                         vtkCellArray* const cells_b)
{
    vtkIdType const n_cells_a{cells_a->GetNumberOfCells()};
    vtkIdType const n_cells_b{cells_b->GetNumberOfCells()};

    if (n_cells_a != n_cells_b)
    {
        std::cerr << "Number of cells in the first mesh is " << n_cells_a
                  << " and differs from the number of cells in the second "
                     "mesh, which is "
                  << n_cells_b << "\n";
        return false;
    }

    vtkIdType n_cell_points_a, n_cell_points_b;
    #if (VTK_MAJOR_VERSION > 8 || VTK_MINOR_VERSION == 90)
        const vtkIdType *cell_points_a, *cell_points_b;
    #else
        vtkIdType *cell_points_a, *cell_points_b;
    #endif
    cells_a->InitTraversal();
    cells_b->InitTraversal();
    int get_next_cell_a = cells_a->GetNextCell(n_cell_points_a, cell_points_a);
    int get_next_cell_b = cells_b->GetNextCell(n_cell_points_b, cell_points_b);
    int cell_number = 0;
    while (get_next_cell_a == 1 && get_next_cell_b == 1)
    {
        if (n_cell_points_a != n_cell_points_b)
        {
            std::cerr << "Cell " << cell_number << " in first input has "
                      << n_cell_points_a << " points but in the second input "
                      << n_cell_points_b << " points.\n";
        }

        for (vtkIdType i = 0; i < n_cell_points_a; ++i)
        {
            if (cell_points_a[i] != cell_points_b[i])
            {
                std::cerr << "Point " << i << " of cell " << cell_number
                          << " has id " << cell_points_a[i]
                          << " in the first input but id " << cell_points_b[i]
                          << " in the second input.\n";
                return false;
            }
        }

        get_next_cell_a = cells_a->GetNextCell(n_cell_points_a, cell_points_a);
        get_next_cell_b = cells_b->GetNextCell(n_cell_points_b, cell_points_b);
        cell_number++;
    }

    if (get_next_cell_a != 0)
    {
        std::cerr << "Unexpected return value (" << get_next_cell_a
                  << ") for cells_a->GetNextCell() call. Expected 0 for "
                     "end-of-list or 1 for no error.\n";
        return false;
    }
    if (get_next_cell_b != 0)
    {
        std::cerr << "Unexpected return value (" << get_next_cell_b
                  << ") for cells_b->GetNextCell() call. Expected 0 for "
                     "end-of-list or 1 for no error.\n";
        return false;
    }

    return true;
}

bool comparePoints(vtkPoints* const points_a, vtkPoints* const points_b,
                   double const eps_squared)
{
    vtkIdType const n_points_a{points_a->GetNumberOfPoints()};
    vtkIdType const n_points_b{points_b->GetNumberOfPoints()};

    if (n_points_a != n_points_b)
    {
        std::cerr << "Number of points in the first mesh is " << n_points_a
                  << " and differst from the number of point in the second "
                     "mesh, which is "
                  << n_points_b << "\n";
        return false;
    }

    for (vtkIdType p = 0; p < n_points_a; ++p)
    {
        auto const a = points_a->GetPoint(p);
        auto const b = points_b->GetPoint(p);
        double const distance2 = vtkMath::Distance2BetweenPoints(a, b);
        if (distance2 >= eps_squared)
        {
            std::cerr << "Point " << p << " with coordinates (" << a[0] << ", "
                      << a[1] << ", " << a[2]
                      << ") from the first mesh is significantly different "
                         "from the same point in the second mesh, which "
                         "has coordinates ("
                      << b[0] << ", " << b[1] << ", " << b[2]
                      << ") with distance between them " << std::sqrt(distance2)
                      << "\n";
            return false;
        }
    }
    return true;
}

int main(int argc, char* argv[])
{
    auto const digits10 = std::numeric_limits<double>::digits10;
    auto const args = parseCommandLine(argc, argv);

    // Setup the standard output and error stream numerical formats.
    std::cout << std::scientific << std::setprecision(digits10);
    std::cerr << std::scientific << std::setprecision(digits10);

    auto meshes = readMeshes(args.vtk_input_a, args.vtk_input_b);

    if (args.meshcheck)
    {
        if (args.vtk_input_a == args.vtk_input_b)
        {
            std::cout << "Will not compare meshes from same input file.\n";
            return EXIT_SUCCESS;
        }
        if (!comparePoints(std::get<0>(meshes)->GetPoints(),
                           std::get<1>(meshes)->GetPoints(),
                           args.abs_err_thr * args.abs_err_thr))
        {
            std::cerr << "Error in mesh points' comparison occured.\n";
            return EXIT_FAILURE;
        }

        if (!compareCellTopology(std::get<0>(meshes)->GetCells(),
                                 std::get<1>(meshes)->GetCells()))
        {
            std::cerr << "Error in cells' topology comparison occured.\n";
            return EXIT_FAILURE;
        }
        return EXIT_SUCCESS;
    }

    // Read arrays from input file.
    bool read_successful;
    vtkSmartPointer<vtkDataArray> a;
    vtkSmartPointer<vtkDataArray> b;

    std::tie(read_successful, a, b) =
        readDataArraysFromMeshes(meshes, args.data_array_a, args.data_array_b);

    if (!read_successful)
        return EXIT_FAILURE;

    if (!args.quiet)
        std::cout << "Comparing data array `" << args.data_array_a
                  << "' from file `" << args.vtk_input_a << "' to data array `"
                  << args.data_array_b << "' from file `" << args.vtk_input_b
                  << "'.\n";

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
    std::vector<double> abs_err_norm_l1(num_components);
    std::vector<double> abs_err_norm_2_2(num_components);
    std::vector<double> abs_err_norm_max(num_components);

    // Relative error and norms.
    std::vector<double> rel_err_norm_l1(num_components);
    std::vector<double> rel_err_norm_2_2(num_components);
    std::vector<double> rel_err_norm_max(num_components);

    for (auto tuple_idx = 0; tuple_idx < num_tuples; ++tuple_idx)
    {
        for (auto component_idx = 0; component_idx < num_components;
             ++component_idx)
        {
            auto const a_comp = a->GetComponent(tuple_idx, component_idx);
            auto const b_comp = b->GetComponent(tuple_idx, component_idx);
            auto const abs_err = std::abs(a_comp - b_comp);

            abs_err_norm_l1[component_idx] += abs_err;
            abs_err_norm_2_2[component_idx] += abs_err * abs_err;
            abs_err_norm_max[component_idx] =
                std::max(abs_err_norm_max[component_idx], abs_err);

            // relative error and its norms:
            double rel_err;

            if (abs_err == 0.0)
            {
                rel_err = 0.0;
            }
            else if (a_comp == 0.0 || b_comp == 0.0)
            {
                rel_err = std::numeric_limits<double>::infinity();
            }
            else
            {
                rel_err =
                    abs_err / std::min(std::abs(a_comp), std::abs(b_comp));
            }

            rel_err_norm_l1[component_idx] += rel_err;
            rel_err_norm_2_2[component_idx] += rel_err * rel_err;
            rel_err_norm_max[component_idx] =
                std::max(rel_err_norm_max[component_idx], rel_err);

            if (abs_err > args.abs_err_thr && rel_err > args.rel_err_thr &&
                args.verbose)
            {
                std::cout << "tuple: " << std::setw(4) << tuple_idx
                          << "component: " << std::setw(2) << component_idx
                          << ": abs err = " << std::setw(digits10 + 7)
                          << abs_err
                          << ", rel err = " << std::setw(digits10 + 7)
                          << rel_err << "\n";
            }
        }
    }

    // Error information
    if (!args.quiet)
    {
        std::cout << "Computed difference between data arrays:\n";
        std::cout << "abs l1 norm      = " << abs_err_norm_l1 << "\n";
        std::cout << "abs l2-norm^2    = " << abs_err_norm_2_2 << "\n";

        // temporary squared norm vector for output.
        std::vector<double> abs_err_norm_2;
        std::transform(std::begin(abs_err_norm_2_2), std::end(abs_err_norm_2_2),
                       std::back_inserter(abs_err_norm_2),
                       [](double x) { return std::sqrt(x); });
        std::cout << "abs l2-norm      = " << abs_err_norm_2 << "\n";

        std::cout << "abs maximum norm = " << abs_err_norm_max << "\n";
        std::cout << "\n";

        std::cout << "rel l1 norm      = " << rel_err_norm_l1 << "\n";
        std::cout << "rel l2-norm^2    = " << rel_err_norm_2_2 << "\n";

        // temporary squared norm vector for output.
        std::vector<double> rel_err_norm_2;
        std::transform(std::begin(rel_err_norm_2_2), std::end(rel_err_norm_2_2),
                       std::back_inserter(rel_err_norm_2),
                       [](double x) { return std::sqrt(x); });
        std::cout << "rel l2-norm      = " << rel_err_norm_2_2 << "\n";

        std::cout << "rel maximum norm = " << rel_err_norm_max << "\n";
    }

    if (*std::max_element(abs_err_norm_max.begin(), abs_err_norm_max.end()) >
            args.abs_err_thr &&
        *std::max_element(rel_err_norm_max.begin(), rel_err_norm_max.end()) >
            args.rel_err_thr)
    {
        if (!args.quiet)
            std::cout << "Absolute and relative error (maximum norm) are larger"
                         " than the corresponding thresholds "
                      << args.abs_err_thr << " and " << args.rel_err_thr
                      << ".\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
