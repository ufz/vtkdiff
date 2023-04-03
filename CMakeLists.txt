cmake_minimum_required(VERSION 3.12)
cmake_policy(SET CMP0074 NEW)

project(vtkdiff)

# C++ standard setup
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --- CPM ---
include(cmake/CPM.cmake)
CPMFindPackage(
    NAME tclap
    GITHUB_REPOSITORY ufz/tclap
    VERSION 1.2.4
    GIT_TAG 098dd0fe07a31618f3c2a9f8727bb01c8c5d61e2
    DOWNLOAD_ONLY YES
)
if(tclap_ADDED)
    add_library(tclap INTERFACE IMPORTED)
    target_include_directories(
        tclap SYSTEM INTERFACE ${tclap_SOURCE_DIR}/include
    )
endif()

set(VTK_COMPONENTS IOXML)
find_package(VTK 9 QUIET COMPONENTS ${VTK_COMPONENTS})

if(NOT VTK_FOUND)
    list(
        APPEND
        VTK_OPTIONS
        "BUILD_SHARED_LIBS OFF"
        "BUILD_TESTING OFF"
        "VTK_BUILD_EXAMPLES OFF"
        "VTK_BUILD_TESTING OFF"
        "VTK_ENABLE_WRAPPING OFF"
        "VTK_GROUP_ENABLE_Rendering DONT_WANT"
        "VTK_GROUP_ENABLE_StandAlone DONT_WANT"
        "VTK_USE_64BIT_IDS ON"
        "VTK_MODULE_ENABLE_VTK_IOXML YES"
    )

    CPMAddPackage(
        NAME VTK
        GITHUB_REPOSITORY kitware/vtk
        VERSION 9.1.0
        OPTIONS ${VTK_OPTIONS}
        EXCLUDE_FROM_ALL YES GIT_SUBMODULES "" # Disable submodules
    )
endif()
# --- CPM end ---

add_executable(vtkdiff vtkdiff.cpp)
target_link_libraries(vtkdiff tclap VTK::IOXML)

# Set compiler helper variables
if(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
    set(COMPILER_IS_CLANG TRUE)
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    set(COMPILER_IS_GCC TRUE)
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
    set(COMPILER_IS_INTEL TRUE)
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "MSVC")
    set(COMPILER_IS_MSVC TRUE)
endif() # CMAKE_CXX_COMPILER_ID

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(COMPILER_IS_GCC OR COMPILER_IS_CLANG)
    target_compile_options(vtkdiff PUBLIC -Wall -Wextra)
endif()

if(COMPILER_IS_INTEL)
    target_compile_options(vtkdiff PUBLIC -xHOST -no-prec-div)
endif()

if(COMPILER_IS_MSVC)
    target_compile_options(
        vtkdiff
        PUBLIC -D_CRT_SECURE_NO_WARNINGS
               -D_CRT_NONSTDC_NO_WARNINGS
               -D_CRT_XNONSTDC_NO_WARNINGS
               -D__restrict__=__restrict # this fixes #5
               -DNOMINMAX # This fixes compile errors with
                          # std::numeric_limits<T>::min() / max()
               /W3
               /wd4290
               /wd4267
               /wd4996
    )
endif()
