cmake_minimum_required(VERSION 3.1)
project(vtkdiff)

# --- CPM ---
include(cmake/CPM.cmake)
CPMAddPackage(
    NAME tclap
    GITHUB_REPOSITORY ufz/tclap
    GIT_TAG 03abc3a3327214137c6ffd5b9a6efe23f0927cc2
    DOWNLOAD_ONLY YES
)
if(tclap_ADDED)
    add_library(tclap INTERFACE IMPORTED)
    target_include_directories(
        tclap SYSTEM INTERFACE ${tclap_SOURCE_DIR}/include
    )
endif()

set(VTK_COMPONENTS vtkIOXML)
find_package(VTK 8.2.0 QUIET COMPONENTS ${VTK_COMPONENTS})

if(VTK_FOUND)
    include(${VTK_USE_FILE})
else()
    list(
        APPEND
        VTK_OPTIONS
        "BUILD_SHARED_LIBS OFF"
        "BUILD_TESTING OFF"
        "VTK_BUILD_EXAMPLES OFF"
        "VTK_BUILD_TESTING OFF"
        "VTK_ENABLE_WRAPPING OFF"
        "VTK_Group_Rendering OFF"
        "VTK_Group_StandAlone OFF"
        "VTK_USE_64BIT_IDS ON"
    )
    foreach(comp ${VTK_COMPONENTS})
        list(APPEND VTK_OPTIONS "Module_${comp} ON")
    endforeach()

    # Workaround for configuration error in [vtk]/CMake/vtkGroups.cmake:43
    set(VTK_Group_Rendering OFF CACHE BOOL "")
    set(VTK_Group_StandAlone OFF CACHE BOOL "")

    CPMAddPackage(
        NAME VTK
        GITHUB_REPOSITORY kitware/vtk
        VERSION 8.2.0
        OPTIONS ${VTK_OPTIONS}
        EXCLUDE_FROM_ALL YES GIT_SUBMODULES "" # Disable submodules
    )
    include(${VTK_BINARY_DIR}/VTKConfig.cmake)
endif()
# --- CPM end ---

add_executable(vtkdiff vtkdiff.cpp)
target_include_directories(vtkdiff SYSTEM PRIVATE ${VTK_INCLUDE_DIRS})
target_link_libraries(vtkdiff tclap ${VTK_LIBRARIES})

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
        PUBLIC _CRT_SECURE_NO_WARNINGS
               _CRT_NONSTDC_NO_WARNINGS
               _CRT_XNONSTDC_NO_WARNINGS
               __restrict__=__restrict # this fixes #5
               NOMINMAX # This fixes compile errors with
                        # std::numeric_limits<T>::min() / max()
               /W3
               /wd4290
               /wd4267
               /wd4996
    )
endif()
