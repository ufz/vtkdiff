cmake_minimum_required(VERSION 3.12)
cmake_policy(SET CMP0074 NEW)

project(vtkdiff)

# C++ standard setup
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --- CPM ---
include(cmake/CPM.cmake)
if(GUIX_BUILD)
    find_path(_tclap_include UnlabeledValueArg.h PATH_SUFFIXES tclap REQUIRED)
    add_library(tclap INTERFACE IMPORTED)
    target_include_directories(
        tclap SYSTEM INTERFACE ${_tclap_include} ${_tclap_include}/..
    )
else()
    CPMFindPackage(
        NAME tclap
        GIT_REPOSITORY https://git.code.sf.net/p/tclap/code
        GIT_TAG 81b3d2a0c47895c22e9bb8c577f5ab521f76e5d2
        VERSION 1.4.0
        DOWNLOAD_ONLY YES
    )
    if(tclap_ADDED)
        add_library(tclap INTERFACE IMPORTED)
        target_include_directories(
            tclap SYSTEM INTERFACE ${tclap_SOURCE_DIR}/include
        )
    endif()
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

# Coverage (enabled by Debug-config)
if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug" AND ${PROJECT_SOURCE_DIR} STREQUAL
                                              ${CMAKE_SOURCE_DIR}
)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        set(CMAKE_CXX_FLAGS_DEBUG "-g -Og --coverage")
        if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
            string(APPEND CMAKE_CXX_FLAGS_DEBUG " -fprofile-abs-path")
        endif()
        set(CMAKE_EXE_LINKER_FLAGS_DEBUG "--coverage")
        set(CMAKE_SHARED_LINKER_FLAGS_DEBUG "--coverage")
        set(CMAKE_MODULE_LINKER_FLAGS_DEBUG "--coverage")
        message(STATUS "Setting coverage flags!")
        message(
            STATUS
                "Run executables and the run:\n     gcovr -r ${PROJECT_SOURCE_DIR}"
        )
    else()
        message(FATAL_ERROR "OGS_COVERAGE requires clang or gcc compiler!")
    endif()
endif()

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
