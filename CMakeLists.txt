cmake_minimum_required(VERSION 3.1)

if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    find_package(VTK REQUIRED)
    include(${VTK_USE_FILE})
endif()

include_directories(SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/tclap/include)

add_executable(vtkdiff vtkdiff.cpp)

target_link_libraries(vtkdiff ${VTK_LIBRARIES})

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
