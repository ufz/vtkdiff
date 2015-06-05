cmake_minimum_required(VERSION 3.1)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})


include_directories(SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/tclap/include)

add_executable(vtkdiff vtkdiff.cpp)

if(VTK_LIBRARIES)
  target_link_libraries(vtkdiff ${VTK_LIBRARIES})
else()
  target_link_libraries(vtkdiff vtkHybrid)
endif()

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
    set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wextra")
endif()

if(COMPILER_IS_INTEL)
	set(CMAKE_CXX_FLAGS "-std=c++11 -xHOST -no-prec-div")
endif()

if (COMPILER_IS_MSVC)
    add_definitions(
        -D_CRT_SECURE_NO_WARNINGS
        -D_CRT_NONSTDC_NO_WARNINGS
        -D_CRT_XNONSTDC_NO_WARNINGS
        -D__restrict__=__restrict   # this fixes #5
        -DNOMINMAX # This fixes compile errors with std::numeric_limits<T>::min() / max()
    )
    # Sets warning level 3 and ignores some warnings
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3 /wd4290 /wd4267 /wd4996")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} /ZI /Od /Ob0")

    DisableCompilerFlag(DEBUG /RTC1)
endif ()
