cmake_minimum_required(VERSION 3.2)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})


include_directories(SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/tclap/include)

add_executable(vtkdiff vtkdiff.cpp)

if(VTK_LIBRARIES)
  target_link_libraries(vtkdiff ${VTK_LIBRARIES})
else()
  target_link_libraries(vtkdiff vtkHybrid)
endif()

set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wextra -Weffc++")
