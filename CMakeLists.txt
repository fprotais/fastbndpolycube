cmake_minimum_required(VERSION 2.8)
project(Polycube)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)



add_subdirectory(ultimaille)
include_directories(ultimaille ultimaille/ext)

add_subdirectory(amgcl)

find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if (WIN32) 
	set(BUILD_SHARED_LIBS OFF)
endif()


add_executable(fastpolycube main.cpp)
target_link_libraries(fastpolycube ultimaille OpenNL_psm)

add_executable(fastpolycube_amgcl main_w_amgcl.cpp)
target_link_libraries(fastpolycube_amgcl ultimaille OpenNL_psm)

target_link_libraries(fastpolycube_amgcl amgcl)

if(OpenMP_CXX_FOUND)
    target_link_libraries(fastpolycube OpenMP::OpenMP_CXX)
	target_link_libraries(fastpolycube_amgcl OpenMP::OpenMP_CXX)
endif()

