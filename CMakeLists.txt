# CMakeList.txt: Algo_Distrib_Robust 的 CMake 项目，在此处包括源代码并定义
# 项目特定的逻辑。
#
cmake_minimum_required (VERSION 3.19)

project ("tp2")

# 将源代码添加到此项目的可执行文件。
add_executable (tp2 "ex01 fibbonacci.cpp" )

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET tp2 PROPERTY CXX_STANDARD 20)
endif()

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
  message("OPENMP FOUND")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# TODO: 如有需要，请添加测试并安装目标。
