cmake_minimum_required(VERSION 3.23)
project(QLDPC_SC_Construction)

set(CMAKE_CXX_STANDARD 14)

add_executable(QLDPC_SC_Construction sc_generate_code.cpp sc_qldpc_optimization.h common_code.h MD.h SC.h tools.h)
