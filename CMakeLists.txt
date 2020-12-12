cmake_minimum_required(VERSION 3.11)
set(CMAKE_CXX_STANDARD 14)             # available options are [98, 11, 14, 17. 20]
project(mesh_averaging)

# OpenVDB only provides a module file (FindOpenVDB.cmake), so we use OpenVDB_DIR to find it
if (NOT DEFINED OpenVDB_DIR OR "${OpenVDB_DIR}" STREQUAL "OpenVDB_DIR-NOTFOUND")
  set(OpenVDB_DIR ${CMAKE_PREFIX_PATH}/lib/cmake/OpenVDB)
endif()
message(STATUS "OpenVDB_DIR set to ${OpenVDB_DIR} (must contain FindOpenVDB.cmake)")
file(TO_CMAKE_PATH ${OpenVDB_DIR} OpenVDB_DIR) # fixes path to use forward slashes on Windows
list(APPEND CMAKE_MODULE_PATH ${OpenVDB_DIR})
find_package(OpenVDB MODULE REQUIRED)



file(GLOB CPP_SRC
        "src/*.cpp"
        )
add_executable(mesh_averaging ${CPP_SRC})
target_include_directories(mesh_averaging PUBLIC src/)


target_link_libraries(mesh_averaging PUBLIC OpenVDB::openvdb)
