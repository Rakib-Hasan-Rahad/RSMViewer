cmake_minimum_required(VERSION 3.16)

project(RNAMotifScanX LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(MSVC)
  add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
  add_compile_options(/W3)
else()
  add_compile_options(-Wall -Wextra -Wno-deprecated-declarations)
endif()

set(RMSX_THREADPOOL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Libraries/threadpool")

find_package(Boost 1.54 REQUIRED COMPONENTS program_options thread iostreams)
find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)

set(RMSX_CORE_SOURCES
  structural_motif.cc
  annotated_structure.cc
  scan_structure.cc
  motif_graph_matching.cc
  simulation.cc
  parameters.cc
  map_PDB_index.cc
)

set(RMSX_COMMON_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${RMSX_THREADPOOL_DIR}
  ${Boost_INCLUDE_DIRS}
)

set(RMSX_COMMON_LIBS
  ${Boost_LIBRARIES}
  Threads::Threads
  ZLIB::ZLIB
)

function(rmsx_add_tool target_name entry_source)
  add_executable(${target_name} ${entry_source} ${RMSX_CORE_SOURCES})
  target_include_directories(${target_name} PRIVATE ${RMSX_COMMON_INCLUDE_DIRS})
  target_link_libraries(${target_name} PRIVATE ${RMSX_COMMON_LIBS})
  if(UNIX AND NOT APPLE)
    target_link_libraries(${target_name} PRIVATE rt)
  endif()
endfunction()

rmsx_add_tool(scan main_scan.cc)
rmsx_add_tool(cut main_cut.cc)
rmsx_add_tool(align main_align.cc)
