# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(neBEM)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()
# ########################################################################################

# ---Build executables-----------------------------------------------------------
add_executable(RPCbyneBEM RPCbyneBEM.C)
target_link_libraries(RPCbyneBEM Garfield::Garfield)

add_executable(MWPCbyneBEM MWPCbyneBEM.C)
target_link_libraries(MWPCbyneBEM Garfield::Garfield)

add_executable(BiconicalGemByneBEM BiconicalGemByneBEM.C)
target_link_libraries(BiconicalGemByneBEM Garfield::Garfield)

add_executable(parallelPlates parallelPlates.C)
target_link_libraries(parallelPlates Garfield::Garfield)

add_executable(crossingWires crossingWires.C)
target_link_libraries(crossingWires Garfield::Garfield)

add_executable(parallelPlatesWeightingField parallelPlatesWeightingField.C)
target_link_libraries(parallelPlatesWeightingField Garfield::Garfield)
