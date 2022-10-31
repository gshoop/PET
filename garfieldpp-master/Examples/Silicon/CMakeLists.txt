cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(Silicon)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()

# ---Build executable------------------------------------------------------------
add_executable(planar planar.C)
target_link_libraries(planar Garfield::Garfield)

add_executable(planar_movie planar_movie.C)
target_link_libraries(planar_movie Garfield::Garfield)

add_executable(example ${CMAKE_CURRENT_SOURCE_DIR}/ProjectCZT/example.C)
target_link_libraries(example Garfield::Garfield)

add_executable(cdTe ${CMAKE_CURRENT_SOURCE_DIR}/ProjectCZT/exampleCdTe.C)
target_link_libraries(cdTe Garfield::Garfield)