cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(Silicon)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()

add_executable(working working.C)
target_link_libraries(working Garfield::Garfield)
