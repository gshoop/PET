# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /mnt/c/Users/greys/Documents/PET/garfieldpp-master/Examples/Silicon

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/c/Users/greys/Documents/PET/garfieldpp-master/Examples/Silicon

# Include any dependencies generated for this target.
include CMakeFiles/planar_movie.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/planar_movie.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/planar_movie.dir/flags.make

CMakeFiles/planar_movie.dir/planar_movie.C.o: CMakeFiles/planar_movie.dir/flags.make
CMakeFiles/planar_movie.dir/planar_movie.C.o: planar_movie.C
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/Users/greys/Documents/PET/garfieldpp-master/Examples/Silicon/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/planar_movie.dir/planar_movie.C.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/planar_movie.dir/planar_movie.C.o -c /mnt/c/Users/greys/Documents/PET/garfieldpp-master/Examples/Silicon/planar_movie.C

CMakeFiles/planar_movie.dir/planar_movie.C.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/planar_movie.dir/planar_movie.C.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/Users/greys/Documents/PET/garfieldpp-master/Examples/Silicon/planar_movie.C > CMakeFiles/planar_movie.dir/planar_movie.C.i

CMakeFiles/planar_movie.dir/planar_movie.C.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/planar_movie.dir/planar_movie.C.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/Users/greys/Documents/PET/garfieldpp-master/Examples/Silicon/planar_movie.C -o CMakeFiles/planar_movie.dir/planar_movie.C.s

# Object files for target planar_movie
planar_movie_OBJECTS = \
"CMakeFiles/planar_movie.dir/planar_movie.C.o"

# External object files for target planar_movie
planar_movie_EXTERNAL_OBJECTS =

planar_movie: CMakeFiles/planar_movie.dir/planar_movie.C.o
planar_movie: CMakeFiles/planar_movie.dir/build.make
planar_movie: /home/swuupie/garfieldpp/install/lib/libGarfield.so.0.3.0
planar_movie: /home/swuupie/root/lib/libGdml.so.6.24.02
planar_movie: /home/swuupie/root/lib/libGeom.so.6.24.02
planar_movie: /home/swuupie/root/lib/libXMLIO.so.6.24.02
planar_movie: /home/swuupie/root/lib/libGraf3d.so.6.24.02
planar_movie: /home/swuupie/root/lib/libGpad.so.6.24.02
planar_movie: /home/swuupie/root/lib/libGraf.so.6.24.02
planar_movie: /home/swuupie/root/lib/libHist.so.6.24.02
planar_movie: /home/swuupie/root/lib/libMatrix.so.6.24.02
planar_movie: /home/swuupie/root/lib/libMathCore.so.6.24.02
planar_movie: /home/swuupie/root/lib/libImt.so.6.24.02
planar_movie: /home/swuupie/root/lib/libMultiProc.so.6.24.02
planar_movie: /home/swuupie/root/lib/libNet.so.6.24.02
planar_movie: /home/swuupie/root/lib/libRIO.so.6.24.02
planar_movie: /home/swuupie/root/lib/libThread.so.6.24.02
planar_movie: /home/swuupie/root/lib/libCore.so.6.24.02
planar_movie: /usr/lib/x86_64-linux-gnu/libgsl.so
planar_movie: /usr/lib/x86_64-linux-gnu/libgslcblas.so
planar_movie: /home/swuupie/garfieldpp/install/lib/libmagboltz.so.11
planar_movie: CMakeFiles/planar_movie.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/c/Users/greys/Documents/PET/garfieldpp-master/Examples/Silicon/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable planar_movie"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/planar_movie.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/planar_movie.dir/build: planar_movie

.PHONY : CMakeFiles/planar_movie.dir/build

CMakeFiles/planar_movie.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/planar_movie.dir/cmake_clean.cmake
.PHONY : CMakeFiles/planar_movie.dir/clean

CMakeFiles/planar_movie.dir/depend:
	cd /mnt/c/Users/greys/Documents/PET/garfieldpp-master/Examples/Silicon && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/c/Users/greys/Documents/PET/garfieldpp-master/Examples/Silicon /mnt/c/Users/greys/Documents/PET/garfieldpp-master/Examples/Silicon /mnt/c/Users/greys/Documents/PET/garfieldpp-master/Examples/Silicon /mnt/c/Users/greys/Documents/PET/garfieldpp-master/Examples/Silicon /mnt/c/Users/greys/Documents/PET/garfieldpp-master/Examples/Silicon/CMakeFiles/planar_movie.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/planar_movie.dir/depend
