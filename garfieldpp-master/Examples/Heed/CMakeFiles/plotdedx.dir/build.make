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
CMAKE_SOURCE_DIR = /home/swuupie/PET/garfieldpp-master/Examples/Heed

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/swuupie/PET/garfieldpp-master/Examples/Heed

# Include any dependencies generated for this target.
include CMakeFiles/plotdedx.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/plotdedx.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/plotdedx.dir/flags.make

CMakeFiles/plotdedx.dir/plotdedx.C.o: CMakeFiles/plotdedx.dir/flags.make
CMakeFiles/plotdedx.dir/plotdedx.C.o: plotdedx.C
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/swuupie/PET/garfieldpp-master/Examples/Heed/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/plotdedx.dir/plotdedx.C.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/plotdedx.dir/plotdedx.C.o -c /home/swuupie/PET/garfieldpp-master/Examples/Heed/plotdedx.C

CMakeFiles/plotdedx.dir/plotdedx.C.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/plotdedx.dir/plotdedx.C.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/swuupie/PET/garfieldpp-master/Examples/Heed/plotdedx.C > CMakeFiles/plotdedx.dir/plotdedx.C.i

CMakeFiles/plotdedx.dir/plotdedx.C.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/plotdedx.dir/plotdedx.C.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/swuupie/PET/garfieldpp-master/Examples/Heed/plotdedx.C -o CMakeFiles/plotdedx.dir/plotdedx.C.s

# Object files for target plotdedx
plotdedx_OBJECTS = \
"CMakeFiles/plotdedx.dir/plotdedx.C.o"

# External object files for target plotdedx
plotdedx_EXTERNAL_OBJECTS =

plotdedx: CMakeFiles/plotdedx.dir/plotdedx.C.o
plotdedx: CMakeFiles/plotdedx.dir/build.make
plotdedx: /home/swuupie/garfieldpp/install/lib/libGarfield.so.0.3.0
plotdedx: /home/swuupie/root/lib/libGdml.so.6.24.02
plotdedx: /home/swuupie/root/lib/libGeom.so.6.24.02
plotdedx: /home/swuupie/root/lib/libXMLIO.so.6.24.02
plotdedx: /home/swuupie/root/lib/libGraf3d.so.6.24.02
plotdedx: /home/swuupie/root/lib/libGpad.so.6.24.02
plotdedx: /home/swuupie/root/lib/libGraf.so.6.24.02
plotdedx: /home/swuupie/root/lib/libHist.so.6.24.02
plotdedx: /home/swuupie/root/lib/libMatrix.so.6.24.02
plotdedx: /home/swuupie/root/lib/libMathCore.so.6.24.02
plotdedx: /home/swuupie/root/lib/libImt.so.6.24.02
plotdedx: /home/swuupie/root/lib/libMultiProc.so.6.24.02
plotdedx: /home/swuupie/root/lib/libNet.so.6.24.02
plotdedx: /home/swuupie/root/lib/libRIO.so.6.24.02
plotdedx: /home/swuupie/root/lib/libThread.so.6.24.02
plotdedx: /home/swuupie/root/lib/libCore.so.6.24.02
plotdedx: /usr/lib/x86_64-linux-gnu/libgsl.so
plotdedx: /usr/lib/x86_64-linux-gnu/libgslcblas.so
plotdedx: /home/swuupie/garfieldpp/install/lib/libmagboltz.so.11
plotdedx: CMakeFiles/plotdedx.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/swuupie/PET/garfieldpp-master/Examples/Heed/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable plotdedx"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/plotdedx.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/plotdedx.dir/build: plotdedx

.PHONY : CMakeFiles/plotdedx.dir/build

CMakeFiles/plotdedx.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/plotdedx.dir/cmake_clean.cmake
.PHONY : CMakeFiles/plotdedx.dir/clean

CMakeFiles/plotdedx.dir/depend:
	cd /home/swuupie/PET/garfieldpp-master/Examples/Heed && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/swuupie/PET/garfieldpp-master/Examples/Heed /home/swuupie/PET/garfieldpp-master/Examples/Heed /home/swuupie/PET/garfieldpp-master/Examples/Heed /home/swuupie/PET/garfieldpp-master/Examples/Heed /home/swuupie/PET/garfieldpp-master/Examples/Heed/CMakeFiles/plotdedx.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/plotdedx.dir/depend
