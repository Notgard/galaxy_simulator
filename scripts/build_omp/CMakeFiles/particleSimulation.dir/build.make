# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.30

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /gpfs/app/2025/spack_install/linux-rhel9-zen4/linux-rhel9-zen4/aocc-5.0.0/cmake-3.30.5-py755ebrh2kf43z2wzqblkcztqcsiied/bin/cmake

# The command to remove a file.
RM = /gpfs/app/2025/spack_install/linux-rhel9-zen4/linux-rhel9-zen4/aocc-5.0.0/cmake-3.30.5-py755ebrh2kf43z2wzqblkcztqcsiied/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/garouba/galaxy_simulator

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/garouba/galaxy_simulator/scripts/build_omp

# Include any dependencies generated for this target.
include CMakeFiles/particleSimulation.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/particleSimulation.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/particleSimulation.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/particleSimulation.dir/flags.make

CMakeFiles/particleSimulation.dir/src/Simulation.cpp.o: CMakeFiles/particleSimulation.dir/flags.make
CMakeFiles/particleSimulation.dir/src/Simulation.cpp.o: /home/garouba/galaxy_simulator/src/Simulation.cpp
CMakeFiles/particleSimulation.dir/src/Simulation.cpp.o: CMakeFiles/particleSimulation.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/garouba/galaxy_simulator/scripts/build_omp/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/particleSimulation.dir/src/Simulation.cpp.o"
	/opt/softwares/amd/aocc-compiler-5.0.0/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/particleSimulation.dir/src/Simulation.cpp.o -MF CMakeFiles/particleSimulation.dir/src/Simulation.cpp.o.d -o CMakeFiles/particleSimulation.dir/src/Simulation.cpp.o -c /home/garouba/galaxy_simulator/src/Simulation.cpp

CMakeFiles/particleSimulation.dir/src/Simulation.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/particleSimulation.dir/src/Simulation.cpp.i"
	/opt/softwares/amd/aocc-compiler-5.0.0/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/garouba/galaxy_simulator/src/Simulation.cpp > CMakeFiles/particleSimulation.dir/src/Simulation.cpp.i

CMakeFiles/particleSimulation.dir/src/Simulation.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/particleSimulation.dir/src/Simulation.cpp.s"
	/opt/softwares/amd/aocc-compiler-5.0.0/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/garouba/galaxy_simulator/src/Simulation.cpp -o CMakeFiles/particleSimulation.dir/src/Simulation.cpp.s

CMakeFiles/particleSimulation.dir/src/main.cpp.o: CMakeFiles/particleSimulation.dir/flags.make
CMakeFiles/particleSimulation.dir/src/main.cpp.o: /home/garouba/galaxy_simulator/src/main.cpp
CMakeFiles/particleSimulation.dir/src/main.cpp.o: CMakeFiles/particleSimulation.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/garouba/galaxy_simulator/scripts/build_omp/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/particleSimulation.dir/src/main.cpp.o"
	/opt/softwares/amd/aocc-compiler-5.0.0/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/particleSimulation.dir/src/main.cpp.o -MF CMakeFiles/particleSimulation.dir/src/main.cpp.o.d -o CMakeFiles/particleSimulation.dir/src/main.cpp.o -c /home/garouba/galaxy_simulator/src/main.cpp

CMakeFiles/particleSimulation.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/particleSimulation.dir/src/main.cpp.i"
	/opt/softwares/amd/aocc-compiler-5.0.0/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/garouba/galaxy_simulator/src/main.cpp > CMakeFiles/particleSimulation.dir/src/main.cpp.i

CMakeFiles/particleSimulation.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/particleSimulation.dir/src/main.cpp.s"
	/opt/softwares/amd/aocc-compiler-5.0.0/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/garouba/galaxy_simulator/src/main.cpp -o CMakeFiles/particleSimulation.dir/src/main.cpp.s

# Object files for target particleSimulation
particleSimulation_OBJECTS = \
"CMakeFiles/particleSimulation.dir/src/Simulation.cpp.o" \
"CMakeFiles/particleSimulation.dir/src/main.cpp.o"

# External object files for target particleSimulation
particleSimulation_EXTERNAL_OBJECTS =

particleSimulation: CMakeFiles/particleSimulation.dir/src/Simulation.cpp.o
particleSimulation: CMakeFiles/particleSimulation.dir/src/main.cpp.o
particleSimulation: CMakeFiles/particleSimulation.dir/build.make
particleSimulation: /opt/softwares/amd/aocc-compiler-5.0.0/lib/libomp.so
particleSimulation: /usr/lib64/libpthread.a
particleSimulation: CMakeFiles/particleSimulation.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/garouba/galaxy_simulator/scripts/build_omp/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable particleSimulation"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/particleSimulation.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/particleSimulation.dir/build: particleSimulation
.PHONY : CMakeFiles/particleSimulation.dir/build

CMakeFiles/particleSimulation.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/particleSimulation.dir/cmake_clean.cmake
.PHONY : CMakeFiles/particleSimulation.dir/clean

CMakeFiles/particleSimulation.dir/depend:
	cd /home/garouba/galaxy_simulator/scripts/build_omp && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/garouba/galaxy_simulator /home/garouba/galaxy_simulator /home/garouba/galaxy_simulator/scripts/build_omp /home/garouba/galaxy_simulator/scripts/build_omp /home/garouba/galaxy_simulator/scripts/build_omp/CMakeFiles/particleSimulation.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/particleSimulation.dir/depend

