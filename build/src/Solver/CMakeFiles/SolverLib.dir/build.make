# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.19

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = C:\cmake-3.19.6-win64-x64\bin\cmake.exe

# The command to remove a file.
RM = C:\cmake-3.19.6-win64-x64\bin\cmake.exe -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = I:\code\git\PPP\PPPLib

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = I:\code\git\PPP\PPPLib\build

# Include any dependencies generated for this target.
include src/Solver/CMakeFiles/SolverLib.dir/depend.make

# Include the progress variables for this target.
include src/Solver/CMakeFiles/SolverLib.dir/progress.make

# Include the compile flags for this target's objects.
include src/Solver/CMakeFiles/SolverLib.dir/flags.make

src/Solver/CMakeFiles/SolverLib.dir/Solver.cc.obj: src/Solver/CMakeFiles/SolverLib.dir/flags.make
src/Solver/CMakeFiles/SolverLib.dir/Solver.cc.obj: src/Solver/CMakeFiles/SolverLib.dir/includes_CXX.rsp
src/Solver/CMakeFiles/SolverLib.dir/Solver.cc.obj: ../src/Solver/Solver.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=I:\code\git\PPP\PPPLib\build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/Solver/CMakeFiles/SolverLib.dir/Solver.cc.obj"
	cd /d I:\code\git\PPP\PPPLib\build\src\Solver && C:\mingw64\bin\x86_64-w64-mingw32-g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\SolverLib.dir\Solver.cc.obj -c I:\code\git\PPP\PPPLib\src\Solver\Solver.cc

src/Solver/CMakeFiles/SolverLib.dir/Solver.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SolverLib.dir/Solver.cc.i"
	cd /d I:\code\git\PPP\PPPLib\build\src\Solver && C:\mingw64\bin\x86_64-w64-mingw32-g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E I:\code\git\PPP\PPPLib\src\Solver\Solver.cc > CMakeFiles\SolverLib.dir\Solver.cc.i

src/Solver/CMakeFiles/SolverLib.dir/Solver.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SolverLib.dir/Solver.cc.s"
	cd /d I:\code\git\PPP\PPPLib\build\src\Solver && C:\mingw64\bin\x86_64-w64-mingw32-g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S I:\code\git\PPP\PPPLib\src\Solver\Solver.cc -o CMakeFiles\SolverLib.dir\Solver.cc.s

# Object files for target SolverLib
SolverLib_OBJECTS = \
"CMakeFiles/SolverLib.dir/Solver.cc.obj"

# External object files for target SolverLib
SolverLib_EXTERNAL_OBJECTS =

../bin/libSolverLib.a: src/Solver/CMakeFiles/SolverLib.dir/Solver.cc.obj
../bin/libSolverLib.a: src/Solver/CMakeFiles/SolverLib.dir/build.make
../bin/libSolverLib.a: src/Solver/CMakeFiles/SolverLib.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=I:\code\git\PPP\PPPLib\build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library ..\..\..\bin\libSolverLib.a"
	cd /d I:\code\git\PPP\PPPLib\build\src\Solver && $(CMAKE_COMMAND) -P CMakeFiles\SolverLib.dir\cmake_clean_target.cmake
	cd /d I:\code\git\PPP\PPPLib\build\src\Solver && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\SolverLib.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/Solver/CMakeFiles/SolverLib.dir/build: ../bin/libSolverLib.a

.PHONY : src/Solver/CMakeFiles/SolverLib.dir/build

src/Solver/CMakeFiles/SolverLib.dir/clean:
	cd /d I:\code\git\PPP\PPPLib\build\src\Solver && $(CMAKE_COMMAND) -P CMakeFiles\SolverLib.dir\cmake_clean.cmake
.PHONY : src/Solver/CMakeFiles/SolverLib.dir/clean

src/Solver/CMakeFiles/SolverLib.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" I:\code\git\PPP\PPPLib I:\code\git\PPP\PPPLib\src\Solver I:\code\git\PPP\PPPLib\build I:\code\git\PPP\PPPLib\build\src\Solver I:\code\git\PPP\PPPLib\build\src\Solver\CMakeFiles\SolverLib.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : src/Solver/CMakeFiles/SolverLib.dir/depend

