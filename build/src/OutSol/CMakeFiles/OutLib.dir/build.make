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
include src/OutSol/CMakeFiles/OutLib.dir/depend.make

# Include the progress variables for this target.
include src/OutSol/CMakeFiles/OutLib.dir/progress.make

# Include the compile flags for this target's objects.
include src/OutSol/CMakeFiles/OutLib.dir/flags.make

src/OutSol/CMakeFiles/OutLib.dir/OutSol.cc.obj: src/OutSol/CMakeFiles/OutLib.dir/flags.make
src/OutSol/CMakeFiles/OutLib.dir/OutSol.cc.obj: src/OutSol/CMakeFiles/OutLib.dir/includes_CXX.rsp
src/OutSol/CMakeFiles/OutLib.dir/OutSol.cc.obj: ../src/OutSol/OutSol.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=I:\code\git\PPP\PPPLib\build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/OutSol/CMakeFiles/OutLib.dir/OutSol.cc.obj"
	cd /d I:\code\git\PPP\PPPLib\build\src\OutSol && C:\mingw64\bin\x86_64-w64-mingw32-g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\OutLib.dir\OutSol.cc.obj -c I:\code\git\PPP\PPPLib\src\OutSol\OutSol.cc

src/OutSol/CMakeFiles/OutLib.dir/OutSol.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/OutLib.dir/OutSol.cc.i"
	cd /d I:\code\git\PPP\PPPLib\build\src\OutSol && C:\mingw64\bin\x86_64-w64-mingw32-g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E I:\code\git\PPP\PPPLib\src\OutSol\OutSol.cc > CMakeFiles\OutLib.dir\OutSol.cc.i

src/OutSol/CMakeFiles/OutLib.dir/OutSol.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/OutLib.dir/OutSol.cc.s"
	cd /d I:\code\git\PPP\PPPLib\build\src\OutSol && C:\mingw64\bin\x86_64-w64-mingw32-g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S I:\code\git\PPP\PPPLib\src\OutSol\OutSol.cc -o CMakeFiles\OutLib.dir\OutSol.cc.s

# Object files for target OutLib
OutLib_OBJECTS = \
"CMakeFiles/OutLib.dir/OutSol.cc.obj"

# External object files for target OutLib
OutLib_EXTERNAL_OBJECTS =

../bin/libOutLib.a: src/OutSol/CMakeFiles/OutLib.dir/OutSol.cc.obj
../bin/libOutLib.a: src/OutSol/CMakeFiles/OutLib.dir/build.make
../bin/libOutLib.a: src/OutSol/CMakeFiles/OutLib.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=I:\code\git\PPP\PPPLib\build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library ..\..\..\bin\libOutLib.a"
	cd /d I:\code\git\PPP\PPPLib\build\src\OutSol && $(CMAKE_COMMAND) -P CMakeFiles\OutLib.dir\cmake_clean_target.cmake
	cd /d I:\code\git\PPP\PPPLib\build\src\OutSol && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\OutLib.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/OutSol/CMakeFiles/OutLib.dir/build: ../bin/libOutLib.a

.PHONY : src/OutSol/CMakeFiles/OutLib.dir/build

src/OutSol/CMakeFiles/OutLib.dir/clean:
	cd /d I:\code\git\PPP\PPPLib\build\src\OutSol && $(CMAKE_COMMAND) -P CMakeFiles\OutLib.dir\cmake_clean.cmake
.PHONY : src/OutSol/CMakeFiles/OutLib.dir/clean

src/OutSol/CMakeFiles/OutLib.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" I:\code\git\PPP\PPPLib I:\code\git\PPP\PPPLib\src\OutSol I:\code\git\PPP\PPPLib\build I:\code\git\PPP\PPPLib\build\src\OutSol I:\code\git\PPP\PPPLib\build\src\OutSol\CMakeFiles\OutLib.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : src/OutSol/CMakeFiles/OutLib.dir/depend

