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
include src/CmnFunc/CMakeFiles/CmnLib.dir/depend.make

# Include the progress variables for this target.
include src/CmnFunc/CMakeFiles/CmnLib.dir/progress.make

# Include the compile flags for this target's objects.
include src/CmnFunc/CMakeFiles/CmnLib.dir/flags.make

src/CmnFunc/CMakeFiles/CmnLib.dir/CmnFunc.cc.obj: src/CmnFunc/CMakeFiles/CmnLib.dir/flags.make
src/CmnFunc/CMakeFiles/CmnLib.dir/CmnFunc.cc.obj: src/CmnFunc/CMakeFiles/CmnLib.dir/includes_CXX.rsp
src/CmnFunc/CMakeFiles/CmnLib.dir/CmnFunc.cc.obj: ../src/CmnFunc/CmnFunc.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=I:\code\git\PPP\PPPLib\build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CmnFunc/CMakeFiles/CmnLib.dir/CmnFunc.cc.obj"
	cd /d I:\code\git\PPP\PPPLib\build\src\CmnFunc && C:\mingw64\bin\x86_64-w64-mingw32-g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\CmnLib.dir\CmnFunc.cc.obj -c I:\code\git\PPP\PPPLib\src\CmnFunc\CmnFunc.cc

src/CmnFunc/CMakeFiles/CmnLib.dir/CmnFunc.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CmnLib.dir/CmnFunc.cc.i"
	cd /d I:\code\git\PPP\PPPLib\build\src\CmnFunc && C:\mingw64\bin\x86_64-w64-mingw32-g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E I:\code\git\PPP\PPPLib\src\CmnFunc\CmnFunc.cc > CMakeFiles\CmnLib.dir\CmnFunc.cc.i

src/CmnFunc/CMakeFiles/CmnLib.dir/CmnFunc.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CmnLib.dir/CmnFunc.cc.s"
	cd /d I:\code\git\PPP\PPPLib\build\src\CmnFunc && C:\mingw64\bin\x86_64-w64-mingw32-g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S I:\code\git\PPP\PPPLib\src\CmnFunc\CmnFunc.cc -o CMakeFiles\CmnLib.dir\CmnFunc.cc.s

# Object files for target CmnLib
CmnLib_OBJECTS = \
"CMakeFiles/CmnLib.dir/CmnFunc.cc.obj"

# External object files for target CmnLib
CmnLib_EXTERNAL_OBJECTS =

../bin/libCmnLib.a: src/CmnFunc/CMakeFiles/CmnLib.dir/CmnFunc.cc.obj
../bin/libCmnLib.a: src/CmnFunc/CMakeFiles/CmnLib.dir/build.make
../bin/libCmnLib.a: src/CmnFunc/CMakeFiles/CmnLib.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=I:\code\git\PPP\PPPLib\build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library ..\..\..\bin\libCmnLib.a"
	cd /d I:\code\git\PPP\PPPLib\build\src\CmnFunc && $(CMAKE_COMMAND) -P CMakeFiles\CmnLib.dir\cmake_clean_target.cmake
	cd /d I:\code\git\PPP\PPPLib\build\src\CmnFunc && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\CmnLib.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CmnFunc/CMakeFiles/CmnLib.dir/build: ../bin/libCmnLib.a

.PHONY : src/CmnFunc/CMakeFiles/CmnLib.dir/build

src/CmnFunc/CMakeFiles/CmnLib.dir/clean:
	cd /d I:\code\git\PPP\PPPLib\build\src\CmnFunc && $(CMAKE_COMMAND) -P CMakeFiles\CmnLib.dir\cmake_clean.cmake
.PHONY : src/CmnFunc/CMakeFiles/CmnLib.dir/clean

src/CmnFunc/CMakeFiles/CmnLib.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" I:\code\git\PPP\PPPLib I:\code\git\PPP\PPPLib\src\CmnFunc I:\code\git\PPP\PPPLib\build I:\code\git\PPP\PPPLib\build\src\CmnFunc I:\code\git\PPP\PPPLib\build\src\CmnFunc\CMakeFiles\CmnLib.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : src/CmnFunc/CMakeFiles/CmnLib.dir/depend

