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
include exe/CMakeFiles/PPPLibMain.dir/depend.make

# Include the progress variables for this target.
include exe/CMakeFiles/PPPLibMain.dir/progress.make

# Include the compile flags for this target's objects.
include exe/CMakeFiles/PPPLibMain.dir/flags.make

exe/CMakeFiles/PPPLibMain.dir/PPPLibMain.cc.obj: exe/CMakeFiles/PPPLibMain.dir/flags.make
exe/CMakeFiles/PPPLibMain.dir/PPPLibMain.cc.obj: exe/CMakeFiles/PPPLibMain.dir/includes_CXX.rsp
exe/CMakeFiles/PPPLibMain.dir/PPPLibMain.cc.obj: ../exe/PPPLibMain.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=I:\code\git\PPP\PPPLib\build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object exe/CMakeFiles/PPPLibMain.dir/PPPLibMain.cc.obj"
	cd /d I:\code\git\PPP\PPPLib\build\exe && C:\mingw64\bin\x86_64-w64-mingw32-g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\PPPLibMain.dir\PPPLibMain.cc.obj -c I:\code\git\PPP\PPPLib\exe\PPPLibMain.cc

exe/CMakeFiles/PPPLibMain.dir/PPPLibMain.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PPPLibMain.dir/PPPLibMain.cc.i"
	cd /d I:\code\git\PPP\PPPLib\build\exe && C:\mingw64\bin\x86_64-w64-mingw32-g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E I:\code\git\PPP\PPPLib\exe\PPPLibMain.cc > CMakeFiles\PPPLibMain.dir\PPPLibMain.cc.i

exe/CMakeFiles/PPPLibMain.dir/PPPLibMain.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PPPLibMain.dir/PPPLibMain.cc.s"
	cd /d I:\code\git\PPP\PPPLib\build\exe && C:\mingw64\bin\x86_64-w64-mingw32-g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S I:\code\git\PPP\PPPLib\exe\PPPLibMain.cc -o CMakeFiles\PPPLibMain.dir\PPPLibMain.cc.s

# Object files for target PPPLibMain
PPPLibMain_OBJECTS = \
"CMakeFiles/PPPLibMain.dir/PPPLibMain.cc.obj"

# External object files for target PPPLibMain
PPPLibMain_EXTERNAL_OBJECTS =

../bin/PPPLibMain.exe: exe/CMakeFiles/PPPLibMain.dir/PPPLibMain.cc.obj
../bin/PPPLibMain.exe: exe/CMakeFiles/PPPLibMain.dir/build.make
../bin/PPPLibMain.exe: ../bin/libPPPLib.a
../bin/PPPLibMain.exe: ../bin/libSolverLib.a
../bin/PPPLibMain.exe: ../bin/libReadLib.a
../bin/PPPLibMain.exe: ../bin/libDecodeRawLib.a
../bin/PPPLibMain.exe: ../bin/libGnssLib.a
../bin/PPPLibMain.exe: ../bin/libInsLib.a
../bin/PPPLibMain.exe: ../bin/libAdjLib.a
../bin/PPPLibMain.exe: ../bin/libOutLib.a
../bin/PPPLibMain.exe: ../bin/libCmnLib.a
../bin/PPPLibMain.exe: ../bin/libLogLib.a
../bin/PPPLibMain.exe: exe/CMakeFiles/PPPLibMain.dir/linklibs.rsp
../bin/PPPLibMain.exe: exe/CMakeFiles/PPPLibMain.dir/objects1.rsp
../bin/PPPLibMain.exe: exe/CMakeFiles/PPPLibMain.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=I:\code\git\PPP\PPPLib\build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ..\..\bin\PPPLibMain.exe"
	cd /d I:\code\git\PPP\PPPLib\build\exe && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\PPPLibMain.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
exe/CMakeFiles/PPPLibMain.dir/build: ../bin/PPPLibMain.exe

.PHONY : exe/CMakeFiles/PPPLibMain.dir/build

exe/CMakeFiles/PPPLibMain.dir/clean:
	cd /d I:\code\git\PPP\PPPLib\build\exe && $(CMAKE_COMMAND) -P CMakeFiles\PPPLibMain.dir\cmake_clean.cmake
.PHONY : exe/CMakeFiles/PPPLibMain.dir/clean

exe/CMakeFiles/PPPLibMain.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" I:\code\git\PPP\PPPLib I:\code\git\PPP\PPPLib\exe I:\code\git\PPP\PPPLib\build I:\code\git\PPP\PPPLib\build\exe I:\code\git\PPP\PPPLib\build\exe\CMakeFiles\PPPLibMain.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : exe/CMakeFiles/PPPLibMain.dir/depend

