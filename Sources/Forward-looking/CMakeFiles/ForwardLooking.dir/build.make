# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

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
CMAKE_SOURCE_DIR = "/home/oberlan/Google Drive UFES/bromro2/Sources/Forward-looking"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/oberlan/Google Drive UFES/bromro2/Sources/Forward-looking"

# Include any dependencies generated for this target.
include CMakeFiles/ForwardLooking.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/ForwardLooking.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/ForwardLooking.dir/flags.make

CMakeFiles/ForwardLooking.dir/main.cpp.o: CMakeFiles/ForwardLooking.dir/flags.make
CMakeFiles/ForwardLooking.dir/main.cpp.o: main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/oberlan/Google Drive UFES/bromro2/Sources/Forward-looking/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/ForwardLooking.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ForwardLooking.dir/main.cpp.o -c "/home/oberlan/Google Drive UFES/bromro2/Sources/Forward-looking/main.cpp"

CMakeFiles/ForwardLooking.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ForwardLooking.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/oberlan/Google Drive UFES/bromro2/Sources/Forward-looking/main.cpp" > CMakeFiles/ForwardLooking.dir/main.cpp.i

CMakeFiles/ForwardLooking.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ForwardLooking.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/oberlan/Google Drive UFES/bromro2/Sources/Forward-looking/main.cpp" -o CMakeFiles/ForwardLooking.dir/main.cpp.s

CMakeFiles/ForwardLooking.dir/main.cpp.o.requires:

.PHONY : CMakeFiles/ForwardLooking.dir/main.cpp.o.requires

CMakeFiles/ForwardLooking.dir/main.cpp.o.provides: CMakeFiles/ForwardLooking.dir/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/ForwardLooking.dir/build.make CMakeFiles/ForwardLooking.dir/main.cpp.o.provides.build
.PHONY : CMakeFiles/ForwardLooking.dir/main.cpp.o.provides

CMakeFiles/ForwardLooking.dir/main.cpp.o.provides.build: CMakeFiles/ForwardLooking.dir/main.cpp.o


CMakeFiles/ForwardLooking.dir/scpp_assert.cpp.o: CMakeFiles/ForwardLooking.dir/flags.make
CMakeFiles/ForwardLooking.dir/scpp_assert.cpp.o: scpp_assert.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/oberlan/Google Drive UFES/bromro2/Sources/Forward-looking/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/ForwardLooking.dir/scpp_assert.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ForwardLooking.dir/scpp_assert.cpp.o -c "/home/oberlan/Google Drive UFES/bromro2/Sources/Forward-looking/scpp_assert.cpp"

CMakeFiles/ForwardLooking.dir/scpp_assert.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ForwardLooking.dir/scpp_assert.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/oberlan/Google Drive UFES/bromro2/Sources/Forward-looking/scpp_assert.cpp" > CMakeFiles/ForwardLooking.dir/scpp_assert.cpp.i

CMakeFiles/ForwardLooking.dir/scpp_assert.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ForwardLooking.dir/scpp_assert.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/oberlan/Google Drive UFES/bromro2/Sources/Forward-looking/scpp_assert.cpp" -o CMakeFiles/ForwardLooking.dir/scpp_assert.cpp.s

CMakeFiles/ForwardLooking.dir/scpp_assert.cpp.o.requires:

.PHONY : CMakeFiles/ForwardLooking.dir/scpp_assert.cpp.o.requires

CMakeFiles/ForwardLooking.dir/scpp_assert.cpp.o.provides: CMakeFiles/ForwardLooking.dir/scpp_assert.cpp.o.requires
	$(MAKE) -f CMakeFiles/ForwardLooking.dir/build.make CMakeFiles/ForwardLooking.dir/scpp_assert.cpp.o.provides.build
.PHONY : CMakeFiles/ForwardLooking.dir/scpp_assert.cpp.o.provides

CMakeFiles/ForwardLooking.dir/scpp_assert.cpp.o.provides.build: CMakeFiles/ForwardLooking.dir/scpp_assert.cpp.o


# Object files for target ForwardLooking
ForwardLooking_OBJECTS = \
"CMakeFiles/ForwardLooking.dir/main.cpp.o" \
"CMakeFiles/ForwardLooking.dir/scpp_assert.cpp.o"

# External object files for target ForwardLooking
ForwardLooking_EXTERNAL_OBJECTS =

ForwardLooking: CMakeFiles/ForwardLooking.dir/main.cpp.o
ForwardLooking: CMakeFiles/ForwardLooking.dir/scpp_assert.cpp.o
ForwardLooking: CMakeFiles/ForwardLooking.dir/build.make
ForwardLooking: /opt/ibm/ILOG/CPLEX_Studio1210/concert/lib/x86-64_linux/static_pic/libconcert.a
ForwardLooking: /opt/ibm/ILOG/CPLEX_Studio1210/cplex/lib/x86-64_linux/static_pic/libilocplex.a
ForwardLooking: /opt/ibm/ILOG/CPLEX_Studio1210/cplex/lib/x86-64_linux/static_pic/libcplex.a
ForwardLooking: /opt/ibm/ILOG/CPLEX_Studio1210/concert/lib/x86-64_linux/static_pic/libconcert.a
ForwardLooking: /opt/ibm/ILOG/CPLEX_Studio1210/cplex/lib/x86-64_linux/static_pic/libilocplex.a
ForwardLooking: /opt/ibm/ILOG/CPLEX_Studio1210/cplex/lib/x86-64_linux/static_pic/libcplex.a
ForwardLooking: CMakeFiles/ForwardLooking.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/home/oberlan/Google Drive UFES/bromro2/Sources/Forward-looking/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable ForwardLooking"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ForwardLooking.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/ForwardLooking.dir/build: ForwardLooking

.PHONY : CMakeFiles/ForwardLooking.dir/build

CMakeFiles/ForwardLooking.dir/requires: CMakeFiles/ForwardLooking.dir/main.cpp.o.requires
CMakeFiles/ForwardLooking.dir/requires: CMakeFiles/ForwardLooking.dir/scpp_assert.cpp.o.requires

.PHONY : CMakeFiles/ForwardLooking.dir/requires

CMakeFiles/ForwardLooking.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/ForwardLooking.dir/cmake_clean.cmake
.PHONY : CMakeFiles/ForwardLooking.dir/clean

CMakeFiles/ForwardLooking.dir/depend:
	cd "/home/oberlan/Google Drive UFES/bromro2/Sources/Forward-looking" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/oberlan/Google Drive UFES/bromro2/Sources/Forward-looking" "/home/oberlan/Google Drive UFES/bromro2/Sources/Forward-looking" "/home/oberlan/Google Drive UFES/bromro2/Sources/Forward-looking" "/home/oberlan/Google Drive UFES/bromro2/Sources/Forward-looking" "/home/oberlan/Google Drive UFES/bromro2/Sources/Forward-looking/CMakeFiles/ForwardLooking.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/ForwardLooking.dir/depend

