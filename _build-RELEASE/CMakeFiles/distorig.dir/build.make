# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Produce verbose output by default.
VERBOSE = 1

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
CMAKE_SOURCE_DIR = /opt/openrobots/src/tsid-fork

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /opt/openrobots/src/tsid-fork/_build-RELEASE

# Utility rule file for distorig.

# Include the progress variables for this target.
include CMakeFiles/distorig.dir/progress.make

CMakeFiles/distorig:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/opt/openrobots/src/tsid-fork/_build-RELEASE/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating orig tarball..."
	cmake -E copy tsid-1.1.0-2-ga574-dirty.tar.gz tsid-1.1.0-2-ga574-dirty.orig.tar.gz

distorig: CMakeFiles/distorig
distorig: CMakeFiles/distorig.dir/build.make

.PHONY : distorig

# Rule to build all files generated by this target.
CMakeFiles/distorig.dir/build: distorig

.PHONY : CMakeFiles/distorig.dir/build

CMakeFiles/distorig.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/distorig.dir/cmake_clean.cmake
.PHONY : CMakeFiles/distorig.dir/clean

CMakeFiles/distorig.dir/depend:
	cd /opt/openrobots/src/tsid-fork/_build-RELEASE && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /opt/openrobots/src/tsid-fork /opt/openrobots/src/tsid-fork /opt/openrobots/src/tsid-fork/_build-RELEASE /opt/openrobots/src/tsid-fork/_build-RELEASE /opt/openrobots/src/tsid-fork/_build-RELEASE/CMakeFiles/distorig.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/distorig.dir/depend

