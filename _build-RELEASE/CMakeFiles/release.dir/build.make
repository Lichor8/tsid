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

# Utility rule file for release.

# Include the progress variables for this target.
include CMakeFiles/release.dir/progress.make

CMakeFiles/release:
	! test x$$VERSION = x || ( echo Please\ set\ a\ version\ for\ this\ release && false ) && cd /opt/openrobots/src/tsid-fork && /usr/bin/git tag -s v$$VERSION -m Release\ of\ version\ $$VERSION. && cd /opt/openrobots/src/tsid-fork/_build-RELEASE && cmake /opt/openrobots/src/tsid-fork && make distcheck || ( echo Please\ fix\ distcheck\ first. && cd /opt/openrobots/src/tsid-fork && /usr/bin/git tag -d v$$VERSION && cd /opt/openrobots/src/tsid-fork/_build-RELEASE && cmake /opt/openrobots/src/tsid-fork && false ) && make dist && make distclean && echo Please,\ run\ 'git\ push\ --tags'\ and\ upload\ the\ tarball\ to\ github\ to\ finalize\ this\ release.

release: CMakeFiles/release
release: CMakeFiles/release.dir/build.make

.PHONY : release

# Rule to build all files generated by this target.
CMakeFiles/release.dir/build: release

.PHONY : CMakeFiles/release.dir/build

CMakeFiles/release.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/release.dir/cmake_clean.cmake
.PHONY : CMakeFiles/release.dir/clean

CMakeFiles/release.dir/depend:
	cd /opt/openrobots/src/tsid-fork/_build-RELEASE && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /opt/openrobots/src/tsid-fork /opt/openrobots/src/tsid-fork /opt/openrobots/src/tsid-fork/_build-RELEASE /opt/openrobots/src/tsid-fork/_build-RELEASE /opt/openrobots/src/tsid-fork/_build-RELEASE/CMakeFiles/release.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/release.dir/depend

