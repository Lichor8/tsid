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

# Utility rule file for dist_tarbz2.

# Include the progress variables for this target.
include CMakeFiles/dist_tarbz2.dir/progress.make

CMakeFiles/dist_tarbz2:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/opt/openrobots/src/tsid-fork/_build-RELEASE/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating tar.bz2 tarball and its signature..."
	/bin/tar -cjf tsid-1.1.0-9-ge81c-dirty.tar.bz2 tsid-1.1.0-9-ge81c-dirty/ && /usr/bin/gpg --detach-sign --armor -o /opt/openrobots/src/tsid-fork/_build-RELEASE/tsid-1.1.0-9-ge81c-dirty.tar.bz2.sig /opt/openrobots/src/tsid-fork/_build-RELEASE/tsid-1.1.0-9-ge81c-dirty.tar.bz2

dist_tarbz2: CMakeFiles/dist_tarbz2
dist_tarbz2: CMakeFiles/dist_tarbz2.dir/build.make

.PHONY : dist_tarbz2

# Rule to build all files generated by this target.
CMakeFiles/dist_tarbz2.dir/build: dist_tarbz2

.PHONY : CMakeFiles/dist_tarbz2.dir/build

CMakeFiles/dist_tarbz2.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/dist_tarbz2.dir/cmake_clean.cmake
.PHONY : CMakeFiles/dist_tarbz2.dir/clean

CMakeFiles/dist_tarbz2.dir/depend:
	cd /opt/openrobots/src/tsid-fork/_build-RELEASE && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /opt/openrobots/src/tsid-fork /opt/openrobots/src/tsid-fork /opt/openrobots/src/tsid-fork/_build-RELEASE /opt/openrobots/src/tsid-fork/_build-RELEASE /opt/openrobots/src/tsid-fork/_build-RELEASE/CMakeFiles/dist_tarbz2.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/dist_tarbz2.dir/depend

