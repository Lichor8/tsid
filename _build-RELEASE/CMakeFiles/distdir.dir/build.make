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

# Utility rule file for distdir.

# Include the progress variables for this target.
include CMakeFiles/distdir.dir/progress.make

CMakeFiles/distdir:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/opt/openrobots/src/tsid-fork/_build-RELEASE/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating dist directory..."
	cd /opt/openrobots/src/tsid-fork && rm -f /tmp/tsid.tar && /opt/openrobots/src/tsid-fork/cmake/git-archive-all.sh --prefix tsid-1.1.0-7-gb18e-dirty/ tsid.tar && cd /opt/openrobots/src/tsid-fork/_build-RELEASE/ && ( test -d tsid-1.1.0-7-gb18e-dirty && find tsid-1.1.0-7-gb18e-dirty/ -type d -print0 | xargs -0 chmod a+w || true ) && rm -rf tsid-1.1.0-7-gb18e-dirty/ && /bin/tar xf /opt/openrobots/src/tsid-fork/tsid.tar && echo 1.1.0-7-gb18e-dirty > /opt/openrobots/src/tsid-fork/_build-RELEASE/tsid-1.1.0-7-gb18e-dirty/.version && /opt/openrobots/src/tsid-fork/cmake/gitlog-to-changelog > /opt/openrobots/src/tsid-fork/_build-RELEASE/tsid-1.1.0-7-gb18e-dirty/ChangeLog && rm -f /opt/openrobots/src/tsid-fork/tsid.tar

distdir: CMakeFiles/distdir
distdir: CMakeFiles/distdir.dir/build.make

.PHONY : distdir

# Rule to build all files generated by this target.
CMakeFiles/distdir.dir/build: distdir

.PHONY : CMakeFiles/distdir.dir/build

CMakeFiles/distdir.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/distdir.dir/cmake_clean.cmake
.PHONY : CMakeFiles/distdir.dir/clean

CMakeFiles/distdir.dir/depend:
	cd /opt/openrobots/src/tsid-fork/_build-RELEASE && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /opt/openrobots/src/tsid-fork /opt/openrobots/src/tsid-fork /opt/openrobots/src/tsid-fork/_build-RELEASE /opt/openrobots/src/tsid-fork/_build-RELEASE /opt/openrobots/src/tsid-fork/_build-RELEASE/CMakeFiles/distdir.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/distdir.dir/depend

