# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

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
CMAKE_COMMAND = /home/wzy/cmake-3.20.0-rc2-linux-x86_64/bin/cmake

# The command to remove a file.
RM = /home/wzy/cmake-3.20.0-rc2-linux-x86_64/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/wzy/Dumpy

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/wzy/Dumpy/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/Dumpy.dir/depend.make
# Include the progress variables for this target.
include CMakeFiles/Dumpy.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Dumpy.dir/flags.make

CMakeFiles/Dumpy.dir/src/main.cpp.o: CMakeFiles/Dumpy.dir/flags.make
CMakeFiles/Dumpy.dir/src/main.cpp.o: ../src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/wzy/Dumpy/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Dumpy.dir/src/main.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Dumpy.dir/src/main.cpp.o -c /home/wzy/Dumpy/src/main.cpp

CMakeFiles/Dumpy.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Dumpy.dir/src/main.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/wzy/Dumpy/src/main.cpp > CMakeFiles/Dumpy.dir/src/main.cpp.i

CMakeFiles/Dumpy.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Dumpy.dir/src/main.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/wzy/Dumpy/src/main.cpp -o CMakeFiles/Dumpy.dir/src/main.cpp.s

CMakeFiles/Dumpy.dir/src/Utils/FileUtil.cpp.o: CMakeFiles/Dumpy.dir/flags.make
CMakeFiles/Dumpy.dir/src/Utils/FileUtil.cpp.o: ../src/Utils/FileUtil.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/wzy/Dumpy/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/Dumpy.dir/src/Utils/FileUtil.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Dumpy.dir/src/Utils/FileUtil.cpp.o -c /home/wzy/Dumpy/src/Utils/FileUtil.cpp

CMakeFiles/Dumpy.dir/src/Utils/FileUtil.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Dumpy.dir/src/Utils/FileUtil.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/wzy/Dumpy/src/Utils/FileUtil.cpp > CMakeFiles/Dumpy.dir/src/Utils/FileUtil.cpp.i

CMakeFiles/Dumpy.dir/src/Utils/FileUtil.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Dumpy.dir/src/Utils/FileUtil.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/wzy/Dumpy/src/Utils/FileUtil.cpp -o CMakeFiles/Dumpy.dir/src/Utils/FileUtil.cpp.s

CMakeFiles/Dumpy.dir/src/Utils/TimeSeriesUtil.cpp.o: CMakeFiles/Dumpy.dir/flags.make
CMakeFiles/Dumpy.dir/src/Utils/TimeSeriesUtil.cpp.o: ../src/Utils/TimeSeriesUtil.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/wzy/Dumpy/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/Dumpy.dir/src/Utils/TimeSeriesUtil.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Dumpy.dir/src/Utils/TimeSeriesUtil.cpp.o -c /home/wzy/Dumpy/src/Utils/TimeSeriesUtil.cpp

CMakeFiles/Dumpy.dir/src/Utils/TimeSeriesUtil.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Dumpy.dir/src/Utils/TimeSeriesUtil.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/wzy/Dumpy/src/Utils/TimeSeriesUtil.cpp > CMakeFiles/Dumpy.dir/src/Utils/TimeSeriesUtil.cpp.i

CMakeFiles/Dumpy.dir/src/Utils/TimeSeriesUtil.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Dumpy.dir/src/Utils/TimeSeriesUtil.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/wzy/Dumpy/src/Utils/TimeSeriesUtil.cpp -o CMakeFiles/Dumpy.dir/src/Utils/TimeSeriesUtil.cpp.s

CMakeFiles/Dumpy.dir/src/Utils/SaxUtil.cpp.o: CMakeFiles/Dumpy.dir/flags.make
CMakeFiles/Dumpy.dir/src/Utils/SaxUtil.cpp.o: ../src/Utils/SaxUtil.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/wzy/Dumpy/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/Dumpy.dir/src/Utils/SaxUtil.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Dumpy.dir/src/Utils/SaxUtil.cpp.o -c /home/wzy/Dumpy/src/Utils/SaxUtil.cpp

CMakeFiles/Dumpy.dir/src/Utils/SaxUtil.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Dumpy.dir/src/Utils/SaxUtil.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/wzy/Dumpy/src/Utils/SaxUtil.cpp > CMakeFiles/Dumpy.dir/src/Utils/SaxUtil.cpp.i

CMakeFiles/Dumpy.dir/src/Utils/SaxUtil.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Dumpy.dir/src/Utils/SaxUtil.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/wzy/Dumpy/src/Utils/SaxUtil.cpp -o CMakeFiles/Dumpy.dir/src/Utils/SaxUtil.cpp.s

CMakeFiles/Dumpy.dir/src/Utils/MathUtil.cpp.o: CMakeFiles/Dumpy.dir/flags.make
CMakeFiles/Dumpy.dir/src/Utils/MathUtil.cpp.o: ../src/Utils/MathUtil.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/wzy/Dumpy/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/Dumpy.dir/src/Utils/MathUtil.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Dumpy.dir/src/Utils/MathUtil.cpp.o -c /home/wzy/Dumpy/src/Utils/MathUtil.cpp

CMakeFiles/Dumpy.dir/src/Utils/MathUtil.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Dumpy.dir/src/Utils/MathUtil.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/wzy/Dumpy/src/Utils/MathUtil.cpp > CMakeFiles/Dumpy.dir/src/Utils/MathUtil.cpp.i

CMakeFiles/Dumpy.dir/src/Utils/MathUtil.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Dumpy.dir/src/Utils/MathUtil.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/wzy/Dumpy/src/Utils/MathUtil.cpp -o CMakeFiles/Dumpy.dir/src/Utils/MathUtil.cpp.s

CMakeFiles/Dumpy.dir/src/IndexConstruction/GraphConstruction.cpp.o: CMakeFiles/Dumpy.dir/flags.make
CMakeFiles/Dumpy.dir/src/IndexConstruction/GraphConstruction.cpp.o: ../src/IndexConstruction/GraphConstruction.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/wzy/Dumpy/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/Dumpy.dir/src/IndexConstruction/GraphConstruction.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Dumpy.dir/src/IndexConstruction/GraphConstruction.cpp.o -c /home/wzy/Dumpy/src/IndexConstruction/GraphConstruction.cpp

CMakeFiles/Dumpy.dir/src/IndexConstruction/GraphConstruction.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Dumpy.dir/src/IndexConstruction/GraphConstruction.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/wzy/Dumpy/src/IndexConstruction/GraphConstruction.cpp > CMakeFiles/Dumpy.dir/src/IndexConstruction/GraphConstruction.cpp.i

CMakeFiles/Dumpy.dir/src/IndexConstruction/GraphConstruction.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Dumpy.dir/src/IndexConstruction/GraphConstruction.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/wzy/Dumpy/src/IndexConstruction/GraphConstruction.cpp -o CMakeFiles/Dumpy.dir/src/IndexConstruction/GraphConstruction.cpp.s

CMakeFiles/Dumpy.dir/src/PqItemSeries.cpp.o: CMakeFiles/Dumpy.dir/flags.make
CMakeFiles/Dumpy.dir/src/PqItemSeries.cpp.o: ../src/PqItemSeries.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/wzy/Dumpy/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/Dumpy.dir/src/PqItemSeries.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Dumpy.dir/src/PqItemSeries.cpp.o -c /home/wzy/Dumpy/src/PqItemSeries.cpp

CMakeFiles/Dumpy.dir/src/PqItemSeries.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Dumpy.dir/src/PqItemSeries.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/wzy/Dumpy/src/PqItemSeries.cpp > CMakeFiles/Dumpy.dir/src/PqItemSeries.cpp.i

CMakeFiles/Dumpy.dir/src/PqItemSeries.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Dumpy.dir/src/PqItemSeries.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/wzy/Dumpy/src/PqItemSeries.cpp -o CMakeFiles/Dumpy.dir/src/PqItemSeries.cpp.s

CMakeFiles/Dumpy.dir/src/Utils/INIReader.cpp.o: CMakeFiles/Dumpy.dir/flags.make
CMakeFiles/Dumpy.dir/src/Utils/INIReader.cpp.o: ../src/Utils/INIReader.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/wzy/Dumpy/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/Dumpy.dir/src/Utils/INIReader.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Dumpy.dir/src/Utils/INIReader.cpp.o -c /home/wzy/Dumpy/src/Utils/INIReader.cpp

CMakeFiles/Dumpy.dir/src/Utils/INIReader.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Dumpy.dir/src/Utils/INIReader.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/wzy/Dumpy/src/Utils/INIReader.cpp > CMakeFiles/Dumpy.dir/src/Utils/INIReader.cpp.i

CMakeFiles/Dumpy.dir/src/Utils/INIReader.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Dumpy.dir/src/Utils/INIReader.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/wzy/Dumpy/src/Utils/INIReader.cpp -o CMakeFiles/Dumpy.dir/src/Utils/INIReader.cpp.s

CMakeFiles/Dumpy.dir/src/Utils/ini.c.o: CMakeFiles/Dumpy.dir/flags.make
CMakeFiles/Dumpy.dir/src/Utils/ini.c.o: ../src/Utils/ini.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/wzy/Dumpy/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building C object CMakeFiles/Dumpy.dir/src/Utils/ini.c.o"
	/usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/Dumpy.dir/src/Utils/ini.c.o -c /home/wzy/Dumpy/src/Utils/ini.c

CMakeFiles/Dumpy.dir/src/Utils/ini.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/Dumpy.dir/src/Utils/ini.c.i"
	/usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/wzy/Dumpy/src/Utils/ini.c > CMakeFiles/Dumpy.dir/src/Utils/ini.c.i

CMakeFiles/Dumpy.dir/src/Utils/ini.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/Dumpy.dir/src/Utils/ini.c.s"
	/usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/wzy/Dumpy/src/Utils/ini.c -o CMakeFiles/Dumpy.dir/src/Utils/ini.c.s

CMakeFiles/Dumpy.dir/src/Const.cpp.o: CMakeFiles/Dumpy.dir/flags.make
CMakeFiles/Dumpy.dir/src/Const.cpp.o: ../src/Const.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/wzy/Dumpy/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/Dumpy.dir/src/Const.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Dumpy.dir/src/Const.cpp.o -c /home/wzy/Dumpy/src/Const.cpp

CMakeFiles/Dumpy.dir/src/Const.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Dumpy.dir/src/Const.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/wzy/Dumpy/src/Const.cpp > CMakeFiles/Dumpy.dir/src/Const.cpp.i

CMakeFiles/Dumpy.dir/src/Const.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Dumpy.dir/src/Const.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/wzy/Dumpy/src/Const.cpp -o CMakeFiles/Dumpy.dir/src/Const.cpp.s

CMakeFiles/Dumpy.dir/src/IndexConstruction/DumpyNode.cpp.o: CMakeFiles/Dumpy.dir/flags.make
CMakeFiles/Dumpy.dir/src/IndexConstruction/DumpyNode.cpp.o: ../src/IndexConstruction/DumpyNode.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/wzy/Dumpy/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object CMakeFiles/Dumpy.dir/src/IndexConstruction/DumpyNode.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Dumpy.dir/src/IndexConstruction/DumpyNode.cpp.o -c /home/wzy/Dumpy/src/IndexConstruction/DumpyNode.cpp

CMakeFiles/Dumpy.dir/src/IndexConstruction/DumpyNode.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Dumpy.dir/src/IndexConstruction/DumpyNode.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/wzy/Dumpy/src/IndexConstruction/DumpyNode.cpp > CMakeFiles/Dumpy.dir/src/IndexConstruction/DumpyNode.cpp.i

CMakeFiles/Dumpy.dir/src/IndexConstruction/DumpyNode.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Dumpy.dir/src/IndexConstruction/DumpyNode.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/wzy/Dumpy/src/IndexConstruction/DumpyNode.cpp -o CMakeFiles/Dumpy.dir/src/IndexConstruction/DumpyNode.cpp.s

CMakeFiles/Dumpy.dir/src/SearchEngine/DumpySearcher.cpp.o: CMakeFiles/Dumpy.dir/flags.make
CMakeFiles/Dumpy.dir/src/SearchEngine/DumpySearcher.cpp.o: ../src/SearchEngine/DumpySearcher.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/wzy/Dumpy/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object CMakeFiles/Dumpy.dir/src/SearchEngine/DumpySearcher.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Dumpy.dir/src/SearchEngine/DumpySearcher.cpp.o -c /home/wzy/Dumpy/src/SearchEngine/DumpySearcher.cpp

CMakeFiles/Dumpy.dir/src/SearchEngine/DumpySearcher.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Dumpy.dir/src/SearchEngine/DumpySearcher.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/wzy/Dumpy/src/SearchEngine/DumpySearcher.cpp > CMakeFiles/Dumpy.dir/src/SearchEngine/DumpySearcher.cpp.i

CMakeFiles/Dumpy.dir/src/SearchEngine/DumpySearcher.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Dumpy.dir/src/SearchEngine/DumpySearcher.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/wzy/Dumpy/src/SearchEngine/DumpySearcher.cpp -o CMakeFiles/Dumpy.dir/src/SearchEngine/DumpySearcher.cpp.s

CMakeFiles/Dumpy.dir/src/IndexConstruction/DumpyFuzzy.cpp.o: CMakeFiles/Dumpy.dir/flags.make
CMakeFiles/Dumpy.dir/src/IndexConstruction/DumpyFuzzy.cpp.o: ../src/IndexConstruction/DumpyFuzzy.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/wzy/Dumpy/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building CXX object CMakeFiles/Dumpy.dir/src/IndexConstruction/DumpyFuzzy.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Dumpy.dir/src/IndexConstruction/DumpyFuzzy.cpp.o -c /home/wzy/Dumpy/src/IndexConstruction/DumpyFuzzy.cpp

CMakeFiles/Dumpy.dir/src/IndexConstruction/DumpyFuzzy.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Dumpy.dir/src/IndexConstruction/DumpyFuzzy.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/wzy/Dumpy/src/IndexConstruction/DumpyFuzzy.cpp > CMakeFiles/Dumpy.dir/src/IndexConstruction/DumpyFuzzy.cpp.i

CMakeFiles/Dumpy.dir/src/IndexConstruction/DumpyFuzzy.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Dumpy.dir/src/IndexConstruction/DumpyFuzzy.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/wzy/Dumpy/src/IndexConstruction/DumpyFuzzy.cpp -o CMakeFiles/Dumpy.dir/src/IndexConstruction/DumpyFuzzy.cpp.s

# Object files for target Dumpy
Dumpy_OBJECTS = \
"CMakeFiles/Dumpy.dir/src/main.cpp.o" \
"CMakeFiles/Dumpy.dir/src/Utils/FileUtil.cpp.o" \
"CMakeFiles/Dumpy.dir/src/Utils/TimeSeriesUtil.cpp.o" \
"CMakeFiles/Dumpy.dir/src/Utils/SaxUtil.cpp.o" \
"CMakeFiles/Dumpy.dir/src/Utils/MathUtil.cpp.o" \
"CMakeFiles/Dumpy.dir/src/IndexConstruction/GraphConstruction.cpp.o" \
"CMakeFiles/Dumpy.dir/src/PqItemSeries.cpp.o" \
"CMakeFiles/Dumpy.dir/src/Utils/INIReader.cpp.o" \
"CMakeFiles/Dumpy.dir/src/Utils/ini.c.o" \
"CMakeFiles/Dumpy.dir/src/Const.cpp.o" \
"CMakeFiles/Dumpy.dir/src/IndexConstruction/DumpyNode.cpp.o" \
"CMakeFiles/Dumpy.dir/src/SearchEngine/DumpySearcher.cpp.o" \
"CMakeFiles/Dumpy.dir/src/IndexConstruction/DumpyFuzzy.cpp.o"

# External object files for target Dumpy
Dumpy_EXTERNAL_OBJECTS =

Dumpy: CMakeFiles/Dumpy.dir/src/main.cpp.o
Dumpy: CMakeFiles/Dumpy.dir/src/Utils/FileUtil.cpp.o
Dumpy: CMakeFiles/Dumpy.dir/src/Utils/TimeSeriesUtil.cpp.o
Dumpy: CMakeFiles/Dumpy.dir/src/Utils/SaxUtil.cpp.o
Dumpy: CMakeFiles/Dumpy.dir/src/Utils/MathUtil.cpp.o
Dumpy: CMakeFiles/Dumpy.dir/src/IndexConstruction/GraphConstruction.cpp.o
Dumpy: CMakeFiles/Dumpy.dir/src/PqItemSeries.cpp.o
Dumpy: CMakeFiles/Dumpy.dir/src/Utils/INIReader.cpp.o
Dumpy: CMakeFiles/Dumpy.dir/src/Utils/ini.c.o
Dumpy: CMakeFiles/Dumpy.dir/src/Const.cpp.o
Dumpy: CMakeFiles/Dumpy.dir/src/IndexConstruction/DumpyNode.cpp.o
Dumpy: CMakeFiles/Dumpy.dir/src/SearchEngine/DumpySearcher.cpp.o
Dumpy: CMakeFiles/Dumpy.dir/src/IndexConstruction/DumpyFuzzy.cpp.o
Dumpy: CMakeFiles/Dumpy.dir/build.make
Dumpy: CMakeFiles/Dumpy.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/wzy/Dumpy/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Linking CXX executable Dumpy"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Dumpy.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Dumpy.dir/build: Dumpy
.PHONY : CMakeFiles/Dumpy.dir/build

CMakeFiles/Dumpy.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Dumpy.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Dumpy.dir/clean

CMakeFiles/Dumpy.dir/depend:
	cd /home/wzy/Dumpy/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/wzy/Dumpy /home/wzy/Dumpy /home/wzy/Dumpy/cmake-build-debug /home/wzy/Dumpy/cmake-build-debug /home/wzy/Dumpy/cmake-build-debug/CMakeFiles/Dumpy.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Dumpy.dir/depend

