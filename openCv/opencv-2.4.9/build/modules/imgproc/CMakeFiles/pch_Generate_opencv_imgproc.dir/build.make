# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/cbassuser/avnOpticalPointing/opencv-2.4.9

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/cbassuser/avnOpticalPointing/opencv-2.4.9/build

# Utility rule file for pch_Generate_opencv_imgproc.

modules/imgproc/CMakeFiles/pch_Generate_opencv_imgproc: modules/imgproc/precomp.hpp.gch/opencv_imgproc_RELEASE.gch

modules/imgproc/precomp.hpp.gch/opencv_imgproc_RELEASE.gch: ../modules/imgproc/src/precomp.hpp
modules/imgproc/precomp.hpp.gch/opencv_imgproc_RELEASE.gch: modules/imgproc/precomp.hpp
modules/imgproc/precomp.hpp.gch/opencv_imgproc_RELEASE.gch: lib/libopencv_imgproc_pch_dephelp.a
	$(CMAKE_COMMAND) -E cmake_progress_report /home/cbassuser/avnOpticalPointing/opencv-2.4.9/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Generating precomp.hpp.gch/opencv_imgproc_RELEASE.gch"
	cd /home/cbassuser/avnOpticalPointing/opencv-2.4.9/build/modules/imgproc && /usr/bin/cmake -E make_directory /home/cbassuser/avnOpticalPointing/opencv-2.4.9/build/modules/imgproc/precomp.hpp.gch
	cd /home/cbassuser/avnOpticalPointing/opencv-2.4.9/build/modules/imgproc && /usr/bin/c++ -O3 -DNDEBUG -DNDEBUG -fPIC -I"/home/cbassuser/avnOpticalPointing/opencv-2.4.9/modules/core/include" -isystem"/home/cbassuser/avnOpticalPointing/opencv-2.4.9/build/modules/imgproc" -I"/home/cbassuser/avnOpticalPointing/opencv-2.4.9/modules/imgproc/src" -I"/home/cbassuser/avnOpticalPointing/opencv-2.4.9/modules/imgproc/include" -isystem"/home/cbassuser/avnOpticalPointing/opencv-2.4.9/build" -fsigned-char -W -Wall -Werror=return-type -Werror=address -Werror=sequence-point -Wformat -Werror=format-security -Wmissing-declarations -Wundef -Winit-self -Wpointer-arith -Wshadow -Wsign-promo -fdiagnostics-show-option -Wno-long-long -pthread -fomit-frame-pointer -msse -msse2 -msse3 -ffunction-sections -DCVAPI_EXPORTS -x c++-header -o /home/cbassuser/avnOpticalPointing/opencv-2.4.9/build/modules/imgproc/precomp.hpp.gch/opencv_imgproc_RELEASE.gch /home/cbassuser/avnOpticalPointing/opencv-2.4.9/build/modules/imgproc/precomp.hpp

modules/imgproc/precomp.hpp: ../modules/imgproc/src/precomp.hpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/cbassuser/avnOpticalPointing/opencv-2.4.9/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Generating precomp.hpp"
	cd /home/cbassuser/avnOpticalPointing/opencv-2.4.9/build/modules/imgproc && /usr/bin/cmake -E copy /home/cbassuser/avnOpticalPointing/opencv-2.4.9/modules/imgproc/src/precomp.hpp /home/cbassuser/avnOpticalPointing/opencv-2.4.9/build/modules/imgproc/precomp.hpp

pch_Generate_opencv_imgproc: modules/imgproc/CMakeFiles/pch_Generate_opencv_imgproc
pch_Generate_opencv_imgproc: modules/imgproc/precomp.hpp.gch/opencv_imgproc_RELEASE.gch
pch_Generate_opencv_imgproc: modules/imgproc/precomp.hpp
pch_Generate_opencv_imgproc: modules/imgproc/CMakeFiles/pch_Generate_opencv_imgproc.dir/build.make
.PHONY : pch_Generate_opencv_imgproc

# Rule to build all files generated by this target.
modules/imgproc/CMakeFiles/pch_Generate_opencv_imgproc.dir/build: pch_Generate_opencv_imgproc
.PHONY : modules/imgproc/CMakeFiles/pch_Generate_opencv_imgproc.dir/build

modules/imgproc/CMakeFiles/pch_Generate_opencv_imgproc.dir/clean:
	cd /home/cbassuser/avnOpticalPointing/opencv-2.4.9/build/modules/imgproc && $(CMAKE_COMMAND) -P CMakeFiles/pch_Generate_opencv_imgproc.dir/cmake_clean.cmake
.PHONY : modules/imgproc/CMakeFiles/pch_Generate_opencv_imgproc.dir/clean

modules/imgproc/CMakeFiles/pch_Generate_opencv_imgproc.dir/depend:
	cd /home/cbassuser/avnOpticalPointing/opencv-2.4.9/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/cbassuser/avnOpticalPointing/opencv-2.4.9 /home/cbassuser/avnOpticalPointing/opencv-2.4.9/modules/imgproc /home/cbassuser/avnOpticalPointing/opencv-2.4.9/build /home/cbassuser/avnOpticalPointing/opencv-2.4.9/build/modules/imgproc /home/cbassuser/avnOpticalPointing/opencv-2.4.9/build/modules/imgproc/CMakeFiles/pch_Generate_opencv_imgproc.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : modules/imgproc/CMakeFiles/pch_Generate_opencv_imgproc.dir/depend

