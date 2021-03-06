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

# Utility rule file for pch_Generate_opencv_test_gpu.

modules/gpu/CMakeFiles/pch_Generate_opencv_test_gpu: modules/gpu/test_precomp.hpp.gch/opencv_test_gpu_RELEASE.gch

modules/gpu/test_precomp.hpp.gch/opencv_test_gpu_RELEASE.gch: ../modules/gpu/test/test_precomp.hpp
modules/gpu/test_precomp.hpp.gch/opencv_test_gpu_RELEASE.gch: modules/gpu/test_precomp.hpp
modules/gpu/test_precomp.hpp.gch/opencv_test_gpu_RELEASE.gch: lib/libopencv_test_gpu_pch_dephelp.a
	$(CMAKE_COMMAND) -E cmake_progress_report /home/cbassuser/avnOpticalPointing/opencv-2.4.9/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Generating test_precomp.hpp.gch/opencv_test_gpu_RELEASE.gch"
	cd /home/cbassuser/avnOpticalPointing/opencv-2.4.9/build/modules/gpu && /usr/bin/cmake -E make_directory /home/cbassuser/avnOpticalPointing/opencv-2.4.9/build/modules/gpu/test_precomp.hpp.gch
	cd /home/cbassuser/avnOpticalPointing/opencv-2.4.9/build/modules/gpu && /usr/bin/c++ -O3 -DNDEBUG -DNDEBUG -I"/home/cbassuser/avnOpticalPointing/opencv-2.4.9/modules/gpu/test" -I"/home/cbassuser/avnOpticalPointing/opencv-2.4.9/modules/features2d/include" -I"/home/cbassuser/avnOpticalPointing/opencv-2.4.9/modules/highgui/include" -I"/home/cbassuser/avnOpticalPointing/opencv-2.4.9/modules/imgproc/include" -I"/home/cbassuser/avnOpticalPointing/opencv-2.4.9/modules/flann/include" -I"/home/cbassuser/avnOpticalPointing/opencv-2.4.9/modules/core/include" -I"/home/cbassuser/avnOpticalPointing/opencv-2.4.9/modules/ts/include" -I"/home/cbassuser/avnOpticalPointing/opencv-2.4.9/modules/gpu/include" -I"/home/cbassuser/avnOpticalPointing/opencv-2.4.9/modules/photo/include" -I"/home/cbassuser/avnOpticalPointing/opencv-2.4.9/modules/objdetect/include" -I"/home/cbassuser/avnOpticalPointing/opencv-2.4.9/modules/legacy/include" -I"/home/cbassuser/avnOpticalPointing/opencv-2.4.9/modules/video/include" -I"/home/cbassuser/avnOpticalPointing/opencv-2.4.9/modules/ml/include" -I"/home/cbassuser/avnOpticalPointing/opencv-2.4.9/modules/calib3d/include" -isystem"/home/cbassuser/avnOpticalPointing/opencv-2.4.9/build/modules/gpu" -I"/home/cbassuser/avnOpticalPointing/opencv-2.4.9/modules/gpu/src" -I"/home/cbassuser/avnOpticalPointing/opencv-2.4.9/modules/gpu/src/cuda" -isystem"/home/cbassuser/avnOpticalPointing/opencv-2.4.9/build" -fsigned-char -W -Wall -Werror=return-type -Werror=address -Werror=sequence-point -Wformat -Werror=format-security -Wmissing-declarations -Wundef -Winit-self -Wpointer-arith -Wshadow -Wsign-promo -fdiagnostics-show-option -Wno-long-long -pthread -fomit-frame-pointer -msse -msse2 -msse3 -ffunction-sections -x c++-header -o /home/cbassuser/avnOpticalPointing/opencv-2.4.9/build/modules/gpu/test_precomp.hpp.gch/opencv_test_gpu_RELEASE.gch /home/cbassuser/avnOpticalPointing/opencv-2.4.9/build/modules/gpu/test_precomp.hpp

modules/gpu/test_precomp.hpp: ../modules/gpu/test/test_precomp.hpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/cbassuser/avnOpticalPointing/opencv-2.4.9/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Generating test_precomp.hpp"
	cd /home/cbassuser/avnOpticalPointing/opencv-2.4.9/build/modules/gpu && /usr/bin/cmake -E copy /home/cbassuser/avnOpticalPointing/opencv-2.4.9/modules/gpu/test/test_precomp.hpp /home/cbassuser/avnOpticalPointing/opencv-2.4.9/build/modules/gpu/test_precomp.hpp

pch_Generate_opencv_test_gpu: modules/gpu/CMakeFiles/pch_Generate_opencv_test_gpu
pch_Generate_opencv_test_gpu: modules/gpu/test_precomp.hpp.gch/opencv_test_gpu_RELEASE.gch
pch_Generate_opencv_test_gpu: modules/gpu/test_precomp.hpp
pch_Generate_opencv_test_gpu: modules/gpu/CMakeFiles/pch_Generate_opencv_test_gpu.dir/build.make
.PHONY : pch_Generate_opencv_test_gpu

# Rule to build all files generated by this target.
modules/gpu/CMakeFiles/pch_Generate_opencv_test_gpu.dir/build: pch_Generate_opencv_test_gpu
.PHONY : modules/gpu/CMakeFiles/pch_Generate_opencv_test_gpu.dir/build

modules/gpu/CMakeFiles/pch_Generate_opencv_test_gpu.dir/clean:
	cd /home/cbassuser/avnOpticalPointing/opencv-2.4.9/build/modules/gpu && $(CMAKE_COMMAND) -P CMakeFiles/pch_Generate_opencv_test_gpu.dir/cmake_clean.cmake
.PHONY : modules/gpu/CMakeFiles/pch_Generate_opencv_test_gpu.dir/clean

modules/gpu/CMakeFiles/pch_Generate_opencv_test_gpu.dir/depend:
	cd /home/cbassuser/avnOpticalPointing/opencv-2.4.9/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/cbassuser/avnOpticalPointing/opencv-2.4.9 /home/cbassuser/avnOpticalPointing/opencv-2.4.9/modules/gpu /home/cbassuser/avnOpticalPointing/opencv-2.4.9/build /home/cbassuser/avnOpticalPointing/opencv-2.4.9/build/modules/gpu /home/cbassuser/avnOpticalPointing/opencv-2.4.9/build/modules/gpu/CMakeFiles/pch_Generate_opencv_test_gpu.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : modules/gpu/CMakeFiles/pch_Generate_opencv_test_gpu.dir/depend

