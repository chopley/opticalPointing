# Install script for directory: /home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/usr/local")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "Release")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "1")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "samples")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/OpenCV/samples/cpp" TYPE FILE PERMISSIONS OWNER_READ GROUP_READ WORLD_READ FILES
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/kalman.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/ffilldemo.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/imagelist_creator.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/watershed.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/contours2.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/drawing.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/intelperc_capture.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/freak_demo.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/distrans.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/detector_descriptor_matcher_evaluation.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/pca.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/descriptor_extractor_matcher.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/bagofwords_classification.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/detection_based_tracker_sample.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/squares.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/stereo_calib.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/bgfg_gmg.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/linemod.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/matching_to_many_images.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/hybridtrackingsample.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/rgbdodometry.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/houghcircles.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/convexhull.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/fabmap_sample.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/cout_mat.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/starter_imagelist.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/stereo_match.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/stitching_detailed.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/edge.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/retinaDemo.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/generic_descriptor_match.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/houghlines.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/facerec_demo.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/video_homography.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/video_dmtx.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/fback.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/OpenEXRimages_HDR_Retina_toneMapping.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/kmeans.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/inpaint.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/connected_components.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/fitellipse.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/openni_capture.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/bgfg_segm.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/calibration_artificial.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/3calibration.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/filestorage.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/detector_descriptor_evaluation.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/stitching.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/image.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/em.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/camshiftdemo.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/chamfer.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/logpolar_bsm.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/gencolors.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/image_sequence.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/meanshift_segmentation.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/build3dmodel.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/tvl1_optical_flow.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/peopledetect.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/morphology2.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/OpenEXRimages_HDR_Retina_toneMapping_video.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/points_classifier.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/dft.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/simpleflow_demo.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/matcher_simple.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/phase_corr.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/demhist.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/select3dobj.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/grabcut.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/opencv_version.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/videostab.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/minarea.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/delaunay2.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/calibration.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/brief_match_test.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/lkdemo.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/laplace.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/segment_objects.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/latentsvm_multidetect.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/starter_video.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/letter_recog.cpp"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/left05.jpg"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/board.jpg"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/left11.jpg"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/right14.jpg"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/left09.jpg"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/right06.jpg"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/left02.jpg"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/building.jpg"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/right07.jpg"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/left01.jpg"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/lena.jpg"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/left03.jpg"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/right08.jpg"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/right13.jpg"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/left12.jpg"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/left08.jpg"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/left14.jpg"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/right02.jpg"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/right01.jpg"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/fruits.jpg"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/baboon.jpg"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/right12.jpg"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/left13.jpg"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/right09.jpg"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/stuff.jpg"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/right11.jpg"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/left06.jpg"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/right05.jpg"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/left07.jpg"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/left04.jpg"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/right04.jpg"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/right03.jpg"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/pic5.png"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/tsukuba_l.png"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/pic2.png"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/pic4.png"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/pic6.png"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/templ.png"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/tsukuba_r.png"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/logo_in_clutter.png"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/logo.png"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/pic3.png"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/pic1.png"
    "/home/cbassuser/avnOpticalPointing/opencv-2.4.9/samples/cpp/letter-recognition.data"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "samples")

