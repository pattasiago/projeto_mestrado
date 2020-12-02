cmake_minimum_required(VERSION 2.8.3)
project(opencv101)

IF(MSVC)
  ADD_DEFINITIONS(-DNOMINMAX)
  # Another MSVC quirk. http://stackoverflow.com/questions/6832666/lnk2019-when-including-asio-headers-solution-generated-with-cmake
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /MP") # For multi-threaded compilation on MSVC
ELSE()
  IF(APPLE)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++") # C++11 support on Apple clang
  ELSE()
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")# C++11 support on GCC
  ENDIF()
ENDIF()

find_package(catkin REQUIRED COMPONENTS
  cv_bridge
  roscpp
  rospy
  sensor_msgs
  std_msgs
  image_transport
  message_generation
  cmake_modules
)

find_package(OpenCV REQUIRED)
find_package(Eigen3 REQUIRED)

FIND_PACKAGE(OpenMP)
IF(OPENMP_FOUND)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  MESSAGE(STATUS "Found OpenMP")
ENDIF()

add_message_files(
  FILES
  imageParams.msg
  desvioParams.msg
  optimizerParams.msg
  optimizerResult.msg
)

generate_messages(
  DEPENDENCIES
  std_msgs
)

catkin_package(
#  INCLUDE_DIRS include
#  LIBRARIES opencv101
#  CATKIN_DEPENDS cv_bridge opencv2 roscpp rospy sensor_msgs std_msgs
#  DEPENDS system_lib
  CATKIN_DEPENDS message_runtime
)


include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${OpenCV_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
  src/optimizer
)


set(SOURCES src/drawCircle.cpp src/Lane.cpp)
add_executable(${PROJECT_NAME}_node ${SOURCES})
target_link_libraries(${PROJECT_NAME}_node
   ${catkin_LIBRARIES}
   ${OpenCV_LIBS}
)
add_dependencies(${PROJECT_NAME}_node ${${PROJECT_NAME}_EXPORTED_TARGETS})


#set(SOURCES2 src/visualize.cpp)
#add_executable(${PROJECT_NAME}_node_visualize ${SOURCES2})
#target_link_libraries(${PROJECT_NAME}_node_visualize
#   ${catkin_LIBRARIES}
#   ${OpenCV_LIBS}
#)
#add_dependencies(${PROJECT_NAME}_node_visualize ${${PROJECT_NAME}_EXPORTED_TARGETS})


#add_executable(image_converter src/image_converter.cpp)
#target_link_libraries(image_converter ${catkin_LIBRARIES} ${OpenCV_LIBRARIES})

add_executable(optimizerCaller src/optimizerCaller.cpp)
target_link_libraries(optimizerCaller ${catkin_LIBRARIES} )
add_dependencies(optimizerCaller ${PROJECT_NAME}_msgs_generate_messages_cpp optimizer/optimizer.c ${${PROJECT_NAME}_EXPORTED_TARGETS})

FIND_LIBRARY(meschach meschach.a PATHS ${PROJECT_SOURCE_DIR}/src/optimizer)
FIND_LIBRARY(optimizer optimizer.o PATHS ${PROJECT_SOURCE_DIR}/src/optimizer)

target_link_libraries(optimizerCaller ${PROJECT_SOURCE_DIR}/src/optimizer/donlp2.o)
target_link_libraries(optimizerCaller ${PROJECT_SOURCE_DIR}/src/optimizer/newx.o)
target_link_libraries(optimizerCaller ${PROJECT_SOURCE_DIR}/src/optimizer/user_eval.o)
target_link_libraries(optimizerCaller ${PROJECT_SOURCE_DIR}/src/optimizer/optimizer.o ${PROJECT_SOURCE_DIR}/src/optimizer/meschach.a)
