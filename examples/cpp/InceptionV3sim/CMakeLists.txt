cmake_minimum_required(VERSION 3.10)

project(FlexFlowExample_Inceptionsim)
set(project_target inceptionsim)

set(CPU_SRC
  ${FLEXFLOW_CPP_DRV_SRC}
  inception.cc
  inception.h)

set(GPU_SRC
  inception.cu)

cuda_add_executable(${project_target} ${CPU_SRC} ${GPU_SRC})
target_include_directories(${project_target} PRIVATE ${FLEXFLOW_INCLUDE_DIRS} ${CMAKE_INSTALL_INCLUDEDIR})
target_link_libraries(${project_target} -Wl,--whole-archive flexflow -Wl,--no-whole-archive ${FLEXFLOW_EXT_LIBRARIES})

set(BIN_DEST "bin")
install(TARGETS ${project_target} DESTINATION ${BIN_DEST})
