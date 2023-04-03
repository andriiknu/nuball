#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(Nuball)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO)
include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(mynuball mynuball.cc ${sources} ${headers})
target_link_libraries(mynuball ${Geant4_LIBRARIES})
target_link_libraries(mynuball ${ROOT_LIBRARIES})



#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build ED. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
file(GLOB Nuball_MAC ${PROJECT_SOURCE_DIR}/*.mac)
file(GLOB Nuball_DAT ${PROJECT_SOURCE_DIR}/*.dat)

file(COPY ${Nuball_MAC} DESTINATION ${PROJECT_BINARY_DIR})
file(COPY ${Nuball_DAT} DESTINATION ${PROJECT_BINARY_DIR})

set(EXAMPLEED_SCRIPTS
  # init_vis.mac
  # run.mac
  run.png
  # vis.mac
  # CloverBGOPosition.dat
  # CloverPosition.dat
  # PhaseIBGOPosition.dat
  # PhaseIPosition.dat
  # LaBr3Position.dat
  )
# list(APPEND EXAMPLEED_SCRIPTS
#   # ${NuballMACROS}
#   # ${NuballData}
#   )


foreach(_script ${EXAMPLEED_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()




#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
add_custom_target(Nuball DEPENDS mynuball)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS mynuball DESTINATION bin)


