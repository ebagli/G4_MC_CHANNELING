#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(channeling)

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

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/channeling/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc ${PROJECT_SOURCE_DIR}/channeling/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/channeling/include/*.hh)

#----------------------------------------------------------------------------
# Install data 
#
file(GLOB _data RELATIVE ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/data/*.txt)
foreach(_dat ${_data})
     configure_file(
	${PROJECT_SOURCE_DIR}/${_dat}
     	${PROJECT_BINARY_DIR}/${_dat}
	COPYONLY
	)
endforeach()

file(GLOB _mac RELATIVE ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/mac/*.mac)
foreach(_mac ${_mac})
     configure_file(
	${PROJECT_SOURCE_DIR}/${_mac}
     	${PROJECT_BINARY_DIR}/${_mac}
	COPYONLY
	)
endforeach()


#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(channeling channeling.cc ${sources} ${headers})
message(${Geant4_LIBRARIES})
target_link_libraries(channeling ${Geant4_LIBRARIES})


#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
add_custom_target(channeling_custom DEPENDS channeling)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS channeling DESTINATION bin)


