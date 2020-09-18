#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(mantis)
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
option(WITH_GEANT4_ANALYSIS_USE_ROOT "Build example with Geant4 ROOT ANALYSIS" OFF)
#add_compile_options(-Wno-unused-variable -Wno-shadow -Wno-unused-parameter)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Find ROOT (required package)
if(WITH_GEANT4_ANALYSIS_USE_ROOT)
  list(APPEND CMAKE_MODULE_PATH $ENV{ROOTCMAKEPATH})
  find_package(ROOT QUIET)
  if(ROOT_FOUND)
    message(STATUS "ROOT found. --> Compiling with ROOT enabled")
    add_definitions(-DG4ANALYSIS_USE_ROOT)
  else()
    message(STATUS "G4 ROOT NOT Found. --> ROOT disabled")
  endif()
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${Geant4_INCLUDE_DIR}
                    ${ROOT_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(mantis mantis.cc ${sources} ${headers})
target_link_libraries(mantis ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})


#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build mantis. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(mantis_SCRIPTS
    mantis.in
    brems_distributions.root
    vis.mac
    input_spectrum.txt
    init_vis.mac
    vis_save.mac
  )

foreach(_script ${mantis_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS mantis DESTINATION bin)
