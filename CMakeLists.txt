cmake_minimum_required(VERSION 3.20)
set(CMAKE_OSX_DEPLOYMENT_TARGET "10.15" CACHE STRING "Minimum OS X deployment version")

project(CASMcode_configuration VERSION 2.0.0 LANGUAGES CXX)

# set CMAKE_INSTALL_X variables
include(GNUInstallDirs)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# try to use ccache
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
endif()

##############################################
## Find dependencies

# Should find ZLIB::ZLIB
find_package(ZLIB)

# Find CASM
if(NOT DEFINED CASM_PREFIX)
  message(STATUS "CASM_PREFIX not defined")
  # try to find Python
  find_package (Python COMPONENTS Interpreter Development)
  if(DEFINED Python_EXECUTABLE)
    # if Python found, obtain CASM_PREFIX from the libcasm.casmglobal
    message(STATUS "found Python_EXECUTABLE: ${Python_EXECUTABLE}")
    message(STATUS "checking for libcasm-global")
    execute_process(
      COMMAND pip show libcasm-global
      RESULT_VARIABLE EXIT_CODE
      OUTPUT_QUIET
    )
    if (${EXIT_CODE} EQUAL 0)
      message(STATUS "found libcasm-global")
      execute_process(COMMAND ${Python_EXECUTABLE} -m libcasm.casmglobal --prefix
                      OUTPUT_VARIABLE CASM_PREFIX_RAW)
      string(STRIP ${CASM_PREFIX_RAW} CASM_PREFIX)
      message(STATUS "CASM_PREFIX: ${CASM_PREFIX}")
    else()
      message(STATUS "did not find libcasm-global")
    endif()
  endif()
endif()
if(DEFINED CASM_PREFIX)
  set(CASMcode_global_ROOT ${CASM_PREFIX}/share/CASMcode_global/cmake)
  set(CASMcode_crystallography_ROOT ${CASM_PREFIX}/share/CASMcode_crystallography/cmake)
  set(CASMcode_composition_ROOT ${CASM_PREFIX}/share/CASMcode_composition/cmake)
  set(CASMcode_clexulator_ROOT ${CASM_PREFIX}/share/CASMcode_clexulator/cmake)
  set(CASMcode_configuration_ROOT ${CASM_PREFIX}/share/CASMcode_configuration/cmake)
  set(CASMcode_monte_ROOT ${CASM_PREFIX}/share/CASMcode_monte/cmake)
endif()

find_package(CASMcode_global)
if(NOT CASMcode_global_FOUND)
  message(FATAL_ERROR "CMake failed to find CASMcode_global")
endif()
# if successful, we have CASM::casm_global

find_package(CASMcode_crystallography)
if(NOT CASMcode_crystallography_FOUND)
  message(FATAL_ERROR "CMake failed to find CASMcode_crystallography")
endif()
# if successful, we have CASM::casm_crystallography

find_package(CASMcode_composition)
if(NOT CASMcode_composition_FOUND)
  message(FATAL_ERROR "CMake failed to find CASMcode_composition")
endif()
# if successful, we have CASM::casm_composition

find_package(CASMcode_clexulator)
if(NOT CASMcode_clexulator_FOUND)
  message(FATAL_ERROR "CMake failed to find CASMcode_clexulator")
endif()
# if successful, we have CASM::casm_clexulator

find_package(CASMcode_configuration)
if(NOT CASMcode_configuration_FOUND)
  message(FATAL_ERROR "CMake failed to find CASMcode_configuration")
endif()
# if successful, we have CASM::casm_configuration

find_package(CASMcode_monte)
if(NOT CASMcode_monte_FOUND)
  message(FATAL_ERROR "CMake failed to find CASMcode_monte")
endif()
# if successful, we have CASM::casm_monte


# if no user CMAKE_INSTALL_PREFIX, use CASM_PREFIX if it exists
IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  if(DEFINED CASM_PREFIX)
    message(STATUS "CMAKE_INSTALL_PREFIX initialized to default, so updating CMAKE_INSTALL_PREFIX to CASM_PREFIX")
    set(CMAKE_INSTALL_PREFIX ${CASM_PREFIX} CACHE PATH "set CMAKE_INSTALL_PREFIX to CASM_PREFIX" FORCE)
    message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
  endif()
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)


##############################################
## Build libcasm_clexmonte

# create libcasm_clexmonte
set(
  libcasm_clexmonte_HEADERS
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/canonical/canonical.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/canonical/canonical_impl.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/canonical/canonical_json_io.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/definitions.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/events/AllowedEventList.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/events/CompleteEventList.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/events/EventGroup.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/events/ImpactTable.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/events/event_data.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/events/event_methods.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/events/io/json/event_data_json_io.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/events/io/stream/EventState_stream_io.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/events/lotto.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/events/lotto/event_rate_tree.hpp
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/events/lotto/event_rate_tree_impl.hpp
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/events/lotto/event_selector.hpp
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/events/lotto/random.hpp
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/events/lotto/rejection.hpp
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/events/lotto/rejection_free.hpp
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/events/lotto/sum_tree.hpp
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/events/lotto/sum_tree_impl.hpp
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/events/state_graph.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/kinetic/io/json/EventState_json_io.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/kinetic/io/stream/EventState_stream_io.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/kinetic/kinetic.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/kinetic/kinetic_events.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/kinetic/kinetic_impl.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/kinetic/kinetic_json_io.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/methods/kinetic_monte_carlo.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/methods/occupation_metropolis.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/misc/MapLike.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/misc/Matrix3lCompare.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/misc/check_params.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/misc/diffusion_calculations.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/misc/eigen.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/misc/parse_array.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/misc/subparse_from_file.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/misc/to_json.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/monte_calculator/BaseMonteCalculator.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/monte_calculator/BaseMonteEventData.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/monte_calculator/KineticCalculator.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/monte_calculator/MonteCalculator.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/monte_calculator/MonteEventData.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/monte_calculator/StateData.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/monte_calculator/analysis_functions.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/monte_calculator/io/json/MonteCalculator_json_io.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/monte_calculator/kinetic_events.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/monte_calculator/kinetic_sampling_functions.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/monte_calculator/modifying_functions.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/monte_calculator/sampling_functions.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/monte_calculator/selected_event_functions.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/nfold/nfold.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/nfold/nfold_events.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/nfold/nfold_impl.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/nfold/nfold_json_io.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/run/ConfigGenerator.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/run/FixedConfigGenerator.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/run/IncrementalConditionsStateGenerator.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/run/RunData.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/run/StateGenerator.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/run/StateModifyingFunction.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/run/analysis_functions.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/run/covariance_functions.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/run/functions.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/run/io/RunParams.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/run/io/json/ConfigGenerator_json_io.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/run/io/json/RunData_json_io.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/run/io/json/RunParams_json_io.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/run/io/json/RunParams_json_io_impl.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/run/io/json/StateGenerator_json_io.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/run/io/json/StateGenerator_json_io_impl.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/run/io/json/parse_and_run_series.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/semigrand_canonical/calculator.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/semigrand_canonical/calculator_impl.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/semigrand_canonical/conditions.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/semigrand_canonical/event_generator.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/semigrand_canonical/json_io.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/semigrand_canonical/potential.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/state/Conditions.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/state/Configuration.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/state/CorrMatchingPotential.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/state/LocalOrbitCompositionCalculator.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/state/enforce_composition.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/state/io/json/CorrMatchingPotential_json_io.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/state/io/json/State_json_io.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/state/io/json/parse_conditions.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/state/io/json/parse_conditions_impl.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/state/kinetic_sampling_functions.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/state/make_conditions.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/state/modifying_functions.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/state/sampling_functions.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/system/System.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/system/io/json/System_json_io.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/system/io/json/system_data_json_io.hh
  ${PROJECT_SOURCE_DIR}/include/casm/clexmonte/system/system_data.hh
)
set(
  libcasm_clexmonte_SOURCES
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/canonical/canonical.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/events/AllowedEventList.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/events/CompleteEventList.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/events/EventGroup.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/events/ImpactTable.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/events/event_data.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/events/event_methods.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/events/io/json/event_data_json_io.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/events/io/stream/EventState_stream_io.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/events/state_graph.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/kinetic/io/json/EventState_json_io.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/kinetic/io/stream/EventState_stream_io.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/kinetic/kinetic.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/kinetic/kinetic_events.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/monte_calculator/BaseMonteCalculator.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/monte_calculator/BaseMonteEventData.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/monte_calculator/CanonicalCalculator.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/monte_calculator/KineticCalculator.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/monte_calculator/MonteCalculator.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/monte_calculator/MonteEventData.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/monte_calculator/SemiGrandCanonicalCalculator.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/monte_calculator/StateData.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/monte_calculator/analysis_functions.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/monte_calculator/io/json/MonteCalculator_json_io.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/monte_calculator/kinetic_events.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/monte_calculator/kinetic_sampling_functions.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/monte_calculator/sampling_functions.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/monte_calculator/selected_event_functions.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/nfold/nfold.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/nfold/nfold_events.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/run/io/convariance_functions.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/run/io/json/ConfigGenerator_json_io.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/run/io/json/RunParams_json_io.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/run/io/json/StateGenerator_json_io.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/semigrand_canonical/calculator.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/semigrand_canonical/potential.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/state/Conditions.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/state/CorrMatchingPotential.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/state/LocalOrbitCompositionCalculator.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/state/io/json/CorrMatchingPotential_json_io.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/state/io/json/State_json_io.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/state/io/json/parse_conditions.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/state/make_conditions.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/system/System.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/system/io/json/System_json_io.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/system/io/json/system_data_json_io.cc
  ${PROJECT_SOURCE_DIR}/src/casm/clexmonte/system/system_data.cc
)
add_library(casm_clexmonte SHARED ${libcasm_clexmonte_SOURCES})
target_include_directories(casm_clexmonte
  PUBLIC
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/casm/external>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/casm/external/gzstream>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/casm/external>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/casm/external/gzstream>
)
target_compile_options(casm_clexmonte
  PUBLIC
    "-DCASM_CLEXMONTE_TXT_VERSION=\"${CMAKE_PROJECT_VERSION}\""
    -DEIGEN_DEFAULT_DENSE_INDEX_TYPE=long
    -DGZSTREAM_NAMESPACE=gz
)
target_link_libraries(casm_clexmonte
  ZLIB::ZLIB
  ${CMAKE_DL_LIBS}
  CASM::casm_global
  CASM::casm_crystallography
  CASM::casm_composition
  CASM::casm_clexulator
  CASM::casm_configuration
  CASM::casm_monte
)
if(APPLE)
  set_target_properties(
    casm_clexmonte PROPERTIES INSTALL_RPATH "@loader_path")
else()
  set_target_properties(
    casm_clexmonte PROPERTIES INSTALL_RPATH "$ORIGIN")
endif()


##############################################
## Install libcasm_clexmonte

# install header files in <prefix>/libcasm/include/,
# while preserving directory structure
foreach ( filevar ${libcasm_clexmonte_HEADERS} )
  file(RELATIVE_PATH relfile ${PROJECT_SOURCE_DIR}/include/ ${filevar})
  get_filename_component( reldir ${relfile} DIRECTORY )
  install( FILES ${filevar} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${reldir} )
endforeach()

# install libcasm_clexmonte in <prefix>/libcasm/lib/
install(
  TARGETS casm_clexmonte
  EXPORT CASMcode_clexmonteTargets
  DESTINATION lib)

##############################################
## Python extensions

# The CMake package config and target files are installed under the Python
# package root. This is necessary to ensure that all the relative paths in the
# helloTargets.cmake resolve correctly. It also provides encapsulation.
#
# The actual path used must be selected so that consuming projects can locate it
# via `find_package`. To support finding CMake packages in the Python package
# prefix, using `find_package`s default search path of
# `<prefix>/<name>/share/<name>*/cmake/` is reasonable. Adding the Python
# package installation prefix to CMAKE_PREFIX_PATH in combination with this path
# will allow `find_package` to find this package and any other package installed
# via a Python package if the CMake and Python packages are named the same.
set(CASM_CMAKE_PACKAGE_INSTALL_SUBDIR "share/CASMcode_clexmonte/cmake")

install(
  EXPORT CASMcode_clexmonteTargets
  NAMESPACE CASM::
  DESTINATION ${CASM_CMAKE_PACKAGE_INSTALL_SUBDIR})

include(CMakePackageConfigHelpers)

write_basic_package_version_file(
  CASMcode_clexmonteConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMinorVersion)

configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/CASMcode_clexmonteConfig.cmake.in" CASMcode_clexmonteConfig.cmake
  INSTALL_DESTINATION ${CASM_CMAKE_PACKAGE_INSTALL_SUBDIR})

install(FILES "${PROJECT_BINARY_DIR}/CASMcode_clexmonteConfig.cmake"
              "${PROJECT_BINARY_DIR}/CASMcode_clexmonteConfigVersion.cmake"
        DESTINATION ${CASM_CMAKE_PACKAGE_INSTALL_SUBDIR})

# We are using the SKBUILD variable, which is defined when scikit-build is
# running the CMake build, to control building the Python wrapper. This allows
# the C++ project to be installed, standalone, when using the standard CMake
# build flow.
if(DEFINED SKBUILD)

  # call pybind11-config to obtain the root of the cmake package
  execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pybind11 --cmakedir
                  OUTPUT_VARIABLE pybind11_ROOT_RAW)
  string(STRIP ${pybind11_ROOT_RAW} pybind11_ROOT)
  find_package(pybind11)

  # The extension modules must load:
  # - the casm_global library
  # - the casm_crystallography library
  # - the casm_composition library
  # - the casm_clexulator library
  # - the casm_configuration library
  # - the casm_monte library
  # They can be found by setting a relative rpath

  ### libcasm.clexmonte._clexmonte_system ###
  pybind11_add_module(_clexmonte_system MODULE
                      "${PROJECT_SOURCE_DIR}/python/src/clexmonte_system.cpp")
  target_link_libraries(_clexmonte_system PRIVATE
    CASM::casm_global
    CASM::casm_crystallography
    CASM::casm_composition
    CASM::casm_clexulator
    CASM::casm_configuration
    CASM::casm_monte
    casm_clexmonte
  )
  install(TARGETS _clexmonte_system DESTINATION clexmonte)
  if(APPLE)
    set_target_properties(
      _clexmonte_system PROPERTIES INSTALL_RPATH "@loader_path/../lib")
  else()
    set_target_properties(
      _clexmonte_system PROPERTIES INSTALL_RPATH "$ORIGIN/../lib")
  endif()

  ### libcasm.clexmonte._clexmonte_state ###
  pybind11_add_module(_clexmonte_state MODULE
                      "${PROJECT_SOURCE_DIR}/python/src/clexmonte_state.cpp")
  target_link_libraries(_clexmonte_state PRIVATE
    CASM::casm_global
    CASM::casm_crystallography
    CASM::casm_composition
    CASM::casm_clexulator
    CASM::casm_configuration
    CASM::casm_monte
    casm_clexmonte
  )
  install(TARGETS _clexmonte_state DESTINATION clexmonte)
  if(APPLE)
    set_target_properties(
      _clexmonte_state PROPERTIES INSTALL_RPATH "@loader_path/../lib")
  else()
    set_target_properties(
      _clexmonte_state PROPERTIES INSTALL_RPATH "$ORIGIN/../lib")
  endif()

  ### libcasm.clexmonte._clexmonte_run_management ###
  pybind11_add_module(_clexmonte_run_management MODULE
                      "${PROJECT_SOURCE_DIR}/python/src/clexmonte_run_management.cpp")
  target_link_libraries(_clexmonte_run_management PRIVATE
    CASM::casm_global
    CASM::casm_crystallography
    CASM::casm_composition
    CASM::casm_clexulator
    CASM::casm_configuration
    CASM::casm_monte
    casm_clexmonte
  )
  install(TARGETS _clexmonte_run_management DESTINATION clexmonte)
  if(APPLE)
    set_target_properties(
      _clexmonte_run_management PROPERTIES INSTALL_RPATH "@loader_path/../lib")
  else()
    set_target_properties(
      _clexmonte_run_management PROPERTIES INSTALL_RPATH "$ORIGIN/../lib")
  endif()

  ### libcasm.clexmonte._clexmonte_functions ###
  pybind11_add_module(_clexmonte_functions MODULE
                      "${PROJECT_SOURCE_DIR}/python/src/clexmonte_functions.cpp")
  target_link_libraries(_clexmonte_functions PRIVATE
    CASM::casm_global
    CASM::casm_crystallography
    CASM::casm_composition
    CASM::casm_clexulator
    CASM::casm_configuration
    CASM::casm_monte
    casm_clexmonte
  )
  install(TARGETS _clexmonte_functions DESTINATION clexmonte)
  if(APPLE)
    set_target_properties(
      _clexmonte_functions PROPERTIES INSTALL_RPATH "@loader_path/../lib")
  else()
    set_target_properties(
      _clexmonte_functions PROPERTIES INSTALL_RPATH "$ORIGIN/../lib")
  endif()

  ### libcasm.clexmonte._clexmonte_monte_calculator ###
  pybind11_add_module(_clexmonte_monte_calculator MODULE
                      "${PROJECT_SOURCE_DIR}/python/src/clexmonte_monte_calculator.cpp")
  target_link_libraries(_clexmonte_monte_calculator PRIVATE
    CASM::casm_global
    CASM::casm_crystallography
    CASM::casm_composition
    CASM::casm_clexulator
    CASM::casm_configuration
    CASM::casm_monte
    casm_clexmonte
  )
  install(TARGETS _clexmonte_monte_calculator DESTINATION clexmonte)
  if(APPLE)
    set_target_properties(
      _clexmonte_monte_calculator PROPERTIES INSTALL_RPATH "@loader_path/../lib")
  else()
    set_target_properties(
      _clexmonte_monte_calculator PROPERTIES INSTALL_RPATH "$ORIGIN/../lib")
  endif()
endif()
