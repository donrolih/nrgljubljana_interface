# Copy h5 files to binary dir
file(GLOB_RECURSE all_h5_ref_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.ref.h5)
foreach(file ${all_h5_ref_files})
  configure_file(${file} ${file} COPYONLY)
endforeach()

# List of all tests, by CTest label: unit (no NRG run), solver (short NRG runs), long (manual)
#file(GLOB_RECURSE all_tests RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
set(unit_tests refreq_log.cpp hilbert_transform.cpp nrg_libraries.cpp nrg_adapters.cpp parameter_file.cpp manifest.cpp solver_construction.cpp options_and_result.cpp param_rendering.cpp nrg_output.cpp broadening.cpp dmft_causality.cpp dmft_self_energy.cpp dmft_measurements.cpp dmft_lattice.cpp dmft_mixing.cpp dmft_mu_control.cpp dmft_convergence.cpp dmft_loop.cpp)
set(solver_tests siam.cpp solver_hdf5.cpp raw_results.cpp independent_solves.cpp deferred_broadening.cpp)
set(long_tests)
set(all_tests ${unit_tests} ${solver_tests} ${long_tests})

foreach(test ${all_tests})
  get_filename_component(test_name ${test} NAME_WE)
  get_filename_component(test_dir ${test} DIRECTORY)
  add_executable(${test_name} ${test})
  target_link_libraries(${test_name} ${PROJECT_NAME}::${PROJECT_NAME}_c ${PROJECT_NAME}_warnings gtest_main)
  set_property(TARGET ${test_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})
  add_test(NAME ${test_name} COMMAND ${test_name} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})
  set_property(TEST ${test_name} APPEND PROPERTY ENVIRONMENT
    PATH=${nrgljubljana_BINARY_DIR}/tools:$ENV{PATH}
    NRGIF_TEMPLATE_DIR=${PROJECT_SOURCE_DIR}/templates
  )
  # Run clang-tidy if found
  if(CLANG_TIDY_EXECUTABLE)
    set_target_properties(${test_name} PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXECUTABLE}")
  endif()
  # Run cppcheck if found
  if(CPPCHECK_EXECUTABLE)
    add_custom_command(
      TARGET ${test_name}
      COMMAND ${CPPCHECK_EXECUTABLE}
      --enable=warning,style,performance,portability
      --std=c++20
      --template=gcc
      --verbose
      --force
      --quiet
      ${CMAKE_CURRENT_SOURCE_DIR}/${test}
    )
  endif()
endforeach()

foreach(label unit solver long)
  foreach(test ${${label}_tests})
    get_filename_component(test_name ${test} NAME_WE)
    set_property(TEST ${test_name} APPEND PROPERTY LABELS ${label})
  endforeach()
endforeach()

# Tolerance of a reference X.ref.h5, recorded in X.ref.json (as in test/python): its precision is passed to
# test X as REF_PRECISION. Editing the JSON file reconfigures.
foreach(test ${all_tests})
  get_filename_component(test_name ${test} NAME_WE)
  set(ref_json ${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.ref.json)
  if(EXISTS ${ref_json})
    file(READ ${ref_json} ref_json_content)
    string(JSON ref_precision GET "${ref_json_content}" precision)
    target_compile_definitions(${test_name} PRIVATE REF_PRECISION=${ref_precision})
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${ref_json})
  endif()
endforeach()

# The raw results must be the same on every rank (own directory: the serial run writes the same files)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/raw_results_np2)
add_test(NAME raw_results_np2 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:raw_results> ${MPIEXEC_POSTFLAGS}
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/raw_results_np2)
set_property(TEST raw_results_np2 APPEND PROPERTY ENVIRONMENT
  PATH=${nrgljubljana_BINARY_DIR}/tools:$ENV{PATH}
  NRGIF_TEMPLATE_DIR=${PROJECT_SOURCE_DIR}/templates
)
set_property(TEST raw_results_np2 APPEND PROPERTY LABELS solver)

# Independent solvers, one per rank, run at the same time
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/independent_solves_np2)
add_test(NAME independent_solves_np2
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:independent_solves> ${MPIEXEC_POSTFLAGS}
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/independent_solves_np2)
set_property(TEST independent_solves_np2 APPEND PROPERTY ENVIRONMENT
  PATH=${nrgljubljana_BINARY_DIR}/tools:$ENV{PATH}
  NRGIF_TEMPLATE_DIR=${PROJECT_SOURCE_DIR}/templates
)
set_property(TEST independent_solves_np2 APPEND PROPERTY LABELS solver)

# The DMFT loop on 2 ranks: every rank holds rank 0's state (mock solver, no NRG run); own directory, since the serial
# run writes run directories of the same names
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/dmft_loop_np2)
add_test(NAME dmft_loop_np2 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:dmft_loop> --gtest_filter=Loop.AllRanksAgree:Loop.Resume:Loop.ResumeChecksParameters ${MPIEXEC_POSTFLAGS}
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/dmft_loop_np2)
set_property(TEST dmft_loop_np2 APPEND PROPERTY LABELS unit)

# Deferred broadening in a fresh process, from the file the first run writes
add_test(NAME deferred_broadening_fresh COMMAND deferred_broadening --gtest_filter=DeferredBroadeningFromFile.*
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_property(TEST deferred_broadening_fresh APPEND PROPERTY ENVIRONMENT
  PATH=${nrgljubljana_BINARY_DIR}/tools:$ENV{PATH}
  NRGIF_TEMPLATE_DIR=${PROJECT_SOURCE_DIR}/templates
)
set_property(TEST deferred_broadening_fresh APPEND PROPERTY LABELS solver)
set_tests_properties(deferred_broadening PROPERTIES FIXTURES_SETUP deferred_broadening_file)
set_tests_properties(deferred_broadening_fresh PROPERTIES FIXTURES_REQUIRED deferred_broadening_file)

# The broadening test compares with the broaden and kk tools
add_dependencies(broadening broaden kk)

# These tests include NRG Ljubljana headers directly
target_link_libraries(nrg_libraries nrgljubljana::workflow nrgljubljana::numerics)
target_link_libraries(nrg_adapters nrgljubljana::numerics)
target_link_libraries(dmft_causality nrgljubljana::numerics)
