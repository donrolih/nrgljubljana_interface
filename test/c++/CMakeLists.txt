# Copy h5 files to binary dir
file(GLOB_RECURSE all_h5_ref_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.ref.h5)
foreach(file ${all_h5_ref_files})
  configure_file(${file} ${file} COPYONLY)
endforeach()

# List of all tests, by CTest label: unit (no NRG run), solver (short NRG runs), long (manual)
#file(GLOB_RECURSE all_tests RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
set(unit_tests refreq_log.cpp hilbert_transform.cpp nrg_libraries.cpp nrg_adapters.cpp parameter_file.cpp)
set(solver_tests siam.cpp solver_hdf5.cpp)
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

# These tests include NRG Ljubljana headers directly
target_link_libraries(nrg_libraries nrgljubljana::workflow nrgljubljana::numerics)
target_link_libraries(nrg_adapters nrgljubljana::numerics)
