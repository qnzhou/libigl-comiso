if(TARGET CoMISo::CoMISo)
    return()
endif()

message(STATUS "Third-party: creating target 'CoMISo::CoMISo'")

set(COMISO_NO_INSTALL On CACHE BOOL "Do not install CoMISo when installing the parent project" FORCE)
set(COMISO_CXX_STANDARD 17 CACHE STRING "C++ standard to use" FORCE)

include(gmm)

include(CPM)
CPMAddPackage(
    NAME comiso
    GIT_REPOSITORY https://graphics.rwth-aachen.de:9000/CoMISo/CoMISo.git
    GIT_TAG 0f9781a0433684b3b934732fd846fe07e522ca4d
)
set_target_properties(CoMISo PROPERTIES
    CXX_STANDARD ${COMISO_CXX_STANDARD}
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO)
target_include_directories(CoMISo PUBLIC ${comiso_SOURCE_DIR}/include)


add_library(CoMISo::CoMISo ALIAS CoMISo)

# Copy .hh headers into a subfolder `CoMISo/`
file(GLOB_RECURSE INC_FILES "${comiso_SOURCE_DIR}/*.hh")
set(output_folder "${CMAKE_CURRENT_BINARY_DIR}/CoMISo/include/CoMISo/")
message(VERBOSE "Copying CoMISo headers to '${output_folder}'")
foreach(filepath IN ITEMS ${INC_FILES})
    file(RELATIVE_PATH filename "${comiso_SOURCE_DIR}" ${filepath})
    configure_file(${filepath} "${output_folder}/${filename}" COPYONLY)
endforeach()

target_include_directories(CoMISo PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/CoMISo/include)
target_link_libraries(CoMISo PUBLIC gmm::gmm)

set_target_properties(CoMISo PROPERTIES FOLDER ThirdParty)
