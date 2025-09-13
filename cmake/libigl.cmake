if(TARGET igl::core)
    return()
endif()

message(STATUS "Third-party: creating target 'igl::core'")

#set(LIBIGL_COPYLEFT_CGAL On CACHE BOOL "Build CGAL support in libigl")
#set(LIBIGL_COPYLEFT_COMISO On CACHE BOOL "Build Comiso support in libigl")

include(CPM)
CPMAddPackage(
    libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG v2.6.0
    OPTIONS LIBIGL_INSTALL OFF
)

if(LIBIGL_BUILD_TUTORIALS)
    include(${libigl_SOURCE_DIR}/cmake/recipes/external/libigl_tutorial_data.cmake)
endif()

if (LIBIGL_COMISO_BUILD_TESTS)
    include(${libigl_SOURCE_DIR}/cmake/recipes/external/libigl_tests_data.cmake)
endif()
