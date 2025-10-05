cmake_minimum_required(VERSION 3.14)

project(igl_comiso)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
include(comiso)
include(libigl)

file(GLOB INC_FILES ${CMAKE_CURRENT_LIST_DIR}/include/*.h)
file(GLOB SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/include/*.impl)

if(LIBIGL_USE_STATIC_LIBRARY)
    set(igl_comiso_type STATIC)
    set(igl_comiso_scope PUBLIC)
else()
    set(igl_comiso_type INTERFACE)
    set(igl_comiso_scope INTERFACE)
endif()

add_library(igl_comiso ${igl_comiso_type} ${INC_FILES} ${SRC_FILES})
target_compile_features(igl_comiso ${igl_comiso_scope} cxx_std_14)
target_include_directories(igl_comiso ${igl_comiso_scope} ${CMAKE_CURRENT_LIST_DIR}/include)
target_link_libraries(igl_comiso ${igl_comiso_scope}
    igl::core
    CoMISo::CoMISo
)
add_library(igl_copyleft::comiso ALIAS igl_comiso)

if(LIBIGL_POSITION_INDEPENDENT_CODE)
    set_target_properties(igl_comiso PROPERTIES INTERFACE_POSITION_INDEPENDENT_CODE ON)
    if(LIBIGL_USE_STATIC_LIBRARY)
        set_target_properties(igl_comiso PROPERTIES POSITION_INDEPENDENT_CODE ON)
    endif()
endif()

if(LIBIGL_BUILD_TUTORIALS)
    add_executable(505_MIQ ${CMAKE_CURRENT_LIST_DIR}/tutorial/505_MIQ/main.cpp)
    target_link_libraries(505_MIQ PRIVATE igl_copyleft::comiso igl::glfw igl::tutorial_data)

    add_executable(506_FrameField ${CMAKE_CURRENT_LIST_DIR}/tutorial/506_FrameField/main.cpp)
    target_link_libraries(506_FrameField PRIVATE igl_copyleft::comiso igl::glfw igl::tutorial_data)
endif()

option(LIBIGL_COMISO_BUILD_TESTS "Build the libigl comiso tests" OFF)
if (LIBIGL_COMISO_BUILD_TESTS)
    include(CTest)
    enable_testing()
    include(catch2)

    file(GLOB TEST_FILES "${CMAKE_CURRENT_LIST_DIR}/tests/include/igl/copyleft/comiso/*.cpp")
    add_executable(igl_comiso_tests ${TEST_FILES})
    target_include_directories(igl_comiso_tests PRIVATE ${CMAKE_CURRENT_LIST_DIR}/tests)
    target_link_libraries(igl_comiso_tests PRIVATE
        igl_copyleft::comiso
        igl::tests_data
        Catch2::Catch2WithMain)
endif()

if(SKBUILD)
    if(APPLE)
        set(CMAKE_INSTALL_RPATH @loader_path)
    elseif(UNIX)
        set(CMAKE_INSTALL_RPATH $ORIGIN)
    endif()

    include(nanobind)
    set(PY_SRC_FILE "${PROJECT_SOURCE_DIR}/python/binding.cpp")
    nanobind_add_module(pyigl_comiso NB_STATIC ${PY_SRC_FILE})
    target_compile_features(pyigl_comiso PRIVATE cxx_std_20)
    target_link_libraries(pyigl_comiso PUBLIC igl_copyleft::comiso)

    install(TARGETS pyigl_comiso CoMISo LIBRARY
        COMPONENT LIBIGL_COMISO_Python_Runtime
        DESTINATION ${SKBUILD_PLATLIB_DIR}/igl_comiso)
    add_library(igl_copyleft::comiso::python ALIAS pyigl_comiso)

endif()
