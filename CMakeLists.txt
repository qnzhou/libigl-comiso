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
    #igl_copyleft::core
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
