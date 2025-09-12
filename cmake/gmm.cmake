if (TARGET gmm::gmm)
    return()
endif()

message(STATUS "Third-party (external): importing 'gmm'")

include(CPM)
CPMAddPackage(
    NAME gmm
    URL http://download-mirror.savannah.gnu.org/releases/getfem/stable/gmm-5.4.4.tar.gz
)

add_library(gmm INTERFACE)
target_include_directories(gmm INTERFACE ${gmm_SOURCE_DIR}/include)
#target_compile_definitions(gmm INTERFACE -DGMM_USES_SUPERLU)

add_library(gmm::gmm ALIAS gmm)
