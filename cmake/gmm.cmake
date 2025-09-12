if (TARGET gmm::gmm)
    return()
endif()

message(STATUS "Third-party (external): importing 'gmm'")

include(CPM)
CPMAddPackage(
    NAME gmm
    URL http://download-mirror.savannah.gnu.org/releases/getfem/stable/gmm-5.4.4.tar.gz
    URL_HASH SHA256=15eb1943011b92665aab3b02ecf3cede1cf89ea15a9006f81f2ba2cd662aa02b
)

add_library(gmm INTERFACE)
target_include_directories(gmm INTERFACE ${gmm_SOURCE_DIR}/include)
#target_compile_definitions(gmm INTERFACE -DGMM_USES_SUPERLU)

add_library(gmm::gmm ALIAS gmm)
