if(TARGET CoMISo::CoMISo)
    return()
endif()

message(STATUS "Third-party: creating target 'CoMISo::CoMISo'")

set(COMISO_NO_INSTALL On CACHE BOOL "Do not install CoMISo when installing the parent project" FORCE)
set(COMISO_CXX_STANDARD 17 CACHE STRING "C++ standard to use" FORCE)
# We only need gmm's headers (the igl comiso wrappers include <gmm/gmm_kernel.h>),
# not CoMISo's optional GMM features. Disable CoMISo's own GMM lookup and provide
# the headers ourselves below.
set(COMISO_ENABLE_GMM Off CACHE BOOL "Build CoMISo with GMM" FORCE)

include(gmm)

include(CPM)
CPMAddPackage(
    NAME comiso
    GIT_REPOSITORY https://graphics.rwth-aachen.de:9000/CoMISo/CoMISo.git
    GIT_TAG fe7feead486a4a825f5dc95d68b98249366b66b8
)
# NOTE: upstream now defines the `CoMISo::CoMISo` alias itself and exports its own
# PUBLIC include directories (src/ and include/), so the previous manual alias and
# the .hh header-copy hack are no longer needed.
set_target_properties(CoMISo PROPERTIES
    CXX_STANDARD ${COMISO_CXX_STANDARD}
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO)

# CoMISo's cmake-library sets -DINCLUDE_TEMPLATES globally, which pulls in
# ExactConstraintProjection_impl.hh. That file uses assert() without including
# <cassert>, relying on Eigen 3 to leak it. Eigen 5 (conda-forge) does not.
target_compile_options(CoMISo PRIVATE "-include" "cassert")
# Provide gmm headers to CoMISo and, transitively, to igl_comiso.
target_link_libraries(CoMISo PUBLIC gmm::gmm)

set_target_properties(CoMISo PROPERTIES FOLDER ThirdParty)
