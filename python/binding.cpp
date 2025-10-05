#include <igl/copyleft/comiso/frame_field.h>
#include <igl/copyleft/comiso/miq.h>
#include <igl/copyleft/comiso/nrosy.h>

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <Eigen/Core>

namespace nb = nanobind;

NB_MODULE(pyigl_comiso, m)
{
    using namespace nb::literals;

    using VectorXd = Eigen::VectorXd;
    using VectorXi = Eigen::VectorXi;
    using MatrixXd = Eigen::MatrixXd;
    using MatrixXi = Eigen::MatrixXi;

    m.def(
        "nrosy",
        [](const MatrixXd& V,
           const MatrixXi& F,
           const VectorXi& b,
           const MatrixXd& bc,
           const VectorXi b_soft,
           const VectorXd w_soft,
           const MatrixXd bc_soft,
           int N,
           double soft) {
            MatrixXd R;
            VectorXd S;
            igl::copyleft::comiso::nrosy(V, F, b, bc, b_soft, w_soft, bc_soft, N, soft, R, S);
            return std::make_tuple(R, S);
        },
        "V"_a,
        "F"_a,
        "b"_a,
        "bc"_a,
        "b_soft"_a = VectorXi(0),
        "w_soft"_a = VectorXd(0),
        "bc_soft"_a = MatrixXd(0, 3),
        "N"_a = 4,
        "soft"_a = 0.5,
        R"(Generate a N-RoSy field from a sparse set of constraints.

:param V: #V by 3 list of mesh vertex coordinates
:param F: #F by 3 list of mesh faces (must be triangles)
:param b: #B by 1 list of constrained face indices
:param bc: #B by 3 list of representative vectors for the constrained faces
:param b_soft: #S by 1 b for soft constraints
:param w_soft: #S by 1 weight for the soft constraints (0-1)
:param bc_soft: #S by 3 bc for soft constraints
:param N: the degree of the N-RoSy vector field
:param soft: the strength of the soft constraints w.r.t. smoothness (0 -> smoothness only, 1->constraints only)
:returns: R, S where R is #F by 3 the representative vectors of the interpolated field and S is #V by 1 the singularity index for each vertex (0 = regular))");

    // MIQ - Mixed Integer Quadrangulation (simple version)
    m.def(
        "miq",
        [](const MatrixXd& V,
           const MatrixXi& F,
           const MatrixXd& PD1,
           const MatrixXd& PD2,
           double gradientSize,
           double stiffness,
           bool directRound,
           unsigned int iter,
           unsigned int localIter,
           bool doRound,
           bool singularityRound,
           const std::vector<int>& roundVertices,
           const std::vector<std::vector<int>>& hardFeatures) {
            MatrixXd UV;
            MatrixXi FUV;
            igl::copyleft::comiso::miq(
                V,
                F,
                PD1,
                PD2,
                UV,
                FUV,
                gradientSize,
                stiffness,
                directRound,
                iter,
                localIter,
                doRound,
                singularityRound,
                roundVertices,
                hardFeatures);
            return std::make_tuple(UV, FUV);
        },
        "V"_a,
        "F"_a,
        "PD1"_a,
        "PD2"_a,
        "gradientSize"_a = 30.0,
        "stiffness"_a = 5.0,
        "directRound"_a = false,
        "iter"_a = 5,
        "localIter"_a = 5,
        "doRound"_a = true,
        "singularityRound"_a = true,
        "roundVertices"_a = std::vector<int>(),
        "hardFeatures"_a = std::vector<std::vector<int>>(),
        R"(Global seamless parametrization aligned with a given per-face Jacobian.

Based on "Mixed-Integer Quadrangulation" by D. Bommes, H. Zimmer, L. Kobbelt
ACM SIGGRAPH 2009, Article No. 77

:param V: #V by 3 list of mesh vertex 3D positions
:param F: #F by 3 list of faces indices in V
:param PD1: #F by 3 first line of the Jacobian per triangle
:param PD2: #F by 3 second line of the Jacobian per triangle (if empty, will be orthogonal to PD1)
:param gradientSize: global scaling for the gradient (controls the quads resolution)
:param stiffness: weight for the stiffness iterations (reserved but not used)
:param directRound: greedily round all integer variables at once (improves speed but lowers quality)
:param iter: stiffness iterations (0 = no stiffness)
:param localIter: number of local iterations for the integer rounding
:param doRound: enables the integer rounding (disabling useful for debugging)
:param singularityRound: round singularities' coordinates to nearest integers
:param roundVertices: additional vertices that should be snapped to integer coordinates
:param hardFeatures: pairs of vertices that belong to edges that should be snapped to integer coordinates
:returns: UV, FUV where UV is #UV by 2 list of vertices in 2D and FUV is #FUV by 3 list of face indices in UV)");

    // MIQ advanced version with pre-combed bisectors
    m.def(
        "miq_precombed",
        [](const MatrixXd& V,
           const MatrixXi& F,
           const MatrixXd& PD1_combed,
           const MatrixXd& PD2_combed,
           const Eigen::Matrix<int, Eigen::Dynamic, 3>& mismatch,
           const Eigen::Matrix<int, Eigen::Dynamic, 1>& singular,
           const Eigen::Matrix<int, Eigen::Dynamic, 3>& seams,
           double gradientSize,
           double stiffness,
           bool directRound,
           unsigned int iter,
           unsigned int localIter,
           bool doRound,
           bool singularityRound,
           const std::vector<int>& roundVertices,
           const std::vector<std::vector<int>>& hardFeatures) {
            MatrixXd UV;
            MatrixXi FUV;
            igl::copyleft::comiso::miq(
                V,
                F,
                PD1_combed,
                PD2_combed,
                mismatch,
                singular,
                seams,
                UV,
                FUV,
                gradientSize,
                stiffness,
                directRound,
                iter,
                localIter,
                doRound,
                singularityRound,
                roundVertices,
                hardFeatures);
            return std::make_tuple(UV, FUV);
        },
        "V"_a,
        "F"_a,
        "PD1_combed"_a,
        "PD2_combed"_a,
        "mismatch"_a,
        "singular"_a,
        "seams"_a,
        "gradientSize"_a = 30.0,
        "stiffness"_a = 5.0,
        "directRound"_a = false,
        "iter"_a = 5,
        "localIter"_a = 5,
        "doRound"_a = true,
        "singularityRound"_a = true,
        "roundVertices"_a = std::vector<int>(),
        "hardFeatures"_a = std::vector<std::vector<int>>(),
        R"(MIQ helper function with pre-combed bisectors for an already cut mesh.

:param V: #V by 3 list of mesh vertex 3D positions
:param F: #F by 3 list of faces indices in V
:param PD1_combed: #F by 3 first combed Jacobian
:param PD2_combed: #F by 3 second combed Jacobian
:param mismatch: #F by 3 list of per-corner integer PI/2 rotations
:param singular: #V list of flag that denotes if a vertex is singular or not
:param seams: #F by 3 list of per-corner flag that denotes seams
:param gradientSize: global scaling for the gradient (controls the quads resolution)
:param stiffness: weight for the stiffness iterations (reserved but not used)
:param directRound: greedily round all integer variables at once (improves speed but lowers quality)
:param iter: stiffness iterations (0 = no stiffness)
:param localIter: number of local iterations for the integer rounding
:param doRound: enables the integer rounding (disabling useful for debugging)
:param singularityRound: round singularities' coordinates to nearest integers
:param roundVertices: additional vertices that should be snapped to integer coordinates
:param hardFeatures: pairs of vertices that belong to edges that should be snapped to integer coordinates
:returns: UV, FUV where UV is #UV by 2 list of vertices in 2D and FUV is #FUV by 3 list of face indices in UV)");

    // Frame Field
    m.def(
        "frame_field",
        [](const MatrixXd& V,
           const MatrixXi& F,
           const VectorXi& b,
           const MatrixXd& bc1,
           const MatrixXd& bc2) {
            MatrixXd FF1, FF2;
            igl::copyleft::comiso::frame_field(V, F, b, bc1, bc2, FF1, FF2);
            return std::make_tuple(FF1, FF2);
        },
        "V"_a,
        "F"_a,
        "b"_a,
        "bc1"_a,
        "bc2"_a,
        R"(Generate a piecewise-constant frame-field from a sparse set of constraints on faces.

Based on "Frame Fields: Anisotropic and Non-Orthogonal Cross Fields"
by Daniele Panozzo, Enrico Puppo, Marco Tarini, Olga Sorkine-Hornung,
ACM Transactions on Graphics (SIGGRAPH, 2014)

:param V: #V by 3 list of mesh vertex coordinates
:param F: #F by 3 list of mesh faces (must be triangles)
:param b: #B by 1 list of constrained face indices
:param bc1: #B by 3 list of the constrained first representative vector of the frame field
:param bc2: #B by 3 list of the constrained second representative vector of the frame field
:returns: FF1, FF2 where FF1 is #F by 3 the first representative vector of the frame field and FF2 is #F by 3 the second representative vector of the frame field)");
}
