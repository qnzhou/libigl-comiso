#include <test_common.h>
#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/comb_cross_field.h>
#include <igl/comb_frame_field.h>
#include <igl/compute_frame_field_bisectors.h>
#include <igl/cross_field_mismatch.h>
#include <igl/cut_mesh_from_singularities.h>
#include <igl/find_cross_field_singularities.h>
#include <igl/local_basis.h>
#include <igl/readOFF.h>
#include <igl/rotate_vectors.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/copyleft/comiso/miq.h>
#include <igl/copyleft/comiso/nrosy.h>
#include <igl/PI.h>
#include <igl/serialize.h>
#include <sstream>
#include <iostream>
#include <cmath>
#include <igl/writeDMAT.h>

TEST_CASE("miq: 3_holes", "[igl/copyleft/comiso]")
{
using namespace Eigen;

// Input mesh
Eigen::MatrixXd V;
Eigen::MatrixXi F;

// Face barycenters
Eigen::MatrixXd B;

// Cross field
Eigen::MatrixXd X1,X2;

// Bisector field
Eigen::MatrixXd BIS1, BIS2;

// Combed bisector
Eigen::MatrixXd BIS1_combed, BIS2_combed;

// Per-corner, integer mismatches
Eigen::Matrix<int, Eigen::Dynamic, 3> MMatch;

// Field singularities
Eigen::Matrix<int, Eigen::Dynamic, 1> isSingularity, singularityIndex;

// Per corner seams
Eigen::Matrix<int, Eigen::Dynamic, 3> Seams;

// Combed field
Eigen::MatrixXd X1_combed, X2_combed;

// Global parametrization
Eigen::MatrixXd UV;
Eigen::MatrixXi FUV;

// Global parametrization (reference)
Eigen::MatrixXd UV_ref;
Eigen::MatrixXi FUV_ref;

// Load a mesh in OFF format
igl::readOFF(test_common::data_path("3holes.off"), V, F);

double gradient_size = 50;
double iter = 0;
double stiffness = 5.0;
bool direct_round = 0;

// Compute face barycenters
igl::barycenter(V, F, B);

// Contrain one face
VectorXi b(1);
b << 0;
MatrixXd bc(1, 3);
bc << 1, 0, 0;

// Create a smooth 4-RoSy field
VectorXd S;
igl::copyleft::comiso::nrosy(V, F, b, bc, VectorXi(), VectorXd(), MatrixXd(), 4, 0.5, X1, S);

// Find the orthogonal vector
MatrixXd B1, B2, B3;
igl::local_basis(V, F, B1, B2, B3);
X2 = igl::rotate_vectors(X1, VectorXd::Constant(1, igl::PI / 2), B1, B2);

// Always work on the bisectors, it is more general
igl::compute_frame_field_bisectors(V, F, X1, X2, BIS1, BIS2);

// Comb the field, implicitly defining the seams
igl::comb_cross_field(V, F, BIS1, BIS2, BIS1_combed, BIS2_combed);

// Find the integer mismatches
    igl::cross_field_mismatch(V, F, BIS1_combed, BIS2_combed, true, MMatch);

// Find the singularities
igl::find_cross_field_singularities(V, F, MMatch, isSingularity, singularityIndex);

// Cut the mesh, duplicating all vertices on the seams
igl::cut_mesh_from_singularities(V, F, MMatch, Seams);

// Comb the frame-field accordingly
igl::comb_frame_field(V, F, X1, X2, BIS1_combed, BIS2_combed, X1_combed, X2_combed);

// Global parametrization
igl::copyleft::comiso::miq(V,
          F,
          X1_combed,
          X2_combed,
          MMatch,
          isSingularity,
          Seams,
          UV,
          FUV,
          gradient_size,
          stiffness,
          direct_round,
          iter,
          5,
          true);

  // Refresh the test data
  // igl::writeDMAT(test_common::data_path("3holes-miq-UV.dmat"),UV);
  // igl::writeDMAT(test_common::data_path("3holes-miq-FUV.dmat"),FUV);

  igl::readDMAT(test_common::data_path("3holes-miq-UV.dmat"),UV_ref);
  igl::readDMAT(test_common::data_path("3holes-miq-FUV.dmat"),FUV_ref);

  REQUIRE (1e-6 > (UV-UV_ref).array().abs().maxCoeff());
  REQUIRE (1e-6 > (FUV-FUV_ref).array().abs().maxCoeff());
}

TEST_CASE("miq: 3_holes_quantization", "[igl/copyleft/comiso]")
{
  using namespace Eigen;

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::MatrixXd X1, X2;
  Eigen::MatrixXd BIS1, BIS2, BIS1_combed, BIS2_combed;
  Eigen::Matrix<int, Eigen::Dynamic, 3> MMatch;
  Eigen::Matrix<int, Eigen::Dynamic, 1> isSingularity, singularityIndex;
  Eigen::Matrix<int, Eigen::Dynamic, 3> Seams;
  Eigen::MatrixXd X1_combed, X2_combed;
  Eigen::MatrixXd UV;
  Eigen::MatrixXi FUV;

  igl::readOFF(test_common::data_path("3holes.off"), V, F);

  double gradient_size = 50;
  double iter = 0;
  double stiffness = 5.0;
  bool direct_round = 0;

  VectorXi b(1); b << 0;
  MatrixXd bc(1, 3); bc << 1, 0, 0;

  VectorXd S;
  igl::copyleft::comiso::nrosy(V, F, b, bc, VectorXi(), VectorXd(), MatrixXd(), 4, 0.5, X1, S);

  MatrixXd B1, B2, B3;
  igl::local_basis(V, F, B1, B2, B3);
  X2 = igl::rotate_vectors(X1, VectorXd::Constant(1, igl::PI / 2), B1, B2);

  igl::compute_frame_field_bisectors(V, F, X1, X2, BIS1, BIS2);
  igl::comb_cross_field(V, F, BIS1, BIS2, BIS1_combed, BIS2_combed);
  igl::cross_field_mismatch(V, F, BIS1_combed, BIS2_combed, true, MMatch);
  igl::find_cross_field_singularities(V, F, MMatch, isSingularity, singularityIndex);
  igl::cut_mesh_from_singularities(V, F, MMatch, Seams);
  igl::comb_frame_field(V, F, X1, X2, BIS1_combed, BIS2_combed, X1_combed, X2_combed);

  igl::copyleft::comiso::miq(V, F, X1_combed, X2_combed, MMatch, isSingularity, Seams,
      UV, FUV, gradient_size, stiffness, direct_round, iter, 5, true);

  // Triangle-triangle adjacency on the original (uncut) F.
  Eigen::MatrixXi TT, TTi;
  igl::triangle_triangle_adjacency(F, TT, TTi);

  // miq's seam constraint applies a 1<->3 swap on the mismatch when picking
  // the rotation. Mirror that here so we verify the same identity miq enforces.
  auto rotate_for_interval = [](int interval, const Eigen::Vector2d &p) {
    int swapped = interval;
    if (interval == 1)      swapped = 3;
    else if (interval == 3) swapped = 1;
    swapped = ((swapped % 4) + 4) % 4;
    Eigen::Matrix2d R;
    switch (swapped) {
      case 0: R << 1,0, 0,1;   break;
      case 1: R << 0,-1,1,0;   break; // +90
      case 2: R << -1,0, 0,-1; break;
      default:R << 0,1, -1,0;  break; // -90
    }
    return Eigen::Vector2d(R * p);
  };

  // (1) Across every seam edge, UV on side B should equal R(mismatch)*UV on
  // side A plus an integer translation.  Check max deviation from the nearest
  // integer translation across all seam edges.
  double max_seam_int_err = 0.0;
  int seam_edges_checked = 0;
  for (int f0 = 0; f0 < F.rows(); ++f0) {
    for (int k0 = 0; k0 < 3; ++k0) {
      if (Seams(f0, k0) == 0) continue;
      int f1 = TT(f0, k0);
      if (f1 < 0) continue;
      if (f1 < f0) continue; // process each edge once
      int k1 = TTi(f0, k0);

      // cut-mesh vertex indices for the two sides of the seam edge
      int v0a = FUV(f0, k0);
      int v1a = FUV(f0, (k0 + 1) % 3);
      int v1b = FUV(f1, k1);
      int v0b = FUV(f1, (k1 + 1) % 3);

      Eigen::Vector2d uv0a = UV.row(v0a);
      Eigen::Vector2d uv1a = UV.row(v1a);
      Eigen::Vector2d uv0b = UV.row(v0b);
      Eigen::Vector2d uv1b = UV.row(v1b);

      int interval = MMatch(f0, k0);
      Eigen::Vector2d t0 = uv0b - rotate_for_interval(interval, uv0a);
      Eigen::Vector2d t1 = uv1b - rotate_for_interval(interval, uv1a);

      // both translations should equal the same integer 2-vector
      for (int i = 0; i < 2; ++i) {
        max_seam_int_err = std::max(max_seam_int_err, std::abs(t0[i] - std::round(t0[i])));
        max_seam_int_err = std::max(max_seam_int_err, std::abs(t1[i] - std::round(t1[i])));
      }
      // and t0 should equal t1 (same integer per seam edge)
      max_seam_int_err = std::max(max_seam_int_err, (t0 - t1).cwiseAbs().maxCoeff());
      seam_edges_checked++;
    }
  }

  // (2) Singular vertices should have integer UV coordinates.
  double max_singularity_int_err = 0.0;
  int singular_checked = 0;
  for (int v = 0; v < V.rows(); ++v) {
    if (!isSingularity(v)) continue;
    // find an incident face and the corresponding cut-vertex via FUV
    int f_inc = -1, k_inc = -1;
    for (int f = 0; f < F.rows() && f_inc < 0; ++f) {
      for (int k = 0; k < 3; ++k) {
        if (F(f, k) == v) { f_inc = f; k_inc = k; break; }
      }
    }
    if (f_inc < 0) continue;
    int vc = FUV(f_inc, k_inc);
    Eigen::Vector2d uv = UV.row(vc);
    for (int i = 0; i < 2; ++i) {
      max_singularity_int_err = std::max(max_singularity_int_err, std::abs(uv[i] - std::round(uv[i])));
    }
    singular_checked++;
  }

  std::cout << "[miq quantization] seam edges checked: " << seam_edges_checked
            << ", max int-translation error: " << max_seam_int_err << std::endl;
  std::cout << "[miq quantization] singular vertices checked: " << singular_checked
            << ", max integer-coord error: " << max_singularity_int_err << std::endl;

  // Loose tolerance — solver is iterative.
  REQUIRE(max_seam_int_err < 1e-4);
  REQUIRE(max_singularity_int_err < 1e-4);
}

namespace
{
  // Post-hoc check of the loop constraint:
  // walk the loop, accumulate Σ_i n_i · delta_i with the same fan-walk
  // semantics as PoissonSolver::addLoopConstraint, and return the absolute
  // value. Should be ~0 for a constrained loop.
  static double loopOrthogonalSum(
      const Eigen::MatrixXi &F,
      const Eigen::MatrixXi &TT,
      const Eigen::MatrixXi &TTi,
      const Eigen::MatrixXi &FUV,
      const Eigen::Matrix<int, Eigen::Dynamic, 3> &mismatch,
      const Eigen::Matrix<int, Eigen::Dynamic, 3> &seams,
      const Eigen::MatrixXd &UV,
      const std::vector<int> &loop,
      int axis)
  {
    auto axis_components = [](int s) {
      int nx = 0, ny = 0;
      switch (((s % 4) + 4) % 4) {
        case 0: nx =  1; ny =  0; break;
        case 1: nx =  0; ny =  1; break;
        case 2: nx = -1; ny =  0; break;
        case 3: nx =  0; ny = -1; break;
      }
      return std::make_pair(nx, ny);
    };

    auto faceOfHE = [&](int v0, int v1, int &k_out) {
      for (int f = 0; f < F.rows(); ++f)
        for (int k = 0; k < 3; ++k)
          if (F(f, k) == v0 && F(f, (k + 1) % 3) == v1) { k_out = k; return f; }
      k_out = -1; return -1;
    };

    auto fanWalk = [&](int vertex, int f_from, int k_in, int f_to) {
      if (f_from == f_to) return 0;
      int rot = 0, f_curr = f_from, kin = k_in;
      for (int step = 0; step <= F.rows(); ++step) {
        const int k_v = (kin + 1) % 3;
        const int k_exit = k_v;
        if (seams(f_curr, k_exit) != 0) {
          int interval = mismatch(f_curr, k_exit);
          if (interval == 1) interval = 3;
          else if (interval == 3) interval = 1;
          rot = (rot + interval) % 4;
        }
        int fn = TT(f_curr, k_exit);
        if (fn < 0) return rot;
        int kn = TTi(f_curr, k_exit);
        f_curr = fn;
        kin = (kn + 2) % 3;
        if (f_curr == f_to) return rot;
      }
      return rot;
    };

    const int n = static_cast<int>(loop.size());
    std::vector<int> face(n), k_local(n);
    for (int i = 0; i < n; ++i) {
      int k;
      int f = faceOfHE(loop[i], loop[(i + 1) % n], k);
      face[i] = f;
      k_local[i] = k;
    }
    int rot_state = axis;
    double total = 0.0;
    for (int i = 0; i < n; ++i) {
      const int f_i = face[i];
      const int kc = k_local[i];
      const int kn = (kc + 1) % 3;
      const int cc = FUV(f_i, kc);
      const int cn = FUV(f_i, kn);
      auto np = axis_components(rot_state);
      total += np.first  * (UV(cn, 0) - UV(cc, 0));
      total += np.second * (UV(cn, 1) - UV(cc, 1));
      if (i + 1 < n) {
        const int v_shared = loop[(i + 1) % n];
        const int f_next = face[(i + 1) % n];
        if (f_i != f_next) {
          int dr = fanWalk(v_shared, f_i, kc, f_next);
          rot_state = (rot_state + dr) % 4;
        }
      }
    }
    return std::abs(total);
  }

  // Build the closed 1-ring boundary loop of vertex `center`, ordered so
  // that consecutive edges (v_i, v_{i+1}) are oriented half-edges of the
  // mesh (each lying in a face that also contains `center`).
  static std::vector<int> oneRingLoop(
      const Eigen::MatrixXi &F,
      const Eigen::MatrixXi &TT,
      const std::vector<std::vector<int>> &VF,
      int center)
  {
    if (VF[center].empty()) return {};
    int f_start = VF[center][0];
    int k_start = -1;
    for (int k = 0; k < 3; ++k) if (F(f_start, k) == center) { k_start = k; break; }
    // half-edge (F(f_start,(k_start+1)%3), F(f_start,(k_start+2)%3)) is the
    // 1-ring boundary edge in this face.
    std::vector<int> ring;
    int f = f_start, k = k_start;
    do {
      int u = F(f, (k + 1) % 3);
      ring.push_back(u);
      // next face around `center`: cross edge (center, F(f,(k+2)%3)),
      // i.e., local edge (k+2)%3.
      int kn_edge = (k + 2) % 3;
      int f_next = TT(f, kn_edge);
      if (f_next < 0) return {}; // boundary; abort
      // find local index of center in f_next
      int kc = -1;
      for (int kk = 0; kk < 3; ++kk) if (F(f_next, kk) == center) { kc = kk; break; }
      f = f_next;
      k = kc;
    } while (f != f_start);
    return ring;
  }
}

TEST_CASE("miq: 3_holes_loop_trivial", "[igl/copyleft/comiso]")
{
  using namespace Eigen;
  Eigen::MatrixXd V; Eigen::MatrixXi F;
  Eigen::MatrixXd X1, X2, BIS1, BIS2, BIS1_combed, BIS2_combed, X1_combed, X2_combed;
  Eigen::Matrix<int, Eigen::Dynamic, 3> MMatch;
  Eigen::Matrix<int, Eigen::Dynamic, 1> isSingularity, singularityIndex;
  Eigen::Matrix<int, Eigen::Dynamic, 3> Seams;
  Eigen::MatrixXd UV; Eigen::MatrixXi FUV;

  igl::readOFF(test_common::data_path("3holes.off"), V, F);

  VectorXi b(1); b << 0;
  MatrixXd bc(1, 3); bc << 1, 0, 0;
  VectorXd S;
  igl::copyleft::comiso::nrosy(V, F, b, bc, VectorXi(), VectorXd(), MatrixXd(), 4, 0.5, X1, S);
  MatrixXd B1, B2, B3;
  igl::local_basis(V, F, B1, B2, B3);
  X2 = igl::rotate_vectors(X1, VectorXd::Constant(1, igl::PI / 2), B1, B2);
  igl::compute_frame_field_bisectors(V, F, X1, X2, BIS1, BIS2);
  igl::comb_cross_field(V, F, BIS1, BIS2, BIS1_combed, BIS2_combed);
  igl::cross_field_mismatch(V, F, BIS1_combed, BIS2_combed, true, MMatch);
  igl::find_cross_field_singularities(V, F, MMatch, isSingularity, singularityIndex);
  igl::cut_mesh_from_singularities(V, F, MMatch, Seams);
  igl::comb_frame_field(V, F, X1, X2, BIS1_combed, BIS2_combed, X1_combed, X2_combed);

  // Trivial loop: the three vertices of face 0, walked in order.
  std::vector<std::vector<int>> loops = { { F(0, 0), F(0, 1), F(0, 2) } };
  std::vector<int> axes = { 0 };

  igl::copyleft::comiso::miq(V, F, X1_combed, X2_combed, MMatch, isSingularity, Seams,
      UV, FUV, /*gradientSize=*/50, /*stiffness=*/5.0, /*directRound=*/false,
      /*iter=*/0, /*localIter=*/5, /*doRound=*/true, /*singularityRound=*/true,
      std::vector<int>(), std::vector<std::vector<int>>(),
      loops, axes);

  Eigen::MatrixXi TT, TTi;
  igl::triangle_triangle_adjacency(F, TT, TTi);

  // Trivial loop telescopes to zero by construction.
  double s = loopOrthogonalSum(F, TT, TTi, FUV, MMatch, Seams, UV, loops[0], 0);
  std::cout << "[loop trivial] orthogonal sum = " << s << std::endl;
  REQUIRE(s < 1e-9);
}

TEST_CASE("miq: 3_holes_loop_one_ring", "[igl/copyleft/comiso]")
{
  using namespace Eigen;
  Eigen::MatrixXd V; Eigen::MatrixXi F;
  Eigen::MatrixXd X1, X2, BIS1, BIS2, BIS1_combed, BIS2_combed, X1_combed, X2_combed;
  Eigen::Matrix<int, Eigen::Dynamic, 3> MMatch;
  Eigen::Matrix<int, Eigen::Dynamic, 1> isSingularity, singularityIndex;
  Eigen::Matrix<int, Eigen::Dynamic, 3> Seams;
  Eigen::MatrixXd UV; Eigen::MatrixXi FUV;

  igl::readOFF(test_common::data_path("3holes.off"), V, F);

  VectorXi b(1); b << 0;
  MatrixXd bc(1, 3); bc << 1, 0, 0;
  VectorXd S;
  igl::copyleft::comiso::nrosy(V, F, b, bc, VectorXi(), VectorXd(), MatrixXd(), 4, 0.5, X1, S);
  MatrixXd B1, B2, B3;
  igl::local_basis(V, F, B1, B2, B3);
  X2 = igl::rotate_vectors(X1, VectorXd::Constant(1, igl::PI / 2), B1, B2);
  igl::compute_frame_field_bisectors(V, F, X1, X2, BIS1, BIS2);
  igl::comb_cross_field(V, F, BIS1, BIS2, BIS1_combed, BIS2_combed);
  igl::cross_field_mismatch(V, F, BIS1_combed, BIS2_combed, true, MMatch);
  igl::find_cross_field_singularities(V, F, MMatch, isSingularity, singularityIndex);
  igl::cut_mesh_from_singularities(V, F, MMatch, Seams);
  igl::comb_frame_field(V, F, X1, X2, BIS1_combed, BIS2_combed, X1_combed, X2_combed);

  Eigen::MatrixXi TT, TTi;
  igl::triangle_triangle_adjacency(F, TT, TTi);

  // Build vertex-triangle adjacency to construct a 1-ring loop.
  std::vector<std::vector<int>> VF, VFi;
  igl::vertex_triangle_adjacency(V, F, VF, VFi);

  // Pick a singular vertex — its 1-ring has non-trivial cross-field holonomy,
  // so the natural miq solution generally does NOT satisfy our constraint.
  int center = -1;
  for (int v = 0; v < V.rows(); ++v) {
    if (isSingularity(v) && static_cast<int>(VF[v].size()) >= 3) { center = v; break; }
  }
  REQUIRE(center >= 0);

  std::vector<int> ring = oneRingLoop(F, TT, VF, center);
  REQUIRE(ring.size() >= 4);

  // Baseline: solve without the loop constraint and measure orthogonal sum.
  Eigen::MatrixXd UV_base; Eigen::MatrixXi FUV_base;
  igl::copyleft::comiso::miq(V, F, X1_combed, X2_combed, MMatch, isSingularity, Seams,
      UV_base, FUV_base, /*gradientSize=*/50, /*stiffness=*/5.0, /*directRound=*/false,
      /*iter=*/0, /*localIter=*/5, /*doRound=*/true, /*singularityRound=*/true);
  double s_base = loopOrthogonalSum(F, TT, TTi, FUV_base, MMatch, Seams, UV_base, ring, /*axis=*/1);
  std::cout << "[loop 1-ring baseline (no constraint)] orthogonal sum = " << s_base << std::endl;

  // Solve again with the loop constraint.
  std::vector<std::vector<int>> loops = { ring };
  std::vector<int> axes = { 1 }; // V-axis as orthogonal

  igl::copyleft::comiso::miq(V, F, X1_combed, X2_combed, MMatch, isSingularity, Seams,
      UV, FUV, /*gradientSize=*/50, /*stiffness=*/5.0, /*directRound=*/false,
      /*iter=*/0, /*localIter=*/5, /*doRound=*/true, /*singularityRound=*/true,
      std::vector<int>(), std::vector<std::vector<int>>(),
      loops, axes);

  double s = loopOrthogonalSum(F, TT, TTi, FUV, MMatch, Seams, UV, ring, /*axis=*/1);
  std::cout << "[loop 1-ring center=" << center << " size=" << ring.size()
            << "] orthogonal sum (constrained) = " << s << std::endl;
  REQUIRE(s < 1e-4);
  // The constraint must actually constrain: with a 1-ring around a singular
  // vertex, the natural solution does not satisfy the constraint, so the
  // baseline should be substantially larger than the constrained result.
  REQUIRE(s_base > 100.0 * s + 1e-6);
}
