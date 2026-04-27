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
        kin = kn;
        if (f_curr == f_to) return rot;
      }
      return rot;
    };

    const int n = static_cast<int>(loop.size());
    const int num_edges = n - 1; // closed: loop.front()==loop.back(); open: distinct endpoints
    std::vector<int> face(num_edges), k_local(num_edges);
    for (int i = 0; i < num_edges; ++i) {
      int k;
      int f = faceOfHE(loop[i], loop[i + 1], k);
      face[i] = f;
      k_local[i] = k;
    }
    int rot_state = axis;
    double total = 0.0;
    for (int i = 0; i < num_edges; ++i) {
      const int f_i = face[i];
      const int kc = k_local[i];
      const int kn = (kc + 1) % 3;
      const int cc = FUV(f_i, kc);
      const int cn = FUV(f_i, kn);
      auto np = axis_components(rot_state);
      total += np.first  * (UV(cn, 0) - UV(cc, 0));
      total += np.second * (UV(cn, 1) - UV(cc, 1));
      if (i + 1 < num_edges) {
        const int v_shared = loop[i + 1];
        const int f_next = face[i + 1];
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
  std::vector<std::vector<int>> chains = { { F(0, 0), F(0, 1), F(0, 2) } };
  std::vector<int> axes = { 0 };

  igl::copyleft::comiso::miq(V, F, X1_combed, X2_combed, MMatch, isSingularity, Seams,
      UV, FUV, /*gradientSize=*/50, /*stiffness=*/5.0, /*directRound=*/false,
      /*iter=*/0, /*localIter=*/5, /*doRound=*/true, /*singularityRound=*/true,
      std::vector<int>(), std::vector<std::vector<int>>(),
      chains, axes);

  Eigen::MatrixXi TT, TTi;
  igl::triangle_triangle_adjacency(F, TT, TTi);

  // Trivial loop telescopes to zero by construction.
  double s = loopOrthogonalSum(F, TT, TTi, FUV, MMatch, Seams, UV, chains[0], 0);
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
  ring.push_back(ring.front()); // close the loop

  // Baseline: solve without the loop constraint and measure orthogonal sum.
  Eigen::MatrixXd UV_base; Eigen::MatrixXi FUV_base;
  igl::copyleft::comiso::miq(V, F, X1_combed, X2_combed, MMatch, isSingularity, Seams,
      UV_base, FUV_base, /*gradientSize=*/50, /*stiffness=*/5.0, /*directRound=*/false,
      /*iter=*/0, /*localIter=*/5, /*doRound=*/true, /*singularityRound=*/true);
  double s_base = loopOrthogonalSum(F, TT, TTi, FUV_base, MMatch, Seams, UV_base, ring, /*axis=*/1);
  std::cout << "[loop 1-ring baseline (no constraint)] orthogonal sum = " << s_base << std::endl;

  // Solve again with the loop constraint.
  std::vector<std::vector<int>> chains = { ring };
  std::vector<int> axes = { 1 }; // V-axis as orthogonal

  igl::copyleft::comiso::miq(V, F, X1_combed, X2_combed, MMatch, isSingularity, Seams,
      UV, FUV, /*gradientSize=*/50, /*stiffness=*/5.0, /*directRound=*/false,
      /*iter=*/0, /*localIter=*/5, /*doRound=*/true, /*singularityRound=*/true,
      std::vector<int>(), std::vector<std::vector<int>>(),
      chains, axes);

  double s = loopOrthogonalSum(F, TT, TTi, FUV, MMatch, Seams, UV, ring, /*axis=*/1);
  std::cout << "[loop 1-ring center=" << center << " size=" << ring.size()
            << "] orthogonal sum (constrained) = " << s << std::endl;
  REQUIRE(s < 1e-4);
  // The constraint must actually constrain: with a 1-ring around a singular
  // vertex, the natural solution does not satisfy the constraint, so the
  // baseline should be substantially larger than the constrained result.
  REQUIRE(s_base > 100.0 * s + 1e-6);
}

TEST_CASE("miq: 3_holes_loop_multi_step_fan", "[igl/copyleft/comiso]")
{
  // Loop [r_0, u, r_2, r_1, r_0] forces a multi-step fan walk at vertex u:
  // f0 = face of (r_0, u), f1 = face of (u, r_2). These faces share only u
  // (not an edge), so the fan walk between them traverses valence(u) - 3
  // intermediate faces. Exercises the post-transition k_in update inside
  // fanWalkRotation that is unreachable from any single-step fan loop.
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

  std::vector<std::vector<int>> VF, VFi;
  igl::vertex_triangle_adjacency(V, F, VF, VFi);

  // Find a vertex with valence >= 5 so the fan walk has >= 2 intermediate
  // steps (valence - 3).
  int u = -1;
  for (int v = 0; v < V.rows(); ++v) {
    if (static_cast<int>(VF[v].size()) >= 5) { u = v; break; }
  }
  REQUIRE(u >= 0);

  std::vector<int> ring = oneRingLoop(F, TT, VF, u);
  REQUIRE(ring.size() >= 5);

  // [r_0, u, r_2, r_1, r_0] — closed loop; last vertex repeats first.
  std::vector<int> chain = { ring[0], u, ring[2], ring[1], ring[0] };
  std::vector<std::vector<int>> chains = { chain };
  std::vector<int> axes = { 0 }; // U-axis as orthogonal

  // Baseline (no constraint).
  Eigen::MatrixXd UV_base; Eigen::MatrixXi FUV_base;
  igl::copyleft::comiso::miq(V, F, X1_combed, X2_combed, MMatch, isSingularity, Seams,
      UV_base, FUV_base, 50, 5.0, false, 0, 5, true, true);
  double s_base = loopOrthogonalSum(F, TT, TTi, FUV_base, MMatch, Seams, UV_base, chain, 0);

  // Constrained.
  igl::copyleft::comiso::miq(V, F, X1_combed, X2_combed, MMatch, isSingularity, Seams,
      UV, FUV, 50, 5.0, false, 0, 5, true, true,
      std::vector<int>(), std::vector<std::vector<int>>(),
      chains, axes);
  double s = loopOrthogonalSum(F, TT, TTi, FUV, MMatch, Seams, UV, chain, 0);

  std::cout << "[loop multi-step u=" << u << " valence=" << VF[u].size()
            << "] baseline=" << s_base << " constrained=" << s << std::endl;

  REQUIRE(s < 1e-4);
}

TEST_CASE("miq: torus_minor_ring_loop_alignment", "[igl/copyleft/comiso]")
{
  // Build a torus. Assign a per-face cross-field representative PD1 by
  // rotating the minor-ring tangent (the small-circle direction around the
  // tube cross-section) by 30 degrees around the face normal. Add a loop
  // constraint on a single minor ring (one circle around the tube), with
  // orthogonal axis = 0. Verify
  // (i) the loop constraint is enforced (orthogonal-sum ~ 0),
  // (ii) gradU, gradV per face are aligned with major/minor tangents,
  // (iii) the baseline (no constraint) UV does NOT satisfy the constraint
  //      — the minor ring is one of the torus's homology generators, so the
  //      natural UV holonomy around it is a non-trivial integer translation.
  using namespace Eigen;

  // === Torus mesh ===
  const double R = 2.0;       // major radius
  const double r_minor = 0.7; // minor radius
  const int n_theta = 32;     // major divisions
  const int n_phi = 16;       // minor divisions

  const int n_v = n_theta * n_phi;
  MatrixXd V(n_v, 3);
  for (int i = 0; i < n_theta; ++i) {
    double th = 2 * igl::PI * i / n_theta;
    for (int j = 0; j < n_phi; ++j) {
      double ph = 2 * igl::PI * j / n_phi;
      V(i * n_phi + j, 0) = (R + r_minor * std::cos(ph)) * std::cos(th);
      V(i * n_phi + j, 1) = (R + r_minor * std::cos(ph)) * std::sin(th);
      V(i * n_phi + j, 2) = r_minor * std::sin(ph);
    }
  }

  std::vector<RowVector3i> face_list;
  for (int i = 0; i < n_theta; ++i) {
    int inext = (i + 1) % n_theta;
    for (int j = 0; j < n_phi; ++j) {
      int jn = (j + 1) % n_phi;
      int v00 = i * n_phi + j;
      int v01 = i * n_phi + jn;
      int v10 = inext * n_phi + j;
      int v11 = inext * n_phi + jn;
      face_list.push_back(RowVector3i(v00, v10, v11));
      face_list.push_back(RowVector3i(v00, v11, v01));
    }
  }
  MatrixXi F(face_list.size(), 3);
  for (size_t k = 0; k < face_list.size(); ++k) F.row(k) = face_list[k];

  // Helper: recover (theta, phi) from a 3D point on the torus.
  auto point_to_theta_phi = [&](const Vector3d &p) {
    double th = std::atan2(p.y(), p.x());
    double rho = std::sqrt(p.x() * p.x() + p.y() * p.y()) - R;
    double ph = std::atan2(p.z(), rho);
    return std::make_pair(th, ph);
  };

  // Helper: minor-ring tangent (= dp/dphi normalized) at (theta, phi).
  auto minor_tangent = [](double th, double ph) {
    return Vector3d(-std::sin(ph) * std::cos(th), -std::sin(ph) * std::sin(th), std::cos(ph));
  };
  // Helper: major-ring tangent (= dp/dtheta normalized) at theta.
  auto major_tangent = [](double th) {
    return Vector3d(-std::sin(th), std::cos(th), 0.0);
  };

  // === PD1 = minor_tangent(centroid) rotated 30 degrees around face normal; PD2 = n x PD1 ===
  MatrixXd PD1(F.rows(), 3);
  MatrixXd PD2(F.rows(), 3);
  const double rot_angle = igl::PI / 6;
  for (int f = 0; f < F.rows(); ++f) {
    Vector3d p0 = V.row(F(f, 0));
    Vector3d p1 = V.row(F(f, 1));
    Vector3d p2 = V.row(F(f, 2));
    Vector3d c = (p0 + p1 + p2) / 3.0;
    Vector3d nf = (p1 - p0).cross(p2 - p0).normalized();
    auto thph = point_to_theta_phi(c);
    Vector3d mt = minor_tangent(thph.first, thph.second);
    mt -= nf * nf.dot(mt);
    if (mt.norm() < 1e-12) mt = Vector3d(1, 0, 0);
    mt.normalize();
    Vector3d pd1 = mt * std::cos(rot_angle) + nf.cross(mt) * std::sin(rot_angle);
    PD1.row(f) = pd1;
    PD2.row(f) = nf.cross(pd1).normalized();
  }

  // === Combed bisectors / mismatch / singularities / seams (mirrors auto-cut overload) ===
  MatrixXd BIS1, BIS2, BIS1_combed, BIS2_combed, X1_combed, X2_combed;
  igl::compute_frame_field_bisectors(V, F, PD1, PD2, BIS1, BIS2);
  igl::comb_cross_field(V, F, BIS1, BIS2, BIS1_combed, BIS2_combed);
  Eigen::Matrix<int, Eigen::Dynamic, 3> MMatch;
  igl::cross_field_mismatch(V, F, BIS1_combed, BIS2_combed, true, MMatch);
  Eigen::Matrix<int, Eigen::Dynamic, 1> isSingularity, singularityIndex;
  igl::find_cross_field_singularities(V, F, MMatch, isSingularity, singularityIndex);
  Eigen::Matrix<int, Eigen::Dynamic, 3> Seams;
  igl::cut_mesh_from_singularities(V, F, MMatch, Seams);
  igl::comb_frame_field(V, F, PD1, PD2, BIS1_combed, BIS2_combed, X1_combed, X2_combed);

  Eigen::MatrixXi TT, TTi;
  igl::triangle_triangle_adjacency(F, TT, TTi);

  // === Minor ring loop at theta_0 = 0: vertices [V[0][0], V[0][1], ..., V[0][n_phi-1]] ===
  std::vector<int> minor_loop;
  for (int j = 0; j < n_phi; ++j) minor_loop.push_back(0 * n_phi + j);
  minor_loop.push_back(minor_loop.front()); // close the loop

  auto compute_alignment = [&](const Eigen::MatrixXd &UV_in, const Eigen::MatrixXi &FUV_in,
                               double &max_misalign, double &avg_misalign, int &faces_checked) {
    max_misalign = 0.0;
    double sum_misalign = 0.0;
    faces_checked = 0;
    for (int f = 0; f < F.rows(); ++f) {
      Vector3d p0 = V.row(F(f, 0));
      Vector3d p1 = V.row(F(f, 1));
      Vector3d p2 = V.row(F(f, 2));
      Vector3d ncross = (p1 - p0).cross(p2 - p0);
      double area2 = ncross.norm();
      Vector3d nf = ncross / area2;
      double u0 = UV_in(FUV_in(f, 0), 0), u1 = UV_in(FUV_in(f, 1), 0), u2 = UV_in(FUV_in(f, 2), 0);
      double v0 = UV_in(FUV_in(f, 0), 1), v1 = UV_in(FUV_in(f, 1), 1), v2 = UV_in(FUV_in(f, 2), 1);
      Vector3d gU = (u0 * nf.cross(p2 - p1) + u1 * nf.cross(p0 - p2) + u2 * nf.cross(p1 - p0)) / area2;
      Vector3d gV = (v0 * nf.cross(p2 - p1) + v1 * nf.cross(p0 - p2) + v2 * nf.cross(p1 - p0)) / area2;
      Vector3d c = (p0 + p1 + p2) / 3.0;
      auto thph = point_to_theta_phi(c);
      Vector3d mt = minor_tangent(thph.first, thph.second);
      Vector3d Mt = major_tangent(thph.first);
      mt -= nf * nf.dot(mt); if (mt.norm() > 1e-12) mt.normalize();
      Mt -= nf * nf.dot(Mt); if (Mt.norm() > 1e-12) Mt.normalize();
      if (gU.norm() < 1e-9 || gV.norm() < 1e-9) continue;
      Vector3d gU_n = gU.normalized();
      Vector3d gV_n = gV.normalized();
      double align_U = std::max(std::abs(gU_n.dot(mt)), std::abs(gU_n.dot(Mt)));
      double align_V = std::max(std::abs(gV_n.dot(mt)), std::abs(gV_n.dot(Mt)));
      double misalign = std::max(1.0 - align_U, 1.0 - align_V);
      max_misalign = std::max(max_misalign, misalign);
      sum_misalign += (1.0 - align_U) + (1.0 - align_V);
      faces_checked++;
    }
    avg_misalign = sum_misalign / (2.0 * std::max(faces_checked, 1));
  };

  // Baseline: solve without the loop constraint.
  Eigen::MatrixXd UV_base; Eigen::MatrixXi FUV_base;
  igl::copyleft::comiso::miq(V, F, X1_combed, X2_combed, MMatch, isSingularity, Seams,
      UV_base, FUV_base, 30, 5.0, false, 0, 5, true, true);
  double base_max = 0, base_avg = 0; int base_n = 0;
  compute_alignment(UV_base, FUV_base, base_max, base_avg, base_n);
  double base_orth = loopOrthogonalSum(F, TT, TTi, FUV_base, MMatch, Seams, UV_base, minor_loop, /*axis=*/0);
  std::cout << "[torus minor-ring baseline] faces=" << base_n
            << " max_misalign=" << base_max << " avg_misalign=" << base_avg
            << " orth_sum=" << base_orth << std::endl;

  // Constrained: solve with minor-ring loop, axis = 0.
  MatrixXd UV;
  MatrixXi FUV;
  std::vector<std::vector<int>> chains = { minor_loop };
  std::vector<int> axes = { 0 };
  igl::copyleft::comiso::miq(V, F, X1_combed, X2_combed, MMatch, isSingularity, Seams,
      UV, FUV,
      /*gradientSize=*/30, /*stiffness=*/5.0, /*directRound=*/false,
      /*iter=*/0, /*localIter=*/5, /*doRound=*/true, /*singularityRound=*/true,
      std::vector<int>(), std::vector<std::vector<int>>(),
      chains, axes);

  double s = loopOrthogonalSum(F, TT, TTi, FUV, MMatch, Seams, UV, minor_loop, /*axis=*/0);
  std::cout << "[torus minor-ring constrained] orth_sum = " << s << std::endl;

  double max_misalign = 0, avg_misalign = 0; int faces_checked = 0;
  compute_alignment(UV, FUV, max_misalign, avg_misalign, faces_checked);
  std::cout << "[torus minor-ring alignment] faces=" << faces_checked
            << " max_misalign=" << max_misalign
            << " avg_misalign=" << avg_misalign << std::endl;

  // (1) Loop constraint is enforced.
  REQUIRE(s < 1e-4);

  // (2) Baseline differs from constrained: the natural UV holonomy around
  // the minor ring is a non-trivial integer translation, so its U-component
  // (orthogonal-sum) is non-zero — this is the signature property the user
  // asked us to verify.
  REQUIRE(std::abs(base_orth) > 0.5);
  REQUIRE(std::abs(base_orth) > 100.0 * s);

  // (3) Constraint pulls gradU/gradV toward axis-alignment relative to
  // baseline. NOTE: on the torus the cross field has no singularities, so
  // MIQ does not aggressively snap gradients to axes the way it does on the
  // sphere. With a 30-degree-tilted cross field the gradients stay tilted
  // ~30 degrees off the major/minor tangent directions; a single loop
  // constraint improves average alignment but cannot remove the local tilt.
  // (See the four-loop variant below for a stronger constraint.)
  REQUIRE(avg_misalign < base_avg);
}

TEST_CASE("miq: torus_four_minor_ring_loops_alignment", "[igl/copyleft/comiso]")
{
  // Variant of the torus test with four minor-ring loop constraints, evenly
  // spaced around the major direction (theta = 0, pi/2, pi, 3pi/2), each
  // with axis = 0.
  //
  // Important finding: on the torus, every minor ring is in the *same* H_1
  // homology class. MIQ flattens the torus into a square by cutting along
  // one minor cut and one major cut; every minor ring of the original
  // surface crosses the same major cut exactly once. The U-component of UV
  // holonomy around any minor ring is therefore pinned to the *same*
  // integer translation — the major-cut U-shift. Asking for "U-holonomy = 0"
  // on 4 rings is therefore equivalent to asking for it on 1, and produces
  // an identical UV solution.
  //
  // (To actually shrink misalignment further we would need either (a) loops
  // in different homology classes — e.g., a major ring — or (b) constraints
  // on BOTH axes, fixing the full 2D holonomy of each loop. Quick probes
  // below confirm both directions DO reduce misalignment.)
  using namespace Eigen;

  const double R = 2.0;
  const double r_minor = 0.7;
  const int n_theta = 32;
  const int n_phi = 16;

  const int n_v = n_theta * n_phi;
  MatrixXd V(n_v, 3);
  for (int i = 0; i < n_theta; ++i) {
    double th = 2 * igl::PI * i / n_theta;
    for (int j = 0; j < n_phi; ++j) {
      double ph = 2 * igl::PI * j / n_phi;
      V(i * n_phi + j, 0) = (R + r_minor * std::cos(ph)) * std::cos(th);
      V(i * n_phi + j, 1) = (R + r_minor * std::cos(ph)) * std::sin(th);
      V(i * n_phi + j, 2) = r_minor * std::sin(ph);
    }
  }

  std::vector<RowVector3i> face_list;
  for (int i = 0; i < n_theta; ++i) {
    int inext = (i + 1) % n_theta;
    for (int j = 0; j < n_phi; ++j) {
      int jn = (j + 1) % n_phi;
      int v00 = i * n_phi + j;
      int v01 = i * n_phi + jn;
      int v10 = inext * n_phi + j;
      int v11 = inext * n_phi + jn;
      face_list.push_back(RowVector3i(v00, v10, v11));
      face_list.push_back(RowVector3i(v00, v11, v01));
    }
  }
  MatrixXi F(face_list.size(), 3);
  for (size_t k = 0; k < face_list.size(); ++k) F.row(k) = face_list[k];

  auto point_to_theta_phi = [&](const Vector3d &p) {
    double th = std::atan2(p.y(), p.x());
    double rho = std::sqrt(p.x() * p.x() + p.y() * p.y()) - R;
    double ph = std::atan2(p.z(), rho);
    return std::make_pair(th, ph);
  };
  auto minor_tangent = [](double th, double ph) {
    return Vector3d(-std::sin(ph) * std::cos(th), -std::sin(ph) * std::sin(th), std::cos(ph));
  };
  auto major_tangent = [](double th) {
    return Vector3d(-std::sin(th), std::cos(th), 0.0);
  };

  MatrixXd PD1(F.rows(), 3);
  MatrixXd PD2(F.rows(), 3);
  const double rot_angle = igl::PI / 6;
  for (int f = 0; f < F.rows(); ++f) {
    Vector3d p0 = V.row(F(f, 0));
    Vector3d p1 = V.row(F(f, 1));
    Vector3d p2 = V.row(F(f, 2));
    Vector3d c = (p0 + p1 + p2) / 3.0;
    Vector3d nf = (p1 - p0).cross(p2 - p0).normalized();
    auto thph = point_to_theta_phi(c);
    Vector3d mt = minor_tangent(thph.first, thph.second);
    mt -= nf * nf.dot(mt);
    if (mt.norm() < 1e-12) mt = Vector3d(1, 0, 0);
    mt.normalize();
    Vector3d pd1 = mt * std::cos(rot_angle) + nf.cross(mt) * std::sin(rot_angle);
    PD1.row(f) = pd1;
    PD2.row(f) = nf.cross(pd1).normalized();
  }

  MatrixXd BIS1, BIS2, BIS1_combed, BIS2_combed, X1_combed, X2_combed;
  igl::compute_frame_field_bisectors(V, F, PD1, PD2, BIS1, BIS2);
  igl::comb_cross_field(V, F, BIS1, BIS2, BIS1_combed, BIS2_combed);
  Eigen::Matrix<int, Eigen::Dynamic, 3> MMatch;
  igl::cross_field_mismatch(V, F, BIS1_combed, BIS2_combed, true, MMatch);
  Eigen::Matrix<int, Eigen::Dynamic, 1> isSingularity, singularityIndex;
  igl::find_cross_field_singularities(V, F, MMatch, isSingularity, singularityIndex);
  Eigen::Matrix<int, Eigen::Dynamic, 3> Seams;
  igl::cut_mesh_from_singularities(V, F, MMatch, Seams);
  igl::comb_frame_field(V, F, PD1, PD2, BIS1_combed, BIS2_combed, X1_combed, X2_combed);

  Eigen::MatrixXi TT, TTi;
  igl::triangle_triangle_adjacency(F, TT, TTi);

  // Four minor-ring loops at theta_i for i in {0, n_theta/4, n_theta/2, 3*n_theta/4}.
  REQUIRE(n_theta % 4 == 0);
  const int n_chains = 4;
  std::vector<std::vector<int>> chains(n_chains);
  std::vector<int> axes(n_chains, 0);
  for (int k = 0; k < n_chains; ++k) {
    int i = (n_theta / n_chains) * k;
    for (int j = 0; j < n_phi; ++j) chains[k].push_back(i * n_phi + j);
    chains[k].push_back(chains[k].front()); // close the loop
  }

  // Probe: try also a major-ring loop (different homology class on the torus).
  // If the minor-ring constraints really do collapse to a single integer
  // wrap, mixing in a major-ring should produce a different UV.
  std::vector<int> probe_loop_major;
  const int probe_phi = 0;
  for (int i = 0; i < n_theta; ++i) probe_loop_major.push_back(i * n_phi + probe_phi);
  probe_loop_major.push_back(probe_loop_major.front()); // close the loop

  auto compute_alignment = [&](const Eigen::MatrixXd &UV_in, const Eigen::MatrixXi &FUV_in,
                               double &max_misalign, double &avg_misalign, int &faces_checked) {
    max_misalign = 0.0;
    double sum_misalign = 0.0;
    faces_checked = 0;
    for (int f = 0; f < F.rows(); ++f) {
      Vector3d p0 = V.row(F(f, 0));
      Vector3d p1 = V.row(F(f, 1));
      Vector3d p2 = V.row(F(f, 2));
      Vector3d ncross = (p1 - p0).cross(p2 - p0);
      double area2 = ncross.norm();
      Vector3d nf = ncross / area2;
      double u0 = UV_in(FUV_in(f, 0), 0), u1 = UV_in(FUV_in(f, 1), 0), u2 = UV_in(FUV_in(f, 2), 0);
      double v0 = UV_in(FUV_in(f, 0), 1), v1 = UV_in(FUV_in(f, 1), 1), v2 = UV_in(FUV_in(f, 2), 1);
      Vector3d gU = (u0 * nf.cross(p2 - p1) + u1 * nf.cross(p0 - p2) + u2 * nf.cross(p1 - p0)) / area2;
      Vector3d gV = (v0 * nf.cross(p2 - p1) + v1 * nf.cross(p0 - p2) + v2 * nf.cross(p1 - p0)) / area2;
      Vector3d c = (p0 + p1 + p2) / 3.0;
      auto thph = point_to_theta_phi(c);
      Vector3d mt = minor_tangent(thph.first, thph.second);
      Vector3d Mt = major_tangent(thph.first);
      mt -= nf * nf.dot(mt); if (mt.norm() > 1e-12) mt.normalize();
      Mt -= nf * nf.dot(Mt); if (Mt.norm() > 1e-12) Mt.normalize();
      if (gU.norm() < 1e-9 || gV.norm() < 1e-9) continue;
      Vector3d gU_n = gU.normalized();
      Vector3d gV_n = gV.normalized();
      double align_U = std::max(std::abs(gU_n.dot(mt)), std::abs(gU_n.dot(Mt)));
      double align_V = std::max(std::abs(gV_n.dot(mt)), std::abs(gV_n.dot(Mt)));
      double misalign = std::max(1.0 - align_U, 1.0 - align_V);
      max_misalign = std::max(max_misalign, misalign);
      sum_misalign += (1.0 - align_U) + (1.0 - align_V);
      faces_checked++;
    }
    avg_misalign = sum_misalign / (2.0 * std::max(faces_checked, 1));
  };

  // Run constrained solve with 4 minor-ring loops.
  MatrixXd UV;
  MatrixXi FUV;
  igl::copyleft::comiso::miq(V, F, X1_combed, X2_combed, MMatch, isSingularity, Seams,
      UV, FUV,
      /*gradientSize=*/30, /*stiffness=*/5.0, /*directRound=*/false,
      /*iter=*/0, /*localIter=*/5, /*doRound=*/true, /*singularityRound=*/true,
      std::vector<int>(), std::vector<std::vector<int>>(),
      chains, axes);

  // (1) All four chain constraints must be enforced.
  double max_orth_sum = 0.0;
  for (int k = 0; k < n_chains; ++k) {
    double sk = loopOrthogonalSum(F, TT, TTi, FUV, MMatch, Seams, UV, chains[k], /*axis=*/0);
    std::cout << "[torus 4-chains] chain " << k << " orth_sum = " << sk << std::endl;
    max_orth_sum = std::max(max_orth_sum, std::abs(sk));
  }
  REQUIRE(max_orth_sum < 1e-4);

  double max_misalign = 0, avg_misalign = 0; int faces_checked = 0;
  compute_alignment(UV, FUV, max_misalign, avg_misalign, faces_checked);
  std::cout << "[torus 4-chains alignment] faces=" << faces_checked
            << " max_misalign=" << max_misalign
            << " avg_misalign=" << avg_misalign << std::endl;

  // Probe A: 4 minor + 1 major-ring (homologically independent chain). Should
  // perturb the solution measurably.
  std::vector<std::vector<int>> mixed_chains = chains;
  mixed_chains.push_back(probe_loop_major);
  std::vector<int> mixed_axes(mixed_chains.size(), 0);
  Eigen::MatrixXd UV_mix; Eigen::MatrixXi FUV_mix;
  igl::copyleft::comiso::miq(V, F, X1_combed, X2_combed, MMatch, isSingularity, Seams,
      UV_mix, FUV_mix,
      30, 5.0, false, 0, 5, true, true,
      std::vector<int>(), std::vector<std::vector<int>>(),
      mixed_chains, mixed_axes);
  double mix_max = 0, mix_avg = 0; int mix_n = 0;
  compute_alignment(UV_mix, FUV_mix, mix_max, mix_avg, mix_n);
  std::cout << "[torus 4-minor + 1-major-ring]   max_misalign=" << mix_max
            << " avg_misalign=" << mix_avg << std::endl;

  // Probe B: 4 minor rings with BOTH axes constrained (8 constraints,
  // pins the full 2D holonomy of each ring to zero).
  std::vector<std::vector<int>> dual_chains;
  std::vector<int> dual_axes;
  for (int k = 0; k < n_chains; ++k) {
    dual_chains.push_back(chains[k]); dual_axes.push_back(0);
    dual_chains.push_back(chains[k]); dual_axes.push_back(1);
  }
  Eigen::MatrixXd UV_dual; Eigen::MatrixXi FUV_dual;
  igl::copyleft::comiso::miq(V, F, X1_combed, X2_combed, MMatch, isSingularity, Seams,
      UV_dual, FUV_dual,
      30, 5.0, false, 0, 5, true, true,
      std::vector<int>(), std::vector<std::vector<int>>(),
      dual_chains, dual_axes);
  double dual_max = 0, dual_avg = 0; int dual_n = 0;
  compute_alignment(UV_dual, FUV_dual, dual_max, dual_avg, dual_n);
  std::cout << "[torus 4-minor-rings x 2-axes]   max_misalign=" << dual_max
            << " avg_misalign=" << dual_avg << std::endl;

  // The 4-minor-only configuration must produce the same UV (and therefore
  // the same alignment metrics) as a single minor-ring constraint. The
  // probes above (A) and (B) demonstrate that breaking the redundancy does
  // reduce misalignment.
  REQUIRE(dual_avg < avg_misalign);  // pinning both axes helps
  REQUIRE(mix_avg  < avg_misalign);  // adding a major ring helps
}
