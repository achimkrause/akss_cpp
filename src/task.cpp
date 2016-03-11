#include "task.h"
#include "spectral_sequence.h"
#include "morphisms.h"

Task::Task(Session& session) : session_(session)
{
}

GroupTask::GroupTask(Session& session, const std::size_t p, const std::size_t q)
    : Task(session), p_(p), q_(q)
{
}

bool GroupTask::autosolve()
{
  SpectralSequence& sequence = session_.get_sequence();
  std::pair<std::size_t, std::size_t> bounds = sequence.get_bounds(q_);

  for (std::size_t s = bounds.first; s <= bounds.second; s++) {
    AbelianGroup null_q_s = sequence.get_e_2(TrigradedIndex(0, q_, s));

    int mon_rank;  // compute rank of Z[l_i] in degree p here!

    AbelianGroup result(null_q_s.free_rank() * mon_rank,
                        null_q_s.tor_rank() * mon_rank);
    for (int i = 0; i < null_q_s.tor_rank(); i++) {
      for (int j = 0; j < mon_rank; j++) {
        result(i * mon_rank + j) = null_q_s(i);
      }
    }
    sequence.set_e2(TrigradedIndex(p_, q_, s), result);
  }
  return true;
}

ExtensionTask::ExtensionTask(Session& session, int q, int s)
    : Task(session), q_(q), s_(s)
{
}

bool ExtensionTask::autosolve()
{
  SpectralSequence& sequence = session_.get_sequence();
  for (int n = 1; n <= q_; n++) {
    // take the term at (n,q-n+1,s-1) taking into account differentials up to
    // n-1 leaving,
    // and differentials up to q-n+2 entering.
    GroupWithMorphisms eab = sequence.get_e_ab(
        TrigradedIndex(n, q_ - n + 1, s_ - 1), n - 1, q_ - n + 2);
    if (eab.group.rank() != 0) {
      list_groups_.emplace(n, eab.group);
      list_maps_.emplace(n, eab.maps_to.front());
    }
  }

  // now, for n=q_+1, s=1, we have to compute e_ab mod the v_n.
  if (s_ == 1) {
    AbelianGroup iterated_kernel =
        sequence.get_kernel(TrigradedIndex(q_ + 1, 0, 0), q_ + 1);
    MatrixQ inclusion =
        sequence.get_inclusion(TrigradedIndex(q_ + 1, 0, 0), q_ + 1);

    MatrixQ v_i_map(0, 0);  // TODO:get from nat's table

    MatrixQ matrix = lift_from_free(
        sequence.get_prime(), v_i_map, inclusion,
        iterated_kernel);  // compute lift of v_i_map along inclusion.

    GroupWithMorphisms coker =
        compute_cokernel(sequence.get_prime(), matrix, iterated_kernel,
                         MatrixQRefList(), MatrixQRefList());
    list_groups_.insert(std::pair<int, AbelianGroup>(q_ + 1, coker.group));
  }

  if (list_groups_.size() == 0) {
    sequence.set_e2(TrigradedIndex(0, q_, s_), AbelianGroup(0, 0));
    for (int i = 2; i <= q_ + 1; i++) {
      sequence.set_diff_zero(TrigradedIndex(i, q_ - i + 1, s_ - 1), i);
    }
    return true;
  }
  if (list_groups_.size() == 1) {
    sequence.set_e2(TrigradedIndex(0, q_, s_), list_groups_.begin()->second);
    int r = list_groups_.begin()->first;

    // set the first r-1 differentials to 0
    for (int i = 2; i < r; i++) {
      sequence.set_diff_zero(TrigradedIndex(i, q_ - i + 1, s_ - 1), i);
    }

    // the r'th to the projection from
    // r-1'st kernel -> r-1'st kernel / image of all differentials
    sequence.set_diff(TrigradedIndex(r, q_ - r + 1, s_ - 1), r,
                      list_maps_.begin()->second);

    // and the rest 0.
    for (int i = r + 1; i <= q_ + 1; i++) {
      sequence.set_diff_zero(TrigradedIndex(i, q_ - i + 1, s_ - 1), i);
    }
    return true;
  }
  return false;
}

DifferentialTask::DifferentialTask(Session& session, TrigradedIndex index,
                                   int r)
    : Task(session), index_(index), r_(r)
{
}

bool DifferentialTask::autosolve()
{
  // First, lift the map d_r: (r,q,s) -> (0,q-r+1,s+1) to a map lift between
  // frees
  //(this just means lifting it against the projection E2 -> r-1'st cokernel),
  // and then taking the corresponding matrix)
  // For each sequence I in the appropriate degree, compute the induced map on
  // r-1'st kernels
  //(this is lifting against inclusion, we use nat's matrices for r_I), and then
  // view
  // lift\circ r_I as a map from an integral lift of the r-1_st kernel at
  // (p,q,s) to E2 at (0,..).
  // By filling in the entries of all these matrices correctly, we can build
  // a big matrix for the whole differential. This is a lift of the differential
  // to E2
  // Now defining K as the kernel of the projection E2 -> r-1st cokernel (or
  // rather, the free lift of E2),
  // we can compute the composition K\tensor Z{l_I} -> E2\tensor Z{l_I} -> r-1st
  // cokernel.
  // this is our indeterminacy.
  // The composition of our big matrix with the projection onto the r-1st
  // cokernel is our representative.
  SpectralSequence& sequence = session_.get_sequence();

  AbelianGroup e2_left_img =
      sequence.get_e_2(TrigradedIndex(0, index_.q() + r_ - 1, index_.s() + 1));
  AbelianGroup er_left_img = sequence.get_cokernel(
      TrigradedIndex(0, index_.q() + r_ - 1, index_.s() + 1), r_);
  MatrixQ projection_left_img = sequence.get_projection(
      TrigradedIndex(0, index_.q() + r_ + 1, index_.s() + 1), r_);

  MatrixQ diff_left =
      sequence.get_diff_from(TrigradedIndex(r_, index_.q(), index_.s()), r_);

  //"lift" the differential from (r,q,s) -> (0,q-r+1,s+1) over
  // projection_left_img
  //(actually just a free presentation)
  MatrixQ lift = lift_from_free(sequence.get_prime(), diff_left,
                                projection_left_img, er_left_img);

  AbelianGroup e2_0_q_s =
      sequence.get_e_2(TrigradedIndex(0, index_.q(), index_.s()));
  AbelianGroup ker_right_domain = sequence.get_kernel(index_, r_);

  std::size_t mon_rank = session_.get_monomial_rank(index_.p() - r_);
  MatrixQ result_lift(mon_rank * e2_left_img.rank(), ker_right_domain.rank());

  MatrixQ inclusion_right_domain = sequence.get_inclusion(index_, r_);
  AbelianGroup ker_left_domain =
      sequence.get_kernel(TrigradedIndex(r_, index_.q(), index_.s()), r_);
  MatrixQ inclusion_left_domain =
      sequence.get_inclusion(TrigradedIndex(r_, index_.q(), index_.s()), r_);

  for (int i = 0; i < mon_rank; i++) {
    MatrixQ r_I;  // obtain the i_th operation from degree p to degree r here
                  // from nat's tables.
    MatrixQ r_I_q(r_I.height() * e2_0_q_s.rank(),
                  r_I.width() * e2_0_q_s.rank());  // obtain the tensor product
                                                   // A\otimes r_I, where A is
                                                   // the group e2_0_q_s.
    for (int diag = 0; diag < e2_0_q_s.rank(); diag++) {
      for (int height = 0; height < r_I.height(); height++) {
        for (int width = 0; width < r_I.width(); width++) {
          r_I_q(diag * r_I.height() + height, diag * r_I.width() + width) =
              r_I(height, width);
        }
      }
    }

    // lift r_I_q * inclusion_right_domain against
    // inclusion_left_domain (okay because this is injective).
    MatrixQ r_I_ker =
        lift_from_free(sequence.get_prime(), r_I_q * inclusion_right_domain,
                       inclusion_left_domain, ker_left_domain);

    MatrixQ lift_r_I = lift * r_I_ker;
    for (int h = 0; h < e2_left_img.rank(); h++) {
      for (int w = 0; w < ker_right_domain.rank(); w++) {
        result_lift(h * mon_rank + i, w) = r_I_ker(h, w);
      }
    }
  }

  MatrixQ projection_right_img = sequence.get_projection(
      TrigradedIndex(index_.p() - r_, index_.q() - r_ + 1, index_.s() + 1), r_);

  diff_candidate_ = projection_right_img * result_lift;

  // now also determine indeterminacy:
  MatrixQ id = IdentityMatrix<mpq_class>(projection_right_img.width());
  MatrixQList from_X;
  from_X.emplace_back(id);
  GroupWithMorphisms ker_proj_morphisms =
      compute_kernel(sequence.get_prime(), projection_left_img, e2_left_img,
                     er_left_img, MatrixQRefList(), ref(from_X));
  MatrixQ from_K = *ker_proj_morphisms.maps_from.begin();
  MatrixQ from_K_tensor(from_K.height() * mon_rank, from_K.width() * mon_rank);
  for (int i = 0; i < from_K.height(); i++) {
    for (int j = 0; j < from_K.width(); j++) {
      for (int d = 0; d < mon_rank; d++) {
        from_K_tensor(i * mon_rank + d, j * mon_rank + d) = from_K(i, j);
      }
    }
  }

  AbelianGroup coker_right_img = sequence.get_cokernel(
      TrigradedIndex(index_.p() - r_, index_.q() + r_ - 1, index_.s() + 1), r_);
  MatrixQ indet_map = projection_right_img * from_K_tensor;
  indeterminacy_ = compute_image(
      sequence.get_prime(), indet_map, AbelianGroup(indet_map.width(), 0),
      coker_right_img);  // compute image of indet_map.
  if (indeterminacy_.group.rank() == 0) {
    sequence.set_diff(index_, r_, diff_candidate_);
    return true;
  }

  return false;
}
