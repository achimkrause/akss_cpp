#include <sstream>
#include "morphisms.h"
#include "p_local.h"
#include "spectral_sequence.h"
#include "task.h"

Task::Task(Session& session) : session_(session)
{
}

GroupTask::GroupTask(Session& session, const deg_t p, const deg_t q)
    : Task(session), p_(p), q_(q)
{
}

bool GroupTask::autosolve()
{
  SpectralSequence& sequence = session_.get_sequence();
  std::pair<deg_t, deg_t> bounds = sequence.get_bounds(q_);

  for (deg_t s = bounds.first; s <= bounds.second; s++) {
    AbelianGroup null_q_s = sequence.get_e_2(TrigradedIndex(0, q_, s));

    dim_t mon_rank = session_.get_monomial_rank(p_);

    AbelianGroup result(null_q_s.free_rank() * mon_rank,
                        null_q_s.tor_rank() * mon_rank);
    for (dim_t i = 0; i < null_q_s.tor_rank(); i++) {
      for (dim_t j = 0; j < mon_rank; j++) {
        result(i * mon_rank + j) = null_q_s(i);
      }
    }
    sequence.set_e2(TrigradedIndex(p_, q_, s), result);
  }
  return true;
}

bool GroupTask::usersolve() {
  throw std::logic_error("GroupTask::usersolve(): should not be called.");
}

void GroupTask::display_detail() {
  throw std::logic_error("GroupTask::display_task(): should not be called.");
}

void GroupTask::display_overview() {
  throw std::logic_error("GroupTask::display_overview(): should not be called.");
}

ExtensionTask::ExtensionTask(Session& session, deg_t q, deg_t s)
    : Task(session), q_(q), s_(s)
{
}

bool ExtensionTask::autosolve()
{
  SpectralSequence& sequence = session_.get_sequence();

  for (dim_t n = 2; n <= q_; n++) {
    // take the term at (n,q-n+1,s-1) taking into account differentials up to
    // n-1 leaving,
    // and differentials up to q-n+2 entering.

    GroupWithMorphisms eab = sequence.get_e_ab(
        source(TrigradedIndex(0, q_, s_), n), n, q_ - n + 3);
    if (eab.group.rank() != 0) {
      list_groups_.emplace(n, eab.group);
      list_maps_.emplace(n, eab.maps_to.front());
    }
  }
  // now, for n=q_+1, s=1, we have to compute e_ab mod the v_n.
  if (s_ == 1) {
    AbelianGroup iterated_kernel =
        sequence.get_kernel(TrigradedIndex(q_ + 1, 0, 0), q_+1);
    MatrixQ inclusion =
        sequence.get_inclusion(TrigradedIndex(q_ + 1, 0, 0), q_+1);

    MatrixQ v_i_map = session_.get_v_inclusion(q_ + 1);
    MatrixQ matrix = lift_from_free(
        sequence.get_prime(), v_i_map, inclusion,
        iterated_kernel);  // compute lift of v_i_map along inclusion.
    MatrixQ id = MatrixQ::identity(matrix.height());
    GroupWithMorphisms coker = compute_cokernel(
        sequence.get_prime(), matrix, iterated_kernel, {id}, MatrixQRefList());
    if (coker.group.rank() > 0) {
      list_groups_.insert(std::pair<deg_t, AbelianGroup>(q_ + 1, coker.group));
      list_maps_.insert(
          std::pair<deg_t, MatrixQ>(q_ + 1, coker.maps_to.front()));
    }
  }
  if (list_groups_.size() == 0) {
    sequence.set_e2(TrigradedIndex(0, q_, s_), AbelianGroup(0, 0));
    for (dim_t i = 2; i <= q_+1; i++) {
      sequence.set_diff_zero(source(TrigradedIndex(0, q_, s_), i), i);
    }
    return true;
  }
  if (list_groups_.size() == 1) {
    sequence.set_e2(TrigradedIndex(0, q_, s_), list_groups_.begin()->second);
    dim_t r = list_groups_.begin()->first;

    // set the first r-1 differentials to 0
    for (dim_t i = 2; i < r; i++) {
      sequence.set_diff_zero(source(TrigradedIndex(0, q_, s_), i), i);
    }

    // the r'th to the projection from
    // r-1'st kernel -> r-1'st kernel / image of all differentials
    sequence.set_diff(source(TrigradedIndex(0, q_, s_), r), r,
                      list_maps_.begin()->second);
    // and the rest 0.
    for (dim_t i = r + 1; i <= q_+1; i++) {
      sequence.set_diff_zero(source(TrigradedIndex(0, q_, s_), i), i);
    }
    return true;
  }
  return false;
}

bool ExtensionTask::usersolve() {
  std::cout << "Enter the E^2-group at (q,s)=(" << q_ << "," << s_ << "):\n";
  std::string str_grp;
  std::getline(std::cin,str_grp);
  std::stringstream input(str_grp);

  SpectralSequence& sequence = session_.get_sequence();

  AbelianGroup grp;
  parse_abelian_group(input, grp, sequence.get_prime());

  sequence.set_e2(TrigradedIndex(0,q_,s_),grp);
  auto group_it = list_groups_.begin();

  deg_t coker_r = 2;
  for(; group_it!=list_groups_.end(); group_it++){
    deg_t r = group_it->first;

    while(coker_r < r){
      sequence.set_diff_zero(source(TrigradedIndex(0, q_, s_),coker_r), coker_r);
      coker_r++;
    }
    AbelianGroup coker = sequence.get_cokernel(TrigradedIndex(0,q_,s_),r);
    std::stringstream filename;
    filename << "ext_task_" << q_ << "_" << s_ << "_" << group_it->first << ".dat";
    std::stringstream text;
    text << "Change the Matrix below to the injective d_"
         << r
         <<" from ";
    group_it->second.print(text, sequence.get_prime());
    text << " into ";
    coker.print(text, sequence.get_prime());
    text << "\n";

    session_.matrix_file_dialog(
            coker.rank(),
            group_it->second.rank(),
            filename.str(),
            text.str());
    MatrixQ diff = session_.read_matrix_file(coker.rank(), group_it->second.rank(),filename.str());
    TrigradedIndex pqs(r, q_+1-r, s_-1);
    MatrixQ proj = sequence.get_e_ab(pqs,r, q_-r + 3).maps_to[0];

    session_.get_sequence().set_diff(pqs, r,diff*proj);
    coker_r++;
    list_maps_.emplace(r, diff*proj);
  }
  return true;
}

void ExtensionTask::display_detail() {
  std::cout << "ExtensionTask for s=" << s_ << ", q=" << q_ << "\n";
  std::cout << "The groups are:\n";

  for (auto group_it = list_groups_.begin(); group_it != list_groups_.end();
       group_it++) {
    std::size_t p = group_it->first;
    std::cout << "   Degree (p,q,s) = (" << p << ","
              << static_cast<dim_t>(q_) + 1 - p << "," << s_ - 1 << ")"
              << ": ";
    group_it->second.print(std::cout, session_.get_sequence().get_prime());
    std::cout << "\n";
  }
}

void ExtensionTask::display_overview() {
  std::cout << "ExtensionTask for s="<<s_ <<", q="<<q_<<"\n";
}

DifferentialTask::DifferentialTask(Session& session, TrigradedIndex index,
                                   dim_t r)
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

  if (index_.p() % 2 == 1 || (index_.p() - static_cast<deg_t>(r_)) % 2 == 1) {
    sequence.set_diff_zero(index_, r_);
    return true;
  }

  AbelianGroup ker_right_domain = sequence.get_kernel(index_, r_);
  AbelianGroup coker_right_codomain =
      sequence.get_cokernel(target(index_, r_), r_);

  GroupWithMorphisms e_right_domain = sequence.get_e_ab(index_, r_, index_.q()+2);
  if (e_right_domain.group.rank() == 0 || coker_right_codomain.rank() == 0) {
    sequence.set_diff_zero(index_, r_);
    return true;
  }

  deg_t r_s = static_cast<deg_t>(r_);
  AbelianGroup e2_left_codomain =
      sequence.get_e_2(TrigradedIndex(0, index_.q() + r_s - 1, index_.s() + 1));
  //std::cout << "e2_left_codomain:\n";
  //e2_left_codomain.print(std::cout, sequence.get_prime());
  //std::cout << "\n";
  AbelianGroup er_left_codomain = sequence.get_cokernel(
      TrigradedIndex(0, index_.q() + r_s - 1, index_.s() + 1), r_);
  //std::cout << "er_left_codomain:\n";
  //er_left_codomain.print(std::cout, sequence.get_prime());
  //std::cout << "\n";
  MatrixQ projection_left_img = sequence.get_projection(
      TrigradedIndex(0, index_.q() + r_s - 1, index_.s() + 1), r_);
  //std::cout << "MatrixQ projection_left_img:\n" << projection_left_img << "\n";


  MatrixQ diff_left =
      sequence.get_diff_from(TrigradedIndex(r_s, index_.q(), index_.s()), r_);
  //std::cout << "MatrixQ diff_left:\n" << diff_left << "\n";
  //"lift" the differential from (r,q,s) -> (0,q-r+1,s+1) over
  // projection_left_img
  //(actually just a free presentation)
  MatrixQ lift = lift_from_free(sequence.get_prime(), diff_left,
                                projection_left_img, er_left_codomain);
  //std::cout << "lift:\n" << lift << "\n";

  AbelianGroup e2_0_q_s =
      sequence.get_e_2(TrigradedIndex(0, index_.q(), index_.s()));

  dim_t mon_rank = session_.get_monomial_rank(index_.p() - r_s);
  MatrixQ result_lift(mon_rank * e2_left_codomain.rank(), ker_right_domain.rank());

  MatrixQ inclusion_right_domain = sequence.get_inclusion(index_, r_);
  //std::cout << "inclusion_right_domain:\n" << inclusion_right_domain << "\n";
  AbelianGroup ker_left_domain =
      sequence.get_kernel(TrigradedIndex(r_s, index_.q(), index_.s()), r_);
  MatrixQ inclusion_left_domain =
      sequence.get_inclusion(TrigradedIndex(r_s, index_.q(), index_.s()), r_);
  //std::cout << "inclusion_left_domain:\n" << inclusion_left_domain << "\n";
  for (dim_t i = 0; i < mon_rank; i++) {
    // r_ is both the page number and the p of the transgression
    MatrixQ r_I =
        session_.get_r_operations(index_.p(), static_cast<deg_t>(r_), i);
    //std::cout << "r_I:\n" << r_I << "\n";
    // obtain the tensor product A\otimes r_I, where A is the group e2_0_q_s.
    MatrixQ r_I_q(r_I.height() * e2_0_q_s.rank(),
                  r_I.width() * e2_0_q_s.rank());
    for (dim_t diag = 0; diag < e2_0_q_s.rank(); diag++) {
      for (dim_t height = 0; height < r_I.height(); height++) {
        for (dim_t width = 0; width < r_I.width(); width++) {
          r_I_q(diag * r_I.height() + height, diag * r_I.width() + width) =
              r_I(height, width);
        }
      }
    }
    //std::cout << "r_I_q:\n" << r_I_q << "\n";
    // lift r_I_q * inclusion_right_domain against
    // inclusion_left_domain (okay because this is injective).
    MatrixQ r_I_ker =
        lift_from_free(sequence.get_prime(), r_I_q * inclusion_right_domain,
                       inclusion_left_domain, ker_left_domain);
    //std::cout << "r_I_ker:\n" << r_I_ker << "\n";
    MatrixQ lift_r_I = lift * r_I_ker;
    //std::cout << "lift_r_I:\n" << lift_r_I << "\n";
    for (dim_t h = 0; h < e2_left_codomain.rank(); h++) {
      for (dim_t w = 0; w < ker_right_domain.rank(); w++) {
        result_lift(h * mon_rank + i, w) = lift_r_I(h, w);
      }
    }
  }

  MatrixQ projection_right_img = sequence.get_projection(
      TrigradedIndex(index_.p() - r_s, index_.q() + r_s - 1, index_.s() + 1),
      r_);

  diff_candidate_ = projection_right_img * result_lift;

  // now also determine indeterminacy:
  MatrixQ id = IdentityMatrix<mpq_class>(projection_left_img.width());
  MatrixQList from_X;
  from_X.emplace_back(id);
  GroupWithMorphisms ker_proj_morphisms =
      compute_kernel(sequence.get_prime(), projection_left_img, e2_left_codomain,
                     er_left_codomain, MatrixQRefList(), ref(from_X));
  MatrixQ from_K = *ker_proj_morphisms.maps_from.begin();
  MatrixQ from_K_tensor(from_K.height() * mon_rank, from_K.width() * mon_rank);
  for (dim_t i = 0; i < from_K.height(); i++) {
    for (dim_t j = 0; j < from_K.width(); j++) {
      for (dim_t d = 0; d < mon_rank; d++) {
        from_K_tensor(i * mon_rank + d, j * mon_rank + d) = from_K(i, j);
      }
    }
  }

  AbelianGroup coker_right_img = sequence.get_cokernel(
      TrigradedIndex(index_.p() - r_s, index_.q() + r_s - 1, index_.s() + 1),
      r_);
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

bool DifferentialTask::usersolve() {
  SpectralSequence& sequence = session_.get_sequence();

  AbelianGroup coker = sequence.get_cokernel(target(index_,r_),r_);
  GroupWithMorphisms domain = sequence.get_e_ab(index_, r_, index_.q()+2);
  std::stringstream filename;
  filename << "diff_task_" << index_.p() << "_" << index_.q() << "_" << index_.s() << "_" << r_ << ".dat";
  std::stringstream text;
  text << "Change the Matrix below to the differential d_"
  << r_
  <<" from ";
  domain.group.print(text, sequence.get_prime());
  text << " at " <<index_
  <<" to ";
  coker.print(text, sequence.get_prime());
  text << " at "<<target(index_,r_) << "\n";

  session_.matrix_file_dialog(
          coker.rank(),
          domain.group.rank(),
          filename.str(),
          text.str());
  MatrixQ diff = session_.read_matrix_file(coker.rank(), domain.group.rank(),filename.str());
  MatrixQ diff_from_ker = diff * domain.maps_to[0];
  sequence.set_diff(index_, r_, diff_from_ker);
  return true;
}

void DifferentialTask::display_detail() {
  std::cout << "DifferentialTask for a d_"
  << r_
  << " from (p,q,s) = ("
  << index_.p()
  << ","
  << index_.q()
  << ","
  << index_.s()
  << ")\n";
}

void DifferentialTask::display_overview() {
  std::cout << "DifferentialTask for a d_"
  << r_
  << " from (p,q,s) = ("
  << index_.p()
  << ","
  << index_.q()
  << ","
  << index_.s()
  << ")\n";
}
