#include "session.h"

Session::Session(std::string ranks_path, std::string v_inclusions_path,
                 std::string r_operations_path_prefix, dim_t max_deg)
{
  parse_ranks(ranks_path, max_deg);
  parse_v_inclusions(v_inclusions_path, max_deg);
  for (dim_t i = 2; 2 * i <= max_deg; i++) {
    parse_r_operations(r_operations_path_prefix + std::to_string(2 * i), 2 * i);
  }

  sequence_.set_bounds(0, 0, 0);
  sequence_.set_e2(TrigradedIndex(0, 0, 0), AbelianGroup(1, 0));
  current_q_ = 0;
}

SpectralSequence& Session::get_sequence()
{
  return sequence_;
}

dim_t Session::get_monomial_rank(deg_t p) const
{
  if (p < 0) throw std::logic_error("Session::get_monomial_rank: p<0");

  dim_t p_u = static_cast<dim_t>(p);

  if (p_u % 2 == 0) {
    if (2 * ranks_.size() <= p_u) {
      throw std::logic_error(
          "Session::get_monomial_rank: rank not set. Call parse_ranks before "
          "this.");
    }
    return ranks_[p_u / 2 - 1];
  } else {
    return 0;
  }
}

MatrixQ Session::get_v_inclusion(deg_t p) const
{
  if (p < 0) throw std::logic_error("Session::get_v_inclusion: p<0");

  dim_t p_u = static_cast<dim_t>(p);

  if (p_u % 2 == 0) {
    if (2 * v_inclusions_.size() <= p_u) {
      throw std::logic_error(
          "Session::get_v_inclusions: v_inclusions not set. Call "
          "parse_v_inclusions before "
          "this.");
    }
    return v_inclusions_[p_u / 2];
  } else {
    return MatrixQ(0, 0);
  }
}

void Session::parse_ranks(std::string path, dim_t max_deg)
{
  std::ifstream file;
  file.open(path);
  eat_whitespace(file);
  for (dim_t degree = 2; degree <= max_deg; degree += 2) {
    mpz_class rank;
    if (!parse_mpz_class(file, rank)) {
      throw std::logic_error(
          "Session::parse_ranks: file syntax error in file " + path +
          ". Is max_deg bigger than the amount of ranks provided?");
    }
    eat_whitespace(file);
    ranks_.emplace_back(rank.get_ui());
  }
}

void Session::parse_v_inclusions(std::string path, dim_t max_deg)
{
  std::ifstream file;
  file.open(path);
  eat_whitespace(file);
  for (dim_t degree = 2; degree <= max_deg; degree += 2) {
    MatrixQ matrix;
    if (!parse_matrix(file, matrix)) {
      throw std::logic_error(
          "Session::parse_v_inclusions: file syntax error in file " + path +
          ". Is max_deg bigger than the amount of matrices provided?");
    }
    v_inclusions_.emplace_back(matrix);
  }
}

void Session::parse_r_operations(std::string path, dim_t domain_deg)
{
  if (domain_deg % 2 != 0) {
    throw std::logic_error(
        "Session::parse_r_operations: called with odd degree.");
  }
  std::ifstream file;
  file.open(path);
  if (!file) {
    throw std::logic_error("Session:parse_r_operations: file " + path + " not "
        "found");
  }
  eat_whitespace(file);

  for (dim_t target_deg = 2; target_deg <= domain_deg - 2; target_deg += 2) {
    dim_t rank = get_monomial_rank(static_cast<deg_t>(domain_deg - target_deg));

    for (dim_t j = 0; j < rank; j++) {
      MatrixQ matrix;
      if (!parse_matrix(file, matrix)) {
        throw std::logic_error(
            "Session::parse_r_operations: file syntax error in file " + path);
      }
      r_operations_.emplace(std::make_tuple(domain_deg, target_deg, j), matrix);
    }
  }
}

void Session::step()
{
  generate_group_tasks();
  autosolve_tasks();

  sequence_.set_bounds(current_q_ + 1, 1, current_q_ + 1);

  // current_q_ is always non-negative
  for (dim_t r = 2; r <= static_cast<dim_t>(current_q_ + 1); r++) {
    generate_differential_tasks(r);  // have to do them one page after another
                                     // here since they depend on each other.
    autosolve_tasks();
    user_solve_tasks();
  }

  generate_extension_tasks();
  autosolve_tasks();
  user_solve_tasks();


  current_q_++;
}

void Session::generate_group_tasks()
{
  for (deg_t p = 2; p <= current_q_ + 2; p++) {
    // implicitly, p=1 is already 0! And for p=1, the needed (0,q_+1)
    // isn't known yet anyways.
    task_list_.emplace_back(new GroupTask(*this, p, current_q_ + 2 - p));
  }
}

void Session::generate_differential_tasks(dim_t r)
{
  for (deg_t p = 2; p <= current_q_ + 2; p++) {
    if (r <= static_cast<dim_t>(p - 1)) {
      std::pair<deg_t, deg_t> bounds = sequence_.get_bounds(current_q_ + 2 - p);
      for (deg_t s = bounds.first; s <= bounds.second; s++) {
        task_list_.emplace_back(new DifferentialTask(
            *this, TrigradedIndex(p, current_q_ + 2 - p, s), r));
      }
    }
    if (r <= static_cast<dim_t>(current_q_ + 3 - p)) {
      deg_t r_s = static_cast<deg_t>(r);
      std::pair<deg_t, deg_t> bounds =
          sequence_.get_bounds(current_q_ + 3 - p - r_s);
      for (deg_t s = bounds.first; s <= bounds.second; s++) {
        task_list_.emplace_back(new DifferentialTask(
            *this, TrigradedIndex(p + r_s, current_q_ + 3 - p - r_s, s), r));
      }
    }
  }
}

void Session::generate_extension_tasks()
{
  for (deg_t s = 1; s <= current_q_ + 1; s++) {
    task_list_.emplace_back(new ExtensionTask(*this, current_q_ + 1, s));
  }
}

void Session::autosolve_tasks()
{
  auto task_it = task_list_.begin();
  while (task_it != task_list_.end()) {
    if ((*task_it)->autosolve()) {
      task_it = task_list_.erase(task_it);
    } else {
      ++task_it;
    }
  }
}

void Session::user_solve_tasks()
{
  if (!task_list_.empty()) {
    throw std::logic_error("User-interaction needed.");
  }
  // call shell dialog with task list, which should eventually return a number
  // then call solve on the corresponding task, again triggering a shell dialog
  // this should eventually return a bool. If it is true, the task has been
  // solved and
  // is removed from the list.
}
