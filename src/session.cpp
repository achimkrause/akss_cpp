#include "session.h"
#include <sstream>

Session::Session(mod_t prime, std::string ranks_path, std::string v_inclusions_path,
                 std::string r_operations_path_prefix, dim_t max_deg)
  : sequence_(prime)
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
    return v_inclusions_[p_u / 2 -1];
  } else {
    return MatrixQ(0, 0);
  }
}

MatrixQ Session::get_r_operations(deg_t source, deg_t target, dim_t index) const
{
  auto r_operations_it = r_operations_.find(std::make_tuple(source, target, index));
  if(r_operations_it == r_operations_.end()){
    std::stringstream str;
    str << "Session::get_r_operations: Operation with index "
              << index
              << " from "
              << source
              << " to "
              << target
              << " not set.";
    throw std::logic_error(str.str());
  }
  return r_operations_it->second;
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
  task_list_.emplace_back(new GroupTask(*this, 1, current_q_));
  task_list_.emplace_back(new GroupTask(*this, 2, current_q_));
  if(current_q_ > 0) {
    task_list_.emplace_back(new GroupTask(*this, 3, current_q_-1));
  }
  for (deg_t p = 4; p <= current_q_ + 3; p++) {
    task_list_.emplace_back(new GroupTask(*this, p, current_q_ + 3 - p));
  }
}

void Session::generate_differential_tasks(dim_t r)
{
  generate_differential_tasks_pq_deg(r+1, current_q_+1-r, r);


  for (deg_t p = 2+r; p<= current_q_ + 3; p++) {
    generate_differential_tasks_pq_deg(p, current_q_ + 3 - p, r);
  }
}

//generates differential tasks from (p,q,s) for all s within bounds of source OR target.
void Session::generate_differential_tasks_pq_deg(dim_t p, dim_t q, dim_t r){

  std::pair<deg_t, deg_t> bounds_source = sequence_.get_bounds(q);
  std::pair<deg_t, deg_t> bounds_target = sequence_.get_bounds(q+r-1);
  deg_t min = bounds_source.first;
  if(bounds_target.first-1 < bounds_source.first){
    min = bounds_target.first-1;
  }
  deg_t max = bounds_source.second;
  if(bounds_target.second-1 > bounds_source.second){
    max = bounds_target.second-1;
  }

  for(deg_t s = min; s<=max; s++){
    task_list_.emplace_back(new DifferentialTask(
            *this, TrigradedIndex(p,q,s), r));
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
      task_it++;
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

void Session::display_tasks_overview() {
  std::size_t i=0;
  for(auto task_it = task_list_.begin(); task_it != task_list_.end(); task_it++){
    std::cout << i << ": ";
    (*task_it)->display_overview();
  }
}

void Session::display_task_detail(std::size_t i) {
  auto task_it = task_list_.begin();
  for(std::size_t j=0; j< i; j++){
    task_it++;
    if(task_it == task_list_.end()){
      std::cout << "Index too high.\n";
      return;
    }
  }
  (*task_it)->display_detail();
}

void Session::display_command_overview() {
  std::cout << "details i: displays details for task i.\n";
  std::cout << "e r p q s: displays E^r page entry in tridegree (p,q,s).\n";
  std::cout << "e a b p q s: displays kernel of differentials up to d_{a-1} mod image of differentials up to d_{b-1} "
            << "at (p,q,s). e r r p q s is equivalent to e r p q s.\n";
  std::cout << "d r p q s: displays differential d_r leaving degree (p,q,s).\n";
  std::cout << "solve i: opens an editor with a template for the user input for task i.\n";
}

void Session::interact(){
  std::string str;
  std::getline(std::cin,str);
  std::stringstream input(str);
  if(accept_string(input, "details ")){
    mpz_class i;
    eat_whitespace(input);
    if(parse_mpz_class(input, i)) {
      eat_whitespace(input);
      if(input.eof()){
        display_task_detail(i.get_ui());
        return;
      }
    }
  }
  else if(accept_string(input, "e ")){
    std::size_t numbers[5];
    std::size_t numbers_parsed=0;
    for(;numbers_parsed<5; numbers_parsed++){
      mpz_class n;
      eat_whitespace(input);
      if(!parse_mpz_class(input,n)){
        break;
      }
      numbers[numbers_parsed]=n.get_ui();
    }
    eat_whitespace(input);
    if(numbers_parsed==5 && input.eof()){
      display_eab(numbers[0], numbers[1], numbers[2], numbers[3], numbers[4]);
      return;
    }
    if(numbers_parsed==4 && input.eof()){
      display_e(numbers[0], numbers[1], numbers[2], numbers[3]);
      return;
    }
  }
  else if(accept_string(input, "d ")){
    std::size_t numbers[4];
    std::size_t numbers_parsed=0;
    for(;numbers_parsed<4; numbers_parsed++){
      mpz_class n;
      eat_whitespace(input);
      if(!parse_mpz_class(input,n)){
        break;
      }
      numbers[numbers_parsed]=n.get_ui();
    }
    eat_whitespace(input);
    if(numbers_parsed==4 && input.eof()){
      display_differential(numbers[0],numbers[1],numbers[2],numbers[3]);
      return;
    }
  }
  else if(accept_string(input, "solve ")){
     //placeholder
  }
  std::cout << "Invalid syntax.\n";
  return;
}

void Session::display_eab(std::size_t a, std::size_t b, std::size_t p, std::size_t q, std::size_t s) {
  AbelianGroup eab = sequence_.get_e_ab(TrigradedIndex(p,q,s),a,b).group;
  //std::cout << eab << "\n";
}

void Session::display_e(std::size_t r, std::size_t p, std::size_t q, std::size_t s) {
  AbelianGroup e = sequence_.get_e_ab(TrigradedIndex(p,q,s),r,r).group;
  //std::cout << e << "\n";
}

void Session::display_differential(std::size_t r, std::size_t p, std::size_t q, std::size_t s) {
  MatrixQ d = sequence_.get_diff_from(TrigradedIndex(p,q,s),r);
  std::cout << d;
}