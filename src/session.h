#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include "parser.h"
#include "spectral_sequence.h"
#include "task.h"

class Session
{
 public:
  void init(std::string ranks_path, std::string v_inclusions_path,
            std::string r_operations_path_prefix, std::size_t max_deg);
  void step();

 private:
  SpectralSequence sequence_;
  // shell

  std::size_t current_q_;

  std::map<std::size_t, std::size_t> ranks_;
  std::map<std::tuple<std::size_t, std::size_t, std::size_t>, MatrixQ>
      r_operations_;
  // at (p, k, i), we find the i'th operation from deg p to deg k.
  std::map<std::size_t, MatrixQ> v_inclusions_;

  std::vector<Task> task_list_;
  void parse_ranks(std::string path, std::size_t max_deg);
  void parse_v_inclusions(std::string path, std::size_t max_deg);
  void parse_r_operations(std::string path, std::size_t domain_deg);

  void generate_group_tasks();
  void generate_differential_tasks(std::size_t r);
  void generate_extension_tasks();
  void autosolve_tasks();
  void user_solve_tasks();
};
