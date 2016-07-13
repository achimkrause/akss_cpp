#pragma once

#include <iostream>
#include <fstream>
#include <list>
#include <map>
#include <memory>
#include <string>

#include "parser.h"
#include "spectral_sequence.h"
#include "task.h"
#include "types.h"

class Task;

class Session
{
 public:
  Session(mod_t prime, std::string ranks_path, std::string v_inclusions_path,
          std::string r_operations_path_prefix, dim_t max_deg);
  void step();

  SpectralSequence& get_sequence();
  dim_t get_monomial_rank(deg_t p) const;
  MatrixQ get_v_inclusion(deg_t p) const;
  MatrixQ get_r_operations(deg_t source, deg_t target, dim_t index) const;

 private:
  void parse_ranks(std::string path, dim_t max_deg);
  void parse_v_inclusions(std::string path, dim_t max_deg);
  void parse_r_operations(std::string path, dim_t domain_deg);

  void generate_group_tasks();
  void generate_differential_tasks(dim_t r);
  void generate_differential_tasks_pq_deg(dim_t p, dim_t q, dim_t r);
  void generate_extension_tasks();
  void autosolve_tasks();
  void user_solve_tasks();

  //IO Stuff
  void display_tasks_overview();
  void display_task_detail(std::size_t i);
  void display_command_overview();
  void display_eab(std::size_t a, std::size_t b, std::size_t p, std::size_t q, std::size_t s);
  void display_e(std::size_t r, std::size_t p, std::size_t q, std::size_t s);
  void display_differential(std::size_t r, std::size_t p, std::size_t q, std::size_t s);
  void interact();
  void display_anss_e2();

  // shell
  SpectralSequence sequence_;

  deg_t current_q_;

  std::vector<dim_t> ranks_;
  // at (p, k, i), we find the i'th operation from deg p to deg k.

  std::map<std::tuple<deg_t, deg_t, dim_t>, MatrixQ> r_operations_;// <domain, codomain, number>
  std::vector<MatrixQ> v_inclusions_;

  std::list<std::unique_ptr<Task>> task_list_;
};
