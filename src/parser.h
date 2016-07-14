#pragma once

#include <gmpxx.h>
#include <istream>
#include <string>
#include <list>

#include "matrix.h"
#include "types.h"
#include "abelian_group.h"
#include "p_local.h"

std::string read(std::istream &stream, std::streamsize count);

bool parse_mpz_class(std::istream& str, mpz_class& result);
bool parse_mpq_class(std::istream& str, mpq_class& result);
bool parse_matrix(std::istream& str, MatrixQ& result);
bool parse_matrix_size(dim_t height, dim_t width, std::istream& str, MatrixQ& result);
void eat_whitespace(std::istream& str);
bool accept_string(std::istream&input, std::string string);
bool parse_abelian_group(std::istream& str, AbelianGroup& result, mod_t p);
bool find_blank_line(std::istream& str);