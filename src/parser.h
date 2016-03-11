#pragma once

#include <gmpxx.h>
#include <istream>
#include <string>

#include "matrix.h"

std::string read(std::istream &stream, std::size_t count);

bool parse_mpz_class(std::istream& str, mpz_class& result);
bool parse_mpq_class(std::istream& str, mpq_class& result);
bool parse_matrix(std::istream& str, MatrixQ& result);
void eat_whitespace(std::istream& str);
