#include "parser.h"

std::string read(std::istream& stream, std::streamsize count)
{
  if (count <= 0)
    return "";

  std::string result(static_cast<std::size_t>(count), ' ');
  stream.read(&result[0], count);

  return result;
}

bool parse_mpz_class(std::istream& str, mpz_class& result)
{
  const std::istream::streampos pos = str.tellg();

  std::streamsize n = 0;
  if (str.peek() == '-') {
    str.ignore();
    n++;
  }
  if (48 <= str.peek() && str.peek() <= 57) {
    str.ignore();
    n++;
  } else {
    str.seekg(pos);
    return false;
  }
  while (48 <= str.peek() && str.peek() <= 57) {
    str.ignore();
    n++;
  }

  str.seekg(pos);
  std::string s = read(str, n);
  result = mpz_class(s);
  result.set_str(s, 10);
  return true;
}

bool parse_mpq_class(std::istream& str, mpq_class& result)
{
  const std::istream::streampos pos = str.tellg();
  mpz_class num;
  if (!parse_mpz_class(str, num)) {
    return false;
  }
  if (str.peek() != '/') {
    result = mpq_class(num, 1);
    return true;
  }

  str.ignore();
  mpz_class denom;
  if (!parse_mpz_class(str, denom)) {
    str.seekg(pos);
    return false;
  }
  result = mpq_class(num, denom);
  result.canonicalize();
  return true;
}

void eat_whitespace(std::istream& str)
{
  while (str.peek() == 9 || str.peek() == 10 || str.peek() == 13 ||
         str.peek() == 32) {
    str.ignore();
  }
}

bool parse_matrix(std::istream& input, MatrixQ& result)
{
  const std::istream::streampos pos = input.tellg();
  mpz_class height;
  mpz_class width;
  if (!parse_mpz_class(input, height)) {
    input.seekg(pos);
    return false;
  }
  eat_whitespace(input);

  if (!parse_mpz_class(input, width)) {
    input.seekg(pos);
    return false;
  }
  eat_whitespace(input);

  if (height < 0 || width < 0) {
    throw std::logic_error("parse_matrix: negative dimensions provided.");
  }

  MatrixQ result_tmp(height.get_ui(), width.get_ui());
  mpq_class tmp;
  for (dim_t i = 0; i < height; i++) {
    for (dim_t j = 0; j < width; j++) {
      if (!parse_mpq_class(input, tmp)) {
        input.seekg(pos);
        return false;
      }
      result_tmp(i,j) = tmp;
      eat_whitespace(input);
    }
  }
  result = result_tmp;
  return true;
}

bool parse_matrix_size(dim_t height, dim_t width, std::istream& input, MatrixQ& result){
  MatrixQ result_tmp(height, width);

  const std::istream::streampos pos = input.tellg();

  mpq_class tmp;
  for (dim_t i = 0; i < height; i++) {
    for (dim_t j = 0; j < width; j++) {
      if (!parse_mpq_class(input, tmp)) {
        input.seekg(pos);
        return false;
      }
      result_tmp(i,j) = tmp;
      eat_whitespace(input);
    }
  }
  result = result_tmp;
  return true;
}

bool accept_string(std::istream& input, std::string string){
  const std::istream::streampos pos = input.tellg();

  for(std::size_t i=0; i<string.length(); i++){
    if(input.peek() == string[i]){
      input.ignore();
    }
    else {
      input.seekg(pos);
      return false;
    }
  }
  return true;
}

bool parse_abelian_group(std::istream& input, AbelianGroup& grp, mod_t p){
  std::size_t free_rank = 0;
  std::size_t tor_rank = 0;
  std::list<std::size_t> orders;

  const std::istream::streampos pos = input.tellg();
  while(accept_string(input, "Z")){
    eat_whitespace(input);
    if(accept_string(input,"/")){
      eat_whitespace(input);
      mpz_class order;
      if(!parse_mpz_class(input, order)){
        input.seekg(pos);
        return false;
      }
      eat_whitespace(input);
      orders.emplace_back( static_cast<std::size_t>(p_val_z(p, order)));
      tor_rank++;
    }
    else if(tor_rank > 0){ //expect no Z after torsion.
      input.seekg(pos);
      return false;
    }
    else {
      free_rank++;
    }
    if(!accept_string(input, "+")){
      break;
    }
    eat_whitespace(input);
  }
  if(tor_rank ==0 && free_rank == 0){
    if(!accept_string(input, "0")){
      input.seekg(pos);
      return false;
    }
    else{
      grp = AbelianGroup(0,0);
      return true;
    }
  }
  grp = AbelianGroup(free_rank,tor_rank);
  std::size_t i=0;
  for(auto orders_it = orders.begin(); orders_it != orders.end(); orders_it++){
    grp(i) = (*orders_it);
    i++;
  }
  return true;
}

bool find_blank_line(std::istream& input){
  const std::istream::streampos pos = input.tellg();
  bool one_before= false;
  while(input.peek()!= -1){
    if(input.peek() == '\n'){
      if(one_before){
        input.ignore();
        return true;
      }
      else{
        one_before = true;
        input.ignore();
      }
    }
    else {
      one_before = false;
      input.ignore();
    }
  }
  input.seekg(pos);
  return false;
}