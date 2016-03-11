#include "parser.h"

std::string read(std::istream& stream, std::size_t count)
{
    std::string result(count, ' ');
    stream.read(&result[0], count);

    return result;
}

bool parse_mpz_class(std::istream& str, mpz_class& result) {
	std::istream::streampos pos = str.tellg();

	int n=0;
	if(str.peek()=='-'){
		str.ignore();
		n++;
	}
	if(48<=str.peek() && str.peek()<=57){
		str.ignore();
		n++;
	}
	else {
		str.seekg(pos);
		return false;
	}
	while(48<=str.peek() && str.peek()<=57){
		str.ignore();
		n++;
	}

	str.seekg(pos);
	std::string s = read(str,n);
	result= mpz_class(s);
	result.set_str(s,10);
	return true;
}

bool parse_mpq_class(std::istream& str, mpq_class& result) {
	std::istream::streampos pos = str.tellg();
	mpz_class num;
	if(!parse_mpz_class(str,num)){
		return false;
	}
	if(str.peek()=='/'){
		str.ignore();
	}
	else{
		str.seekg(pos);
		return false;
	}
	mpz_class denom;
	if(!parse_mpz_class(str,denom)){
		str.seekg(pos);
		return false;
	}
	result = mpq_class(num,denom);
	return true;
}

void eat_whitespace(std::istream& str){
	while(str.peek() == 9 || str.peek() == 10 || str.peek() == 13 || str.peek() == 32){
		str.ignore();
	}
}

bool parse_matrix(std::istream& str, MatrixQ& result)
{
	std::istream::streampos pos = str.tellg();
	mpz_class height;
    mpz_class width;
	if(!parse_mpz_class(str, height)){
		str.seekg(pos);
		return false;
	}
	eat_whitespace(str);

	if(!parse_mpz_class(str, width)){
		str.seekg(pos);
		return false;
	}
	eat_whitespace(str);

	if(height < 0 || width < 0) {
		throw std::logic_error("parse_matrix: negative dimensions provided.");
	}
	MatrixQ result_tmp(height.get_ui(),width.get_ui());

	for(int i=0; i<height; i++){
		for(int j=0; j<height; j++){
			if(!parse_mpq_class(str, result_tmp(i,j))){
				str.seekg(pos);
				return false;
			}
			eat_whitespace(str);
		}
	}

	result = result_tmp;
	return true;
}
