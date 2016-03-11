#include "session.h"


void Session::parse_ranks(std::string path, std::size_t max_deg){
	std::ifstream file;
	file.open(path);
	for(int i=1; 2*i<=max_deg; i++){
		mpz_class rank;
		if(!parse_mpz_class(file, rank)){
			throw std::logic_error("Session::parse_ranks: file syntax error in file " + path
					 + ". Is max_deg bigger than the amount of ranks provided?");
		}
		ranks_.emplace(2*i, rank.get_ui());
	}
}

void Session::parse_v_inclusions(std::string path, std::size_t max_deg){
	std::ifstream file;
	file.open(path);
	for(int i=1; 2*i<=max_deg; i++){
		MatrixQ matrix;
		if(!parse_matrix(file, matrix)){
			throw std::logic_error("Session::parse_v_inclusions: file syntax error in file " + path
					 + ". Is max_deg bigger than the amount of matrices provided?");
		}
		v_inclusions_.emplace(2*i, matrix);
	}
}

void Session::parse_r_operations(std::string path, std::size_t domain_deg){
	if(domain_deg % 2 != 0){
		throw std::logic_error("Session::parse_r_operations: called with odd degree.");
	}
	std::ifstream file;
	file.open(path);
	for(int i=1; 2*i<= domain_deg - 2; i++){
		std::map<std::size_t, std::size_t>::iterator
		     rank = ranks_.find(domain_deg - 2*i);
		if(rank == ranks_.end()){
			throw std::logic_error("Session::parse_r_operations: rank not set. Call parse_ranks before this.");
		}
		for(int j=0; j<rank->second; j++){
			MatrixQ matrix;
			if(!parse_matrix(file, matrix)){
						throw std::logic_error("Session::parse_roperations: file syntax error in file " + path);
			}
			v_inclusions_.emplace(std::make_tuple(domain_deg, 2*i, j), matrix);
		}
	}
}

void Session::init(std::string ranks_path, std::string v_inclusions_path,
		           std::string r_operations_path_prefix, std::size_t max_deg){
	parse_ranks(ranks_path, max_deg);
	parse_v_inclusions(v_inclusions_path, max_deg);
	for(int i=1; 2*i <= max_deg; i++){
		parse_r_operations(r_operations_path_prefix+std::to_string(2*i),2*i);
	}
}

void Session::step(){
	generate_group_tasks();
	autosolve_tasks();

	for(int r=2; r<= current_q_+1; r++){
		generate_differential_tasks(r); //have to do them one page after another here since they depend on each other.
		autosolve_tasks();
		user_solve_tasks();
	}

	generate_extension_tasks();
	autosolve_tasks();
	user_solve_tasks();

	sequence_.set_bounds(current_q_+1, 1,current_q_+1);
	current_q_++;

}

void Session::generate_group_tasks(){
	for(int p=2; p<=current_q_+2; p++){ //implicitly, p=1 is already 0! And for p=1, the needed (0,q_+1) isn't known yet anyways.
		GroupTask task(sequence_, p, current_q_+2-p);
		task_list_.emplace_back(task);
	}
}

void Session::generate_differential_tasks(std::size_t r){
	for(int p=2; p<=current_q_+2; p++){

		if(r <= p-1){
			std::pair<std::size_t, std::size_t> bounds = sequence_.get_bounds(current_q_+2-p);
			for(int s=bounds.first; s<=bounds.second; s++){
				DifferentialTask task(sequence_,TrigradedIndex(p, current_q_+2-p, s), r);
				task_list_.emplace_back(task);
			}
		}
		if(current_q_ + 3 - p - r >= 0 ){
			std::pair<std::size_t, std::size_t> bounds = sequence_.get_bounds(current_q_+3 -p - r);
			for(int s=bounds.first; s<=bounds.second; s++){
				DifferentialTask task(sequence_,TrigradedIndex(p+r, current_q_+3-p - r, s), r);
				task_list_.emplace_back(task);
			}
		}
	}
}

void Session::generate_extension_tasks(){
	for(int s = 1; s <= current_q_+1; s++){
		ExtensionTask task(sequence_, current_q_+1, s);
		task_list_.emplace_back(task);
	}
}

void Session::autosolve_tasks(){
	for(Task t : task_list_){
		if(t.autosolve()){
			//todo: remove from list. how do we do this??
		}
	}
}

void Session::user_solve_tasks(){
	//call shell dialog with task list, which should eventually return a number
	//then call solve on the corresponding task, again triggering a shell dialog
	//this should eventually return a bool. If it is true, the task has been solved and
	//is removed from the list.
}
