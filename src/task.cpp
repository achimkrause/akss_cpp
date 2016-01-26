#include "task.h"
#include "spectral_sequence.h"
#include "morphisms.h"

GroupTask::GroupTask(SpectralSequence& sequence, const TrigradedIndex index)
 : sequence_(sequence), index_(index)
{
}


bool GroupTask::solve(){
	AbelianGroup null_q_s = sequence_.get_e_2(TrigradedIndex(0,index_.q(),index_.s()));

	int mon_rank; //compute rank of Z[l_i] in degree p here!

	AbelianGroup result(null_q_s.free_rank()*mon_rank, null_q_s.tor_rank()*mon_rank);
	for(int i=0; i<null_q_s.tor_rank(); i++){
		for(int j=0; j<mon_rank; j++){
			result(i*mon_rank +j) = null_q_s(i);
		}
	}
	sequence_.set_e_2(index_,result);
	return true;
}


ExtensionTask::ExtensionTask(SpectralSequence& sequence, int q, int s)
  : sequence_(sequence),q_(q),s_(s) {
}

bool ExtensionTask::solve(){

	std::map<int, AbelianGroup> list_groups;
	for(int n=1; n<=q_;n++){
		//take the term at (n,q-n+1,s-1) taking into account differentials up to n-1 leaving,
		//and differentials up to q-n+2 entering.
		AbelianGroup eab = sequence_.get_e_ab(TrigradedIndex(n, q_-n+1, s_-1), n-1,q_-n+2);
		if(eab.rank()!=0) {
			list_groups.insert(std::pair<int,AbelianGroup>(n,eab));
		}
	}

	//now, for n=q_+1, s=1, we have to compute e_ab mod the v_n.
	if(s_==1){
	 AbelianGroup iterated_kernel = sequence_.get_kernel(TrigradedIndex(q_+1,0,0),q_+1);
     MatrixQ inclusion = sequence_.get_inclusion(TrigradedIndex(q_+1,0,0),q_+1);

     MatrixQ v_i_map; //get from nat's table

	 MatrixQ matrix = lift_from_free(sequence_.get_prime(),v_i_map, inclusion, iterated_kernel); //compute lift of v_i_map along inclusion.

	 GroupWithMorphisms coker = compute_cokernel(sequence_.get_prime(),matrix, iterated_kernel, MatrixQRefList(), MatrixQRefList());
	 list_groups.insert(std::pair<int,AbelianGroup>(q_+1,coker.group));
	}


	if(list_groups.size()==0){
		sequence_.set_e_2(TrigradedIndex(0,q_,s_),AbelianGroup(0,0));
		//set differentials to 0, or in sparse differential storage, set upper bounds of lists up.
	}
	else if(list_groups.size()==1){
		sequence_.set_e_2(TrigradedIndex(0,q_,s_),list_groups.begin()->second);
		int r=list_groups.begin()->first;
		//set the first r-1 differentials to 0, the r'th to the projection from
		// r-1'st kernel -> r-1'st kernel / image of all differentials
	}
	else {
		//shell: generate dialog for an extension problem with list_groups
		//get a group from the shell, and then ask for differential maps
		//if we encounter a problem somewhere, stop with false.
		//store progress in groups_ and differentials_
		//if we fail, there should be an option to revert back to stage r, or completely.
		//if it works, then write all of it into the spectral sequence.
	}
	return true;
}

DifferentialTask::DifferentialTask(SpectralSequence& sequence, TrigradedIndex index, int r)
 : sequence_(sequence), index_(index), r_(r){
}

bool DifferentialTask::solve(){
	//First, lift the map d_r: (r,q,s) -> (0,q-r+1,s+1) to a map lift between frees
	//(this just means lifting it against the projection E2 -> r-1'st cokernel), and then taking the corresponding matrix)
	//For each sequence I in the appropriate degree, compute the induced map on r-1'st kernels
	//(this is lifting against inclusion, we use nat's matrices for r_I), and then view
	//lift\circ r_I as a map from an integral lift of the r-1_st kernel at (p,q,s) to E2 at (0,..).
	//By filling in the entries of all these matrices correctly, we can build
	//a big matrix for the whole differential. This is a lift of the differential to E2
	//Now defining K as the kernel of the projection E2 -> r-1st cokernel (or rather, the free lift of E2),
	//we can compute the composition K\tensor Z{l_I} -> E2\tensor Z{l_I} -> r-1st cokernel.
	//this is our indeterminacy.
	//The composition of our big matrix with the projection onto the r-1st cokernel is our representative.


	AbelianGroup e2_left_img = sequence_.get_e_2(TrigradedIndex(0,index_.q()+r_-1,index_.s()+1));
	AbelianGroup er_left_img = sequence_.get_cokernel(TrigradedIndex(0,index_.q()+r_-1,index_.s()+1),r_);
	MatrixQ projection_left_img = sequence_.get_projection(TrigradedIndex(0,index_.q()+r_+1,index_.s()+1),r_);

	MatrixQ lift = lift_from_free(sequence_.get_prime(), /*placeholder for differential from (r,q,s) -> (0,q-r+1,s+1)*/,
			                        projection_left_img, er_left_img);
	  //"lift" the differential from (r,q,s) -> (0,q-r+1,s+1) over projection_left_img
      //(actually just a free presentation)

	AbelianGroup e2_0_q_s = sequence_.get_e_2(TrigradedIndex(0,index_.q(), index_.s()));
	AbelianGroup ker_right_domain = sequence_.get_kernel(index_, r_);
	int mon_rank; //compute rank of Z[l_i] in degree p-r here!

	MatrixQ result_lift(mon_rank*e2_left_img.rank(), ker_right_domain.rank());

	MatrixQ inclusion_right_domain = sequence_.get_inclusion(index_, r_);
	AbelianGroup ker_left_domain = sequence_.get_kernel(TrigradedIndex(r_,index_.q(), index_.s()),r_);
	MatrixQ inclusion_left_domain = sequence_.get_inclusion(TrigradedIndex(r_,index_.q(), index_.s()),r_);

	for(int i=0; i<mon_rank; i++){
		MatrixQ r_I; //obtain the i_th operation from degree p to degree r here from nat's tables.
		MatrixQ r_I_q(r_I.height()*e2_0_q_s.rank(),r_I.width()*e2_0_q_s.rank()); //obtain the tensor product A\otimes r_I, where A is the group e2_0_q_s.
		for(int diag=0; diag<e2_0_q_s.rank(); diag++){
			for(int height=0; height<r_I.height();height++){
				for(int width=0; width<r_I.width();width++){
					r_I_q(diag*r_I.height()+height, diag*r_I.width()+width)= r_I(height,width);
				}
			}
		}


		//lift r_I_q * inclusion_right_domain against
		//inclusion_left_domain (okay because this is injective).
		MatrixQ r_I_ker= lift_from_free(sequence_.get_prime(),
				 r_I_q*inclusion_right_domain, inclusion_left_domain, ker_left_domain);

		MatrixQ lift_r_I = lift * r_I_ker;
		for(int h; h<e2_left_img.rank(); h++){
			for(int w; w<ker_right_domain.rank(); w++){
				result_lift(h*mon_rank + i, w)=r_I_ker(h,w);
			}
		}
	}

	MatrixQ projection_right_img = sequence_.get_projection(TrigradedIndex(index_.p()-r_,index_.q()-r_+1,index_.s()+1),r_);

	MatrixQ differential = projection_right_img * result_lift;


	//now also determine indeterminacy:
	MatrixQ id = IdentityMatrix<mpq_class>(projection_right_img.width());
	MatrixQList from_X;
	from_X.emplace_back(id);
	GroupWithMorphisms ker_proj_morphisms = compute_kernel(sequence_.get_prime(), projection_left_img,e2_left_img, er_left_img,MatrixQRefList(), ref(from_X));
	MatrixQ from_K = ker_proj_morphisms.maps_from.begin();
	MatrixQ from_K_tensor(from_K.height * mon_rank, from_K.width*mon_rank);
	for(int i=0; i<from_K.height; i++){
		for(int j=0; j<from_K.width; j++){
			for(int d=0; d<mon_rank; d++){
				from_K_tensor(i*mon_rank + d, j*mon_rank + d) = from_K(i,j);
			}
		}
	}

	AbelianGroup coker_right_img = sequence_.get_cokernel(TrigradedIndex(index_.p()-r_, index_.q()+r_-1, index_.s()+1),r_);
	MatrixQ indet_map = projection_right_img*from_K_tensor;
	GroupWithMorphisms indeterminacy = compute_image(sequence_.get_prime(),indet_map,AbelianGroup(indet_map.width(),0),coker_right_img); //compute image of indet_map.
	if(indeterminacy.group.rank()==0){
		//enter differential as the corresponding differential
	}
	else {
		//indeterminacy.maps_from.begin() columns generate indeterminacy.
		//prompt the user with both the indeterminacy  and the candidate, ask for a matrix from shell.
	}
	return true;
}


