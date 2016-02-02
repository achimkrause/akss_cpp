#pragma once

#include <map>
#include "spectral_sequence.h"
#include "morphisms.h"


class Task {
public:
	bool autosolve();

private:

	//todo: find reasonable parent class for all Tasks.
};


class GroupTask
{
public:
	bool autosolve();
	GroupTask(SpectralSequence& sequence, const std::size_t p, const std::size_t q);
	//computes E^2_{p,q,s} from E^2_{0,q,s} by tensoring with degree p in Z[l_i]
private:
	SpectralSequence& sequence_;
	//reference to shell

	std::size_t p_;
	std::size_t q_;
};

class DifferentialTask
{
public:
	bool autosolve();
	 //tries to compute the differential from (p,q,s) to (p-r,q+r-1,s+1) using operations and
	 //knowledge about the differential from (r,q,s) to (0,q+r-1,s+1).
	 //If that doesn't work, it displays the partial info obtained,
	 //requests a matrix from the shell.
    DifferentialTask(SpectralSequence& sequence, TrigradedIndex index, int r);
private:
	SpectralSequence& sequence_;
	//reference to shell

	TrigradedIndex index_;
	std::size_t r_;
	GroupWithMorphisms indeterminacy_;
	MatrixQ diff_candidate_;
};

class ExtensionTask
{
public:
	bool autosolve();
	 //gathers all the groups in degrees (n, q-n+1, s-1). If at most one of them is nonzero,
	 //fill in that group at (0,q,s) and take the differential to be the identity.
	 //otherwise,
	 //requests a group from the shell, to fill into spot (0,q,s).
	 //then requests differentials as necessary.
	 //the progress in doing that is saved in groups and differentials:
	 //    in groups, we store the subsequent quotients of the group at (0,q,s) as we enter
	 //                                                                 differentials
	 //    in differentials, we store the maps themselves.
	 //Progress could be stored and revisited later, but may only be entered into
	 //the spectral sequence once it is completed, because it is only valid if it gives the
	 //correct E^\infty page entries.

	ExtensionTask(SpectralSequence& sequence, int q, int s)
private:
	SpectralSequence& sequence_;
    //reference to shell

	std::size_t q_;
    std::size_t s_;

    GroupSequence groups_;
    MatrixQList differentials_;

	std::map<int, AbelianGroup> list_groups_; //the groups remaining after present differentials
	                                         //of the form iterated kernel / subgrp
	std::map<int, MatrixQ> list_maps_;        //the maps iterated kernel -> iterated kernel/subgrp.


};
