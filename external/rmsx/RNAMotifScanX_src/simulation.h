#include "structural_motif.h"
#include "motif_graph_matching.h"

#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <time.h>

#ifndef _SIMULATION_H_
#define _SIMULATION_H_

class Simulation	{

 public:
	
	StructuralMotif motif_template;
	Parameters scoring_function;
	
	double mean, std;
	
	//	the default frequency for non-canonial base pairs, stacking interactions, and nucleotides	
	//	assume the frequencys are multiplied by 1000 to make them integer
	std::vector<int> iso_frequency;
	std::vector<int> iso_orientation, iso_edge_1, iso_edge_2;
	
	std::vector<int> stk_frequency;
	std::vector<int> stk_type;
	
	std::vector<int> nuc_frequency;
	std::vector<int> nuc_letter;
	
	
	Simulation(StructuralMotif query, Parameters input_scoring);
	~Simulation(void);
	int random_select(std::vector<int> frequency);
	//	randomly substitute the base interaction of a motif
	void modify_interaction(StructuralMotif &motif_copy);
	//	randomly generate the sequence of a motif
	void modify_sequence(StructuralMotif &motif_copy);
	//	create random insertions/deletions of a motif
	void modify_loop_size(StructuralMotif &motif_copy);
	//	remap the base-interaction after modifying base interactions and sequences
	void reconsolidate_pair(StructuralMotif &motif_copy);
	//	a wrapper function to shuffle a motif
	void shuffle_motif(StructuralMotif &motif_copy);
	//	estimate the parameter (mean and std) by simulation
	void estimate_mean_std(void);
	double chebyshev_inequality(double score);
};

#endif
