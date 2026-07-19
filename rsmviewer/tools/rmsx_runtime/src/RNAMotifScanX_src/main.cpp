#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <cstdlib>
#include <unistd.h>
#include "./inc/simulation.h"

using namespace std;


alignment_parameters set_default;

void print_usage(void)	{
	cout << "Usage: RNAMotifAlignX -i motif_1 -j motif_2 [-s -t -n -a -b -c -d -r -u]" << endl;
	cout << endl;
	cout << "Required parameters:" << endl << endl;
	cout << "	motif_1 (-i):" << "	" << "The first motif to be aligned" << endl;
	cout << "	motif_2 (-j):" << "	" << "The second motif to be aligned" << endl;
	cout << endl;
	cout << "Optional parameters:" << endl << endl;
	cout << "	iso_matrix (-s):" << "	" << "The substitution scoring function of basepairs based on isostericity" << endl;
	cout << "	stacking_matrix (-t):" << "	" << "The substitution scoring function of base stackings" << endl;
	cout << "	nucleotide_matrix (-n):" << "	" << "The substitution scoring function of nucleotides" << endl;
	cout << "	weight_isosteric (-a):" << "	" << "The weight for matching basepairs with isostericity (default 2.0)" << endl;
	cout << "	weight_nonadjacent (-b):" << "	" << "The weight for matching non-adjacent base stackings (default 1.0)" << endl;
	cout << "	weight_adjacent (-c):" << "	" << "The weight for matching adjacent base stackings (default 0.1)" << endl;
	cout << "	weight_sequence (-d):" << "	" << "The weight for matching nucleotides (default 0.1)" << endl;
	cout << "	use heuristic (-r, no parameter):" << "	" << "Use heuristic method CLCL to align motif instances (default NO)" << endl;
	cout << "	specify a p-value cutoff (-u):" << "	" << "RNAMotifAlignX will only output alignment if it passes the cutoff" << endl;
	cout << "	help info (-h):" << "	" << "Print this information" << endl;
	return;
}


void initpar(void)	{

	set_default.isosteric_file_name = "./Isosteric_scoring_function.txt";
	set_default.stacking_file_name = "./Stacking_scoring_function.txt";
	set_default.nucleotide_file_name = "./Nucleotide_scoring_function.txt";

	set_default.gap_opening = -18.0;
	set_default.gap_extension = -16.0;
	set_default.triple_interaction_bonus = 3.0;
	set_default.hbond_match_base = 1.0;
	set_default.hbond_match_bonus = 2.0;
	set_default.weight_isosteric = 2.0;
	set_default.weight_nonadjacent_stacking = 1.0; 
	set_default.weight_adjacent_stacking = 0.2;
	set_default.weight_sequence = 0.1;
	
	set_default.pair_deletion = -12;
	set_default.stacking_deletion = -2;
	
	set_default.asym_nuc = -0.2; 
	set_default.asym_loop = -2; 
	set_default.cons_stacking = 0.5; 
	set_default.stack_to_bulge = -1.5;
	set_default.stack_to_internal_asym = -1;
	set_default.stack_to_internal_sym = -0.5;	
}

string motif_1, motif_2;
bool heuristic;
float pvalue_cutoff = 1.0;
bool use_cutoff = false;
bool search_mode = false;

void initenv(int argc, char **argv)	{
	//initialize environment variable
	int opt = 0;
	while((opt=getopt(argc, argv, "i:j:s:t:n:a:b:c:d:u:hre")) != -1) {
		switch(opt) {
			case 'i': {motif_1 = string(optarg); break;}
			case 'j': {motif_2 = string(optarg); break;}
			case 's': {set_default.isosteric_file_name = string(optarg); break;}
			case 't': {set_default.stacking_file_name = string(optarg); break;}
			case 'n': {set_default.nucleotide_file_name = string(optarg); break;}
			case 'a': {set_default.weight_isosteric = atof(optarg); break;}
			case 'b': {set_default.weight_nonadjacent_stacking = atof(optarg); break;}
			case 'c': {set_default.weight_adjacent_stacking = atof(optarg); break;}
			case 'd': {set_default.weight_sequence = atof(optarg); break;}
			case 'u': {pvalue_cutoff = atof(optarg); use_cutoff = true; break;}
			case 'r': {heuristic = true; break;}
			case 'e': {search_mode = true; break;}
			case 'h': {print_usage();	exit(0);}
			default: {print_usage();	exit(0);}
		}
	}
	if(motif_1.empty() || motif_2.empty()){
		print_usage();
		exit(0);
	}
	return;
}



int main(int argc, char **argv)	{
	
	//	initialize environment variables	
	initpar();
	initenv(argc, argv);
	
	// TODO: initialize RNA structural motifs with a graph, which can be generated automatically with another class
	structural_motif aln_motif_1(motif_1);
	structural_motif aln_motif_2(motif_2);
	
	//	construct parameters object
	parameters scoring_function(set_default);
	
	
	simulation shuffle(aln_motif_1, scoring_function);
	shuffle.estimate_mean_std();
	
	int i;
	
//	for(i = 0; i < aln_motif_1.record_interaction.size(); ++ i)	{
//		cout << aln_motif_1.record_interaction[i][0] << "	" << aln_motif_1.record_interaction[i][1] << "	" << aln_motif_1.record_interaction[i][2] << "	" << aln_motif_1.record_interaction[i][3] << "	" << aln_motif_1.record_interaction[i][4] << "	" << aln_motif_1.record_interaction[i][5] << endl;
//	}
	if(aln_motif_1.record_interaction.size() == 0 || aln_motif_2.record_interaction.size() == 0)	{
		if(!search_mode)	{
			cout << "No structure contained in one of the RNA structure motif instances that are being aligned." << endl;
			cout << "Cannot generate alignment. Execution terminated." << endl;
		}
		return 0;
	}

	if(use_cutoff)	{
		motif_graph_matching align_motifs(aln_motif_1, aln_motif_2, scoring_function, shuffle.mean, shuffle.std, pvalue_cutoff);
		if(heuristic)	{
			align_motifs.align_heuristic();
		}	else	{
			align_motifs.align();
		}
		
		if(search_mode)	{
			align_motifs.print_alignment_search();
		}	else	{
			align_motifs.print_alignment(align_motifs.maximal_clique);
		}
		
	}	else	{
		motif_graph_matching align_motifs(aln_motif_1, aln_motif_2, scoring_function, shuffle.mean, shuffle.std);
		if(heuristic)	{
			align_motifs.align_heuristic();
		}	else	{
			align_motifs.align();
		}
		
		if(search_mode)	{
			align_motifs.print_alignment_search();
		}	else	{
			align_motifs.print_alignment(align_motifs.maximal_clique);
		}
	}
//	for(i = 0; i < aln_motif_1.record_interaction.size(); ++ i)	{
//		cout << aln_motif_1.record_interaction[i][0] << "	" << aln_motif_1.record_interaction[i][1] << "	" << aln_motif_1.record_interaction[i][2] << "	" << aln_motif_1.record_interaction[i][3] << "	" << aln_motif_1.record_interaction[i][4] << "	" << aln_motif_1.record_interaction[i][5] << endl;
//	}
	
	return 0;
}
