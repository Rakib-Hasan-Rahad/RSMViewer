#include "simulation.h"

using namespace std;

Simulation::Simulation(StructuralMotif query, Parameters input_scoring)	{
	
	motif_template = query;
	scoring_function = input_scoring;
	
	srand(time(NULL));
	
	//	initilize base-pair frequency
	iso_frequency = vector<int> (18);
	iso_orientation = vector<int> (18);
	iso_edge_1 = vector<int> (18);
	iso_edge_2 = vector<int> (18);
	
	iso_frequency[0] = 0;	iso_orientation[0] = 1;	iso_edge_1[0] = (int) ('W');	iso_edge_2[0] = (int) ('W');
	iso_frequency[1] = 433;	iso_orientation[1] = 0;	iso_edge_1[1] = (int) ('W');	iso_edge_2[1] = (int) ('W');
	iso_frequency[2] = 475;	iso_orientation[2] = 1;	iso_edge_1[2] = (int) ('W');	iso_edge_2[2] = (int) ('H');
	iso_frequency[3] = 501;	iso_orientation[3] = 0;	iso_edge_1[3] = (int) ('W');	iso_edge_2[3] = (int) ('H');
	iso_frequency[4] = 562;	iso_orientation[4] = 1;	iso_edge_1[4] = (int) ('W');	iso_edge_2[4] = (int) ('S');
	iso_frequency[5] = 584;	iso_orientation[5] = 0;	iso_edge_1[5] = (int) ('W');	iso_edge_2[5] = (int) ('S');
	iso_frequency[6] = 598;	iso_orientation[6] = 1;	iso_edge_1[6] = (int) ('H');	iso_edge_2[6] = (int) ('W');
	iso_frequency[7] = 609;	iso_orientation[7] = 0;	iso_edge_1[7] = (int) ('H');	iso_edge_2[7] = (int) ('W');
	iso_frequency[8] = 639;	iso_orientation[8] = 1;	iso_edge_1[8] = (int) ('H');	iso_edge_2[8] = (int) ('H');
	iso_frequency[9] = 640;	iso_orientation[9] = 0;	iso_edge_1[9] = (int) ('H');	iso_edge_2[9] = (int) ('H');
	iso_frequency[10] = 674;	iso_orientation[10] = 1;	iso_edge_1[10] = (int) ('H');	iso_edge_2[10] = (int) ('S');
	iso_frequency[11] = 675;	iso_orientation[11] = 0;	iso_edge_1[11] = (int) ('H');	iso_edge_2[11] = (int) ('S');
	iso_frequency[12] = 746;	iso_orientation[12] = 1;	iso_edge_1[12] = (int) ('S');	iso_edge_2[12] = (int) ('W');
	iso_frequency[13] = 770;	iso_orientation[13] = 0;	iso_edge_1[13] = (int) ('S');	iso_edge_2[13] = (int) ('W');
	iso_frequency[14] = 800;	iso_orientation[14] = 1;	iso_edge_1[14] = (int) ('S');	iso_edge_2[14] = (int) ('H');
	iso_frequency[15] = 834;	iso_orientation[15] = 0;	iso_edge_1[15] = (int) ('S');	iso_edge_2[15] = (int) ('H');
	iso_frequency[16] = 919;	iso_orientation[16] = 1;	iso_edge_1[16] = (int) ('S');	iso_edge_2[16] = (int) ('S');
	iso_frequency[17] = 935;	iso_orientation[17] = 0;	iso_edge_1[17] = (int) ('S');	iso_edge_2[17] = (int) ('S');
	
	//	initilize base-stacking frequency
	stk_frequency = vector<int> (4);
	stk_type = vector<int> (4);
	
	stk_frequency[0] = 0;	stk_type[0] = (int) ('u');
	stk_frequency[1] = 188;	stk_type[1] = (int) ('d');
	stk_frequency[2] = 277;	stk_type[2] = (int) ('i');
	stk_frequency[3] = 376;	stk_type[3] = (int) ('o');
	
	//	initilize nucleotide frequency
	nuc_frequency = vector<int> (4);
	nuc_letter = vector<int> (4);
	
	nuc_frequency[0] = 0;	nuc_letter[0] = (int) ('A');
	nuc_frequency[1] = 200;	nuc_letter[1] = (int) ('C');
	nuc_frequency[2] = 500;	nuc_letter[2] = (int) ('G');
	nuc_frequency[3] = 800;	nuc_letter[3] = (int) ('U');
	
}

Simulation::~Simulation(void)	{
	return;
}

int Simulation::random_select(vector<int> frequency)	{
	int num = rand() % 1000;
	int i;
	for(i = 0; i < (int) frequency.size() - 1; ++ i)	{
		if(num >= frequency[i] && num < frequency[i + 1])	{
			return i;
		}
	}
	return frequency.size() - 1;
}


void Simulation::modify_interaction(StructuralMotif &motif_copy)	{
	
	int i;
	for(i = 0; i < (int) motif_copy.record_interaction_.size(); ++ i)	{
		if(motif_copy.record_interaction_[i][2] > 0)	{
			//	this interaction is a base-pair
			int k = random_select(iso_frequency);
			motif_copy.record_interaction_[i][3] = iso_edge_1[k];
			motif_copy.record_interaction_[i][4] = iso_edge_2[k];
			motif_copy.record_interaction_[i][5] = iso_orientation[k];
		}	else if(motif_copy.record_interaction_[i][2] < 0)	{
			//	this interaction is a base-stacking
			int k = random_select(stk_frequency);
			motif_copy.record_interaction_[i][3] = stk_type[k];
		}
	}
	return;
}

void Simulation::modify_sequence(StructuralMotif& motif_copy)	{
	
	int i;
	for(i = 0; i < (int) motif_copy.sequence_.length(); ++ i)	{
		int k = random_select(nuc_frequency);
		motif_copy.sequence_[i] = (char) (nuc_letter[k]);
	}
	return;
}

void Simulation::modify_loop_size(StructuralMotif &motif_copy)	{
	//	create insertion
	float insertion_frequency = 0.2;
	int num_insertion = (int) (motif_copy.sequence_.length() * insertion_frequency);
	int i;
	int already_inserted = 0;
	while(already_inserted < num_insertion)	{
		
		int k = (int) (rand() % (motif_copy.sequence_.length() - 1)) + 1;
		int k_nuc = random_select(nuc_frequency);
		char c = (char) nuc_letter[k_nuc];
		motif_copy.sequence_ = motif_copy.sequence_.insert(k, 1, c);
		
		//	update index of record_interaction
		for(i = 0; i < (int) motif_copy.record_interaction_.size(); ++ i)	{
			if(motif_copy.record_interaction_[i][0] >= k)	{
				++ motif_copy.record_interaction_[i][0];
			}
			if(motif_copy.record_interaction_[i][1] >= k)	{
				++ motif_copy.record_interaction_[i][1];
			}
		}
		
		//	update index of break_points
		for(i = 0; i < (int) motif_copy.break_points_.size(); ++ i)	{
			if(motif_copy.break_points_[i] >= k)	{
				++ motif_copy.break_points_[i];
			}
		}

		//	update chain_stacking
		motif_copy.chain_stacking_.resize(motif_copy.chain_stacking_.size() + 1);
		for(i = motif_copy.chain_stacking_.size() - 2; i >= k; -- i)	{
			motif_copy.chain_stacking_[i + 1] = motif_copy.chain_stacking_[i];
		}
		int k_stk = random_select(stk_frequency);
		motif_copy.chain_stacking_[k] = motif_copy.hash_stacking(stk_type[k_stk]);
		
		++ already_inserted;
	}
	return;
}


void Simulation::reconsolidate_pair(StructuralMotif &motif_copy)	{
	int i;
	motif_copy.interaction_adjacency_ = vector< vector<int> > (motif_copy.sequence_.length());
	for(i = 0; i < (int) motif_copy.interaction_adjacency_.size(); ++ i)	{
		motif_copy.interaction_adjacency_[i] = vector<int> (motif_copy.sequence_.length(), 0);
	}
	
	for(i = 0; i < (int) motif_copy.record_interaction_.size(); ++ i)	{
		if(motif_copy.record_interaction_[i][2] > 0)	{
			//	if its a base pair
			motif_copy.record_interaction_[i][2] = motif_copy.hash_basepair(motif_copy.record_interaction_[i][5], 
				motif_copy.record_interaction_[i][3], motif_copy.record_interaction_[i][4],
				motif_copy.sequence_[motif_copy.record_interaction_[i][0]],
				motif_copy.sequence_[motif_copy.record_interaction_[i][1]]);
		}	else if(motif_copy.record_interaction_[i][2] < 0)	{
			//	if its a base-stacking
			motif_copy.record_interaction_[i][2] = motif_copy.hash_stacking(motif_copy.record_interaction_[i][3]);
		}
		//	update the interaction_adjacency array
		motif_copy.interaction_adjacency_[motif_copy.record_interaction_[i][0]][motif_copy.record_interaction_[i][1]]
			= motif_copy.record_interaction_[i][2];
		motif_copy.interaction_adjacency_[motif_copy.record_interaction_[i][1]][motif_copy.record_interaction_[i][0]]
			= motif_copy.record_interaction_[i][2];
	}
	return;
}

void Simulation::shuffle_motif(StructuralMotif &motif_copy)	{
	
	// shuffle the base pairs 
	modify_interaction(motif_copy);
	modify_sequence(motif_copy);
	modify_loop_size(motif_copy);
	reconsolidate_pair(motif_copy);
	return;
}

float _mean_vector(vector<int> numbers)	{
	
	float total = 0.0;
	int i;
	for(i = 0; i < (int) numbers.size(); ++ i)	{
		total += (float) numbers[i];
	}
	float mean = total / (float) numbers.size();
	return mean;
	
}

float _std_vector(vector<int> numbers, float mean)	{
	
	float total = 0.0;
	int i;
	for(i = 0; i < (int) numbers.size(); ++ i)	{
		total += ((float) numbers[i] - mean) * ((float) numbers[i] - mean);
	}
	float std = sqrt(total / (float) numbers.size());
	return std;
}

void Simulation::estimate_mean_std(void)	{
	int i;
	vector<int> simulated_score;
	double prev_mean, prev_std;
	mean = 0.0;
	std = 0.0;
	
	do	{
		prev_mean = mean;
		prev_std = std;
		int index = simulated_score.size();
		simulated_score.resize(simulated_score.size() + 100);
		for(i = 0; i < 100; ++ i)	{
			StructuralMotif random_motif = motif_template;
			shuffle_motif(random_motif);
			MotifGraphMatching simulate_alignment(motif_template, random_motif, scoring_function);
			simulate_alignment.align_simulation();
			simulated_score[index + i] = simulate_alignment.lower_bound;	
		}
		mean = _mean_vector(simulated_score);
		std = _std_vector(simulated_score, mean);
	}	while(abs(mean - prev_mean) > 1000 || abs(std - prev_std) > 1000);
	return;
}

double Simulation::chebyshev_inequality(double score)	{
	
	double estimate_pvalue = 1 / (((score - mean) / std) * ((score - mean) / std));
	if(score < 0 || estimate_pvalue > 1)	{
		estimate_pvalue = 1;
	}
	return estimate_pvalue;
	
}

