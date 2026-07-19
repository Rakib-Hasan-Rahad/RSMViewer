#include "motif_graph_matching.h"

using namespace std;

template<typename T> 
void MotifGraphMatching::ReduceSizeUnsorted(std::vector<T>& array_to_reduce, const int& index_to_delete)  {
  // removes the 'index_to_delete'th element from the array, and reduce the array size by 1
  // note that the order in the original array is no longer preserved
  if((int) array_to_reduce.size() <= index_to_delete)  {
    cerr << "MotifGraphMatching::ReduceSizeUnsorted:index_to_delete out of bound" << array_to_reduce.size() << "  " << index_to_delete << endl;
    return;
  }
  array_to_reduce[index_to_delete] = array_to_reduce.back();
  array_to_reduce.pop_back();
  return;
}

MotifGraphMatching::MotifGraphMatching(
    const StructuralMotif& M1, 
    const StructuralMotif& M2, 
    const Parameters& input_scoring_function
)  {  
  consider_pvalue = false;
  consider_cutoff = false;
  unsigned int i, j;
  lower_bound = max_score = -MAX;
  scoring_function = input_scoring_function;
  M1_seq = M1.sequence_;
  M2_seq = M2.sequence_;
  M1_info = M1.file_info_;
  M2_info = M2.file_info_;
  M1_break_points = M1.break_points_;
  M2_break_points = M2.break_points_;
  M1_chain_stacking = M1.chain_stacking_;
  M2_chain_stacking = M2.chain_stacking_;
  M1_hash_PDBID = M1.hash_PDBID_;
  M2_hash_PDBID = M2.hash_PDBID_;
  M1_num_stacking = M1_num_pairing = M2_num_stacking = M2_num_pairing = 0;
  
  //  copy interactions
  M1_interaction = vector< vector<int> > (M1.record_interaction_.size());
  for(i = 0; i < M1.record_interaction_.size(); ++ i)  {
    M1_interaction[i] = M1.record_interaction_[i];
    if(M1_interaction[i][2] < 0)  {
      ++ M1_num_stacking;
    }  else  {
      ++ M1_num_pairing;
    }
  }
  
  M2_interaction = vector< vector<int> > (M2.record_interaction_.size());
  for(i = 0; i < M2.record_interaction_.size(); ++ i)  {
    M2_interaction[i] = M2.record_interaction_[i];
    if(M2_interaction[i][2] < 0)  {
      ++ M2_num_stacking;
    }  else  {
      ++ M2_num_pairing;
    }
  }
  
  M1_interaction_adjacency = vector< vector<int> > (M1_seq.length());
  for(i = 0; i < M1_seq.length(); ++ i)  {
    M1_interaction_adjacency[i] = vector<int> (M1_seq.length());
    for(j = 0; j < M1_seq.length(); ++ j)  {
      M1_interaction_adjacency[i][j] = M1.interaction_adjacency_[i][j];
    }
  }
  M1_num_interaction = M1_interaction.size();
  M2_interaction_adjacency = vector< vector<int> > (M2_seq.length());
  for(i = 0; i < M2_seq.length(); ++ i)  {
    M2_interaction_adjacency[i] = vector<int> (M2_seq.length());
    for(j = 0; j < M2_seq.length(); ++ j)  {
      M2_interaction_adjacency[i][j] = M2.interaction_adjacency_[i][j];
    }
  }
  M2_num_interaction = M2_interaction.size();  
  return;
}


MotifGraphMatching::MotifGraphMatching(
    const StructuralMotif& M1, 
    const StructuralMotif& M2, 
    const Parameters& input_scoring_function, 
    const float input_mean, const float input_std
)  {
  
  consider_pvalue = true;
  consider_cutoff = false;
  
  unsigned int i, j;
  lower_bound = max_score = -MAX;
  scoring_function = input_scoring_function;
  M1_seq = M1.sequence_;
  M2_seq = M2.sequence_;
  M1_info = M1.file_info_;
  M2_info = M2.file_info_;
  M1_break_points = M1.break_points_;
  M2_break_points = M2.break_points_;
  M1_chain_stacking = M1.chain_stacking_;
  M2_chain_stacking = M2.chain_stacking_;
  M1_hash_PDBID = M1.hash_PDBID_;
  M2_hash_PDBID = M2.hash_PDBID_;
  M1_num_stacking = M1_num_pairing = M2_num_stacking = M2_num_pairing = 0;
  
  //  copy interactions
  M1_interaction = vector< vector<int> > (M1.record_interaction_.size());
  for(i = 0; i < M1.record_interaction_.size(); ++ i)  {
    M1_interaction[i] = M1.record_interaction_[i];
    if(M1_interaction[i][2] < 0)  {
      ++ M1_num_stacking;
    }  else  {
      ++ M1_num_pairing;
    }
  }
  
  M2_interaction = vector< vector<int> > (M2.record_interaction_.size());
  for(i = 0; i < M2.record_interaction_.size(); ++ i)  {
    M2_interaction[i] = M2.record_interaction_[i];
    if(M2_interaction[i][2] < 0)  {
      ++ M2_num_stacking;
    }  else  {
      ++ M2_num_pairing;
    }
  }
  
  M1_interaction_adjacency = vector< vector<int> > (M1_seq.length());
  for(i = 0; i < M1_seq.length(); ++ i)  {
    M1_interaction_adjacency[i] = vector<int> (M1_seq.length());
    for(j = 0; j < M1_seq.length(); ++ j)  {
      M1_interaction_adjacency[i][j] = M1.interaction_adjacency_[i][j];
    }
  }
  M1_num_interaction = M1_interaction.size();
  
  M2_interaction_adjacency = vector< vector<int> > (M2_seq.length());
  for(i = 0; i < M2_seq.length(); ++ i)  {
    M2_interaction_adjacency[i] = vector<int> (M2_seq.length());
    for(j = 0; j < M2_seq.length(); ++ j)  {
      M2_interaction_adjacency[i][j] = M2.interaction_adjacency_[i][j];
    }
  }
  M2_num_interaction = M2_interaction.size();
  
  //  copy simulated data
  mean = input_mean;
  std = input_std;
  return;
}

MotifGraphMatching::MotifGraphMatching(
    const StructuralMotif& M1, 
    const StructuralMotif& M2, 
    const Parameters& input_scoring_function, 
    const float input_mean, const float input_std, 
    const float input_cutoff
)  {
  
  consider_pvalue = true;
  consider_cutoff = true;
  
  unsigned int i, j;
  lower_bound = max_score = -MAX;
  scoring_function = input_scoring_function;
  M1_seq = M1.sequence_;
  M2_seq = M2.sequence_;
  M1_info = M1.file_info_;
  M2_info = M2.file_info_;
  M1_break_points = M1.break_points_;
  M2_break_points = M2.break_points_;
  M1_chain_stacking = M1.chain_stacking_;
  M2_chain_stacking = M2.chain_stacking_;
  M1_hash_PDBID = M1.hash_PDBID_;
  M2_hash_PDBID = M2.hash_PDBID_;
  M1_num_stacking = M1_num_pairing = M2_num_stacking = M2_num_pairing = 0;
  
  //  copy interactions
  M1_interaction = vector< vector<int> > (M1.record_interaction_.size());
  for(i = 0; i < M1.record_interaction_.size(); ++ i)  {
    M1_interaction[i] = M1.record_interaction_[i];
    if(M1_interaction[i][2] < 0)  {
      ++ M1_num_stacking;
    }  else  {
      ++ M1_num_pairing;
    }
  }
  
  M2_interaction = vector< vector<int> > (M2.record_interaction_.size());
  for(i = 0; i < M2.record_interaction_.size(); ++ i)  {
    M2_interaction[i] = M2.record_interaction_[i];
    if(M2_interaction[i][2] < 0)  {
      ++ M2_num_stacking;
    }  else  {
      ++ M2_num_pairing;
    }
  }
  
  M1_interaction_adjacency = vector< vector<int> > (M1_seq.length());
  for(i = 0; i < M1_seq.length(); ++ i)  {
    M1_interaction_adjacency[i] = vector<int> (M1_seq.length());
    for(j = 0; j < M1_seq.length(); ++ j)  {
      M1_interaction_adjacency[i][j] = M1.interaction_adjacency_[i][j];
    }
  }
  M1_num_interaction = M1_interaction.size();
  
  M2_interaction_adjacency = vector< vector<int> > (M2_seq.length());
  for(i = 0; i < M2_seq.length(); ++ i)  {
    M2_interaction_adjacency[i] = vector<int> (M2_seq.length());
    for(j = 0; j < M2_seq.length(); ++ j)  {
      M2_interaction_adjacency[i][j] = M2.interaction_adjacency_[i][j];
    }
  }
  M2_num_interaction = M2_interaction.size();
  
  //  copy simulated data
  mean = input_mean;
  std = input_std;
  pvalue_cutoff = input_cutoff;
  return;
}


MotifGraphMatching::~MotifGraphMatching()  {
  return;
}


void MotifGraphMatching::align(void)  {
  //  sort the interactions
  sort_interaction();
  //  compute the interaction compatibility
  compute_interaction_compatibility();
  //  compute the nucleotide deletion penalties
  compute_nuc_del();
  
  /*
  for(int i = 0; i < M1_seq.length(); ++ i) {
    for(int j = 0; j < M1_seq.length(); ++ j) {
      cout << M1_nuc_del_penalty[i][j] << "\t";
    }
    cout << endl;
  }
  cout << endl;
  
  for(int i = 0; i < M2_seq.length(); ++ i) {
    for(int j = 0; j < M2_seq.length(); ++ j) {
      cout << M2_nuc_del_penalty[i][j] << "\t";
    }
    cout << endl;
  }
  */
  
  //  establish an lower bound first
  if(consider_pvalue && consider_cutoff)  {
    lower_bound = pvalue_lowerbound();
  }
    
  //std::string alignment1, alignment2;
  //loop_align(6, 10, alignment1, alignment2, 11, 16);
  //return;
    
  //find_maximal_clique_heuristic();
  //  compute the optimal alignment
  find_maximal_clique_optimal();
    
  //maximal_clique.push_back(22);
  //maximal_clique.push_back(24);
  //maximal_clique.push_back(48);
  //maximal_clique.push_back(24);
  //maximal_clique.push_back(39);
  //maximal_clique.push_back(54);
  
  //interpret_maximal_clique();
  
  //std::cout << "******************************" << std::endl;
  max_score = compute_alignment_score(maximal_clique);
  //std::cout << "check score: " << max_score << std::endl;
  //  output the alignment
  pvalue = chebyshev_inequality();
  //print_alignment(maximal_clique);
  return;
}

void MotifGraphMatching::align_heuristic()  {
  //cout << "heuristic mode enabled" << endl;
  //  sort the interactions
  sort_interaction();
  //  compute the interaction compatibility
  compute_interaction_compatibility();
  //  compute the nucleotide deletion penalties
  compute_nuc_del();
  //  establish an lower bound first
  find_maximal_clique_heuristic();
  pvalue = chebyshev_inequality();
  //print_alignment(maximal_clique);
  return;
}
  
void MotifGraphMatching::align_simulation()  {
  //  sort the interactions
  sort_interaction();
  //  compute the interaction compatibility
  compute_interaction_compatibility();
  //  compute the nucleotide deletion penalties
  compute_nuc_del();
  //  establish an lower bound first
  find_maximal_clique_heuristic();
}

void MotifGraphMatching::sort_interaction(void)  {
  int i, j;
  vector<int> temp (6);
  for(i = 0; i < M1_num_interaction - 1; ++ i)  {
    for(j = i + 1; j < M1_num_interaction; ++ j)  {
      if(M1_interaction[i][1] > M1_interaction[j][1] || 
        (M1_interaction[i][1] == M1_interaction[j][1] && M1_interaction[i][0] < M1_interaction[j][0]))  {
        temp = M1_interaction[i];
        M1_interaction[i] = M1_interaction[j];
        M1_interaction[j] = temp;
      }
    }
  }
  
  /*
  for(i = 0; i < M1_num_interaction; ++ i)  {
    cout << M1_interaction[i][0] << "  " << M1_interaction[i][1] << "  " << M1_interaction[i][2] << endl;
  }
  */
    
  for(i = 0; i < M2_num_interaction - 1; ++ i)  {
    for(j = i + 1; j < M2_num_interaction; ++ j)  {
      if(M2_interaction[i][1] > M2_interaction[j][1] || 
        (M2_interaction[i][1] == M2_interaction[j][1] && M2_interaction[i][0] < M2_interaction[j][0]))  {
        temp = M2_interaction[i];
        M2_interaction[i] = M2_interaction[j];
        M2_interaction[j] = temp;
      }
    }
  }
  /*
  for(i = 0; i < M2_num_interaction; ++ i)  {
    cout << M2_interaction[i][0] << "  " << M2_interaction[i][1] << "  " << M2_interaction[i][2] << endl;
  }
  cout << "Interaction in motif 2" << endl;
  */
  return;
}

void MotifGraphMatching::compute_interaction_compatibility(void)  {
  
  int i, j, k, l;
  
  //  allocate space
  interaction_compatibility = vector< vector<int> > (M1_num_interaction * M2_num_interaction);
  for(i = 0; i < M1_num_interaction * M2_num_interaction; ++ i)  {
    interaction_compatibility[i] = vector<int> (M1_num_interaction * M2_num_interaction);
    for(j = 0; j < M1_num_interaction * M2_num_interaction; ++ j)  {
      interaction_compatibility[i][j] = 0;
    }
  }
  
  //  construct the compatibility graph
  for(i = 0; i < M1_num_interaction - 1; ++ i)  {
    for(j = 0; j < M2_num_interaction - 1; ++ j)  {
    
      int index_1 = i * M2_num_interaction + j;
      //  to ensure that a basepair cannot match with a stacking
      if((M1_interaction[i][2] > 0 && M2_interaction[j][2] < 0) || 
        (M1_interaction[i][2] < 0 && M2_interaction[j][2] > 0))  {
        continue;
      }
      
      for(k = i + 1; k < M1_num_interaction; ++ k)  {
        for(l = j + 1; l < M2_num_interaction; ++ l)  {
          
          int index_2 = k * M2_num_interaction + l;
          //  to ensure that a basepair cannot match with a stacking
          if((M1_interaction[k][2] > 0 && M2_interaction[l][2] < 0) || 
            (M1_interaction[k][2] < 0 && M2_interaction[l][2] > 0))  {
            continue;
          }
          
          //  check whether the matching of (i, j) is compatible with the matching (k, l)
          if(  (M1_interaction[k][1] < M1_interaction[i][0] && M2_interaction[l][1] < M2_interaction[j][0]) ||  //  case 1
            (M1_interaction[k][0] > M1_interaction[i][1] && M2_interaction[l][0] > M2_interaction[j][1]) ||  //  case 2
            (M1_interaction[k][1] == M1_interaction[i][0] && M2_interaction[l][1] == M2_interaction[j][0]) ||  //  case 3
            (M1_interaction[k][0] == M1_interaction[i][1] && M2_interaction[l][0] == M2_interaction[j][1]) ||  //  case 4
            (M1_interaction[k][0] == M1_interaction[i][0] && M1_interaction[k][1] < M1_interaction[i][1] && 
             M2_interaction[l][0] == M2_interaction[j][0] && M2_interaction[l][1] < M2_interaction[j][1]) ||  //  case 5
            (M1_interaction[k][0] == M1_interaction[i][0] && M1_interaction[k][1] > M1_interaction[i][1] && 
             M2_interaction[l][0] == M2_interaction[j][0] && M2_interaction[l][1] > M2_interaction[j][1]) ||  //  case 6
            (M1_interaction[k][1] == M1_interaction[i][1] && M1_interaction[k][0] < M1_interaction[i][0] && 
             M2_interaction[l][1] == M2_interaction[j][1] && M2_interaction[l][0] < M2_interaction[j][0]) ||  //  case 7 
            (M1_interaction[k][1] == M1_interaction[i][1] && M1_interaction[k][0] > M1_interaction[i][0] && 
             M2_interaction[l][1] == M2_interaction[j][1] && M2_interaction[l][0] > M2_interaction[j][0]) ||  //  case 8
            (M1_interaction[k][0] > M1_interaction[i][0] && M1_interaction[k][1] < M1_interaction[i][1] && 
             M2_interaction[l][0] > M2_interaction[j][0] && M2_interaction[l][1] < M2_interaction[j][1]) ||  //  case 9
            (M1_interaction[k][0] < M1_interaction[i][0] && M1_interaction[k][1] > M1_interaction[i][1] && 
             M2_interaction[l][0] < M2_interaction[j][0] && M2_interaction[l][1] > M2_interaction[j][1]) ||  //  case 10
            (M1_interaction[k][0] < M1_interaction[i][0] && M1_interaction[k][1] > M1_interaction[i][0] && 
             M1_interaction[k][1] < M1_interaction[i][1] && M2_interaction[l][0] < M2_interaction[j][0] && 
             M2_interaction[l][1] > M2_interaction[j][0] && M2_interaction[l][1] < M2_interaction[j][1]) ||  //  case 11
            (M1_interaction[k][0] > M1_interaction[i][0] && M1_interaction[k][0] < M1_interaction[i][1] && 
             M1_interaction[k][1] > M1_interaction[i][1] && M2_interaction[l][0] > M2_interaction[j][0] && 
             M2_interaction[l][0] < M2_interaction[j][1] && M2_interaction[l][1] > M2_interaction[j][1])  //  case 12
            )  {
            //  if either of these case is satisfied, then the matching (i, j) and (k, l) are compatible
            interaction_compatibility[index_1][index_2] = 1;
            interaction_compatibility[index_2][index_1] = 1;

          }          

        }
      }
    }
  }
  
  //  verify that all compatible interactions are of the same type
  /*
  for(i = 0; i < M1_num_interaction * M2_num_interaction; ++ i)  {
    for(j = 0; j < M1_num_interaction * M2_num_interaction; ++ j)  {
      cout << i << "  " << j << " " << interaction_compatibility[i][j] << endl;
    }
  }
  */
  return;
}

void MotifGraphMatching::interpret_interaction_compatibility(void)  {
  int i, j;
  int idx_i, idx_j, idx_k, idx_l;
  for(i = 0; i < M1_num_interaction * M2_num_interaction - 1; ++ i)  {
    for(j = i + 1; j < M1_num_interaction * M2_num_interaction; ++ j)  {
      if(interaction_compatibility[i][j] == 1)  {
        idx_i = i / M2_num_interaction;
        idx_j = i % M2_num_interaction;
        idx_k = j / M2_num_interaction;
        idx_l = j % M2_num_interaction;
        
        cout << "*******" << endl;
        cout << "The matching of (" << M1_interaction[idx_i][0] << "," << M1_interaction[idx_i][1] << "," << M1_interaction[idx_i][2] << ") and (" << M2_interaction[idx_j][0] << "," << M2_interaction[idx_j][1] << "," << M2_interaction[idx_j][2] << ")" << endl;
        cout << "is compatible with the matching of (" << M1_interaction[idx_k][0] << "," << M1_interaction[idx_k][1] << "," << M1_interaction[idx_k][2] << ") and (" << M2_interaction[idx_l][0] << "," << M2_interaction[idx_l][1] << "," << M2_interaction[idx_l][2] << ")" << endl;
      }
    }
  }
  return;
}

void MotifGraphMatching::interpret_maximal_clique(void)  {
  int idx_i, idx_j;
  for(int i = 0; i < (int) maximal_clique.size(); ++ i)  {
    idx_i = maximal_clique[i] / M2_num_interaction;
    idx_j = maximal_clique[i] % M2_num_interaction;
    cout << "*******: " << maximal_clique[i] << endl;
    cout << "The matching of (" << M1_interaction[idx_i][0] << "," << M1_interaction[idx_i][1] << "," << M1_interaction[idx_i][2] << ") and (" << M2_interaction[idx_j][0] << "," << M2_interaction[idx_j][1] << "," << M2_interaction[idx_j][2] << ")" << endl;
  }
  return;
}

void MotifGraphMatching::interpret_candidates(const vector<bool>& candidates)  {
  int idx_i, idx_j;
  for(int i = 0; i < (int) candidates.size(); ++ i)  {
    idx_i = candidates[i] / M2_num_interaction;
    idx_j = candidates[i] % M2_num_interaction;
    cout << "*******" << endl;
    cout << "The matching of (" << M1_interaction[idx_i][0] << "," << M1_interaction[idx_i][1] << "," << M1_interaction[idx_i][2] << ") and (" << M2_interaction[idx_j][0] << "," << M2_interaction[idx_j][1] << "," << M2_interaction[idx_j][2] << ")" << endl;
  }
  return;
}

int MotifGraphMatching::pvalue_lowerbound(void)  {
  
  if(consider_pvalue && consider_cutoff)  {  
    return mean + std * sqrt(1 / pvalue_cutoff);
  }  else  {
    return 0;
  }
}

void MotifGraphMatching::find_maximal_clique_heuristic(void)  {

  int i, j;
  //  initialize the array maximal_clique
  vector<int> heuristic_clique;
  vector<int> candidates;
  for(i = 0; i < M1_num_interaction; ++ i) {
    for(j = 0; j < M2_num_interaction; ++ j) {
      if(M1_interaction[i][2] * M2_interaction[j][2] >= 0)  { 
        // meaning that they are either both pairing. or both stacking, or one pairing/stacking matched with an H-bond
        candidates.push_back(i * M2_num_interaction + j);
      }
    }
  }
  
  do  {
    /*
    cout << "heristic clique" << endl;
    for(auto it_i = heuristic_clique.begin(); it_i != heuristic_clique.end(); ++ it_i) {
      cout << *it_i << "  ";
    }
    cout << endl;
    
    cout << "candidate clique" << endl;
    for(auto it_i = candidates.begin(); it_i != candidates.end(); ++ it_i) {
      cout << *it_i << "  ";
    }
    cout << endl;
    */
    int max_degree = -1;
    int max_degree_index = -1;
    //  identify the mostly connected vertex
    for(i = 0; i < (int) candidates.size(); ++ i)  {
      int degree = 0;
      for(j = 0; j < (int) candidates.size(); ++ j)  {
        if(i != j && interaction_compatibility[candidates[i]][candidates[j]] == 1)  {
          ++ degree;
        }      
      }
      //  record the maximum degree
      if(degree > max_degree)  {
        max_degree = degree;
        max_degree_index = i;
      }
    }
    //  take the vertex with maximal degree as a part of the clique
    heuristic_clique.push_back(candidates[max_degree_index]);
    //  remove vertex that are not connected to the vertex with maximum degree
    ReduceSizeUnsorted<int>(candidates, max_degree_index);
    list<int> index_to_remove;
    for(i = 0; i < (int) candidates.size(); ++ i)  {
      if(interaction_compatibility[heuristic_clique.back()][candidates[i]] == 0)  {
        index_to_remove.push_back(i);
      }
    }
    for(auto it = index_to_remove.rbegin(); it != index_to_remove.rend(); ++ it) {
      ReduceSizeUnsorted<int>(candidates, *it);
    }
  }  while(candidates.size() > 0);
  //  estimate the lower bound that can be achieved by the heuristic solution
  int current_score = compute_alignment_score(heuristic_clique);
  if(current_score > max_score)  {
    max_score = current_score;
    maximal_clique = heuristic_clique;
  }  
  if(max_score > lower_bound)  {
    lower_bound = max_score;
  }
  return;
}  

void MotifGraphMatching::find_maximal_clique_optimal(void)  {

  vector<int> initial_clique;
  vector<int> initial_candidates;
  //  remove incompatible base interactions from candidates
  int i, j;
  /*
  for(i = 0; i < initial_candidates.size(); ++ i)  {
    int idx_i = i / M2_num_interaction;
    int idx_j = i % M2_num_interaction;
    if(M1_interaction[idx_i][2] * M2_interaction[idx_j][2] < 0)  {
      initial_candidates[i] = false;
      -- initial_candidate_size;
    }
  }
  */
  for(i = 0; i < M1_num_interaction; ++ i) {
    for(j = 0; j < M2_num_interaction; ++ j) {
      if(M1_interaction[i][2] * M2_interaction[j][2] >= 0)  { 
        // meaning that they are either both pairing. or both stacking, or one pairing/stacking matched with an H-bond
        initial_candidates.push_back(i * M2_num_interaction + j);
      }
    }
  }
  
  //cout << "!!!!!!  " << initial_clique_size << "  " << initial_candidate_size << endl;
  extend_clique_recursive(initial_clique, initial_candidates);
  //cout << "******last call of compute_alignment_score******" << endl;
  //int temp_score = compute_alignment_score(maximal_clique);
  //cout << temp_score << endl;
  return;
}

void MotifGraphMatching::extend_clique_recursive(const vector<int>& current_clique, const vector<int>& candidates)  {
  
  /*
  cout << "****************begin recursion*************************" << endl;
  cout << "current clique" << endl;
  for(auto it_i = current_clique.begin(); it_i != current_clique.end(); ++ it_i) {
    cout << *it_i << "  ";
  }
  cout << endl;
    
  cout << "candidate clique" << endl;
  for(auto it_i = candidates.begin(); it_i != candidates.end(); ++ it_i) {
    cout << *it_i << "  ";
  }
  cout << endl;
  */
  
  int alignment_score = compute_alignment_score(current_clique);
  if(alignment_score > max_score)  {
    max_score = alignment_score;
    maximal_clique = current_clique;
  }
  if(max_score > lower_bound)  {
    lower_bound = max_score;
  }  
  //  check boundary condition
  if(candidates.empty())  {
    //cout << "return condition 1" << endl;
    return;
  }  
  int alignment_score_upper_bound = compute_alignment_score_upper_bound(current_clique, candidates);  
  if(alignment_score_upper_bound <= lower_bound)  {
    //cout << "return condition 2" << endl;
    return;
  }  
  
  /*
  cout << "current clique" << endl;
  for(auto it_i = current_clique.begin(); it_i != current_clique.end(); ++ it_i) {
    cout << *it_i << "  ";
  }
  cout << endl;
    
  cout << "candidate clique" << endl;
  for(auto it_i = candidates.begin(); it_i != candidates.end(); ++ it_i) {
    cout << *it_i << "  ";
  }
  cout << endl;
  */
  
  int i, j;
  
  for(i = 0; i < (int) candidates.size(); ++ i)  {
    vector<int> next_clique = current_clique;
    next_clique.push_back(candidates[i]);
    //  pruning candidates
    vector<int> next_candidates;
    for(j = 0; j < (int) candidates.size(); ++ j)  {
      if(j != i && interaction_compatibility[candidates[i]][candidates[j]] == 1)  {
        next_candidates.push_back(candidates[j]);
      }
    }
    /*
    cout << "next current clique" << endl;
    for(auto it_i = next_clique.begin(); it_i != next_clique.end(); ++ it_i) {
      cout << *it_i << "  ";
    }
    cout << endl;
    
    cout << "next candidate clique" << endl;
    for(auto it_i = next_candidates.begin(); it_i != next_candidates.end(); ++ it_i) {
      cout << *it_i << "  ";
    }
    cout << endl;
    */
    extend_clique_recursive(next_clique, next_candidates);
  }
  return;
}

void MotifGraphMatching::compute_nuc_del(void)  {
  unsigned int i, j, k;
  //  initialize space
  M1_nuc_del_penalty = vector< vector<int> > (M1_seq.length());
  M1_effective_stacking = vector< vector<int> > (M1_seq.length());
  for(i = 0; i < M1_seq.length(); ++ i)  {
    M1_nuc_del_penalty[i] = vector<int> (M1_seq.length(), 0);
    M1_effective_stacking[i] = vector<int> (M1_seq.length(), 0);
  }
  for(i = 0; i < M1_seq.length() - 1; ++ i)  {
    M1_effective_stacking[i][i + 1] = M1_chain_stacking[i];
    M1_effective_stacking[i + 1][i] = M1_chain_stacking[i];
  }
  for(i = 0; i < M1_seq.length() - 2; ++ i)  {
    for(j = i + 2; j < M1_seq.length(); ++ j)  {
      bool between_break_point = false;
      //  check for breakpoints
      for(k = 0; k < M1_break_points.size(); ++ k)  {
        if(i < (unsigned int) M1_break_points[k] && j >= (unsigned int) M1_break_points[k])  {
          //  penalty equals to 0, no need to update the score
          between_break_point = true;
          break;
        }
      }
      if(between_break_point)  {
        continue;
      } else  {
        //  evaluate the score for deleting nucleotides between i and j
        M1_nuc_del_penalty[i][j] += scoring_function.gap_opening + scoring_function.gap_extension * (j - i - 1);
        M1_nuc_del_penalty[j][i] = M1_nuc_del_penalty[i][j];  
        //  evaluate effective stacking. see structural_motif::hash_stacking for indication of numbers
        if(  (M1_chain_stacking[i] == -1 || M1_chain_stacking[i] == -4) &&
          (M1_chain_stacking[j] == -1 || M1_chain_stacking[j] == -4))  {
          M1_effective_stacking[i][j] = -1;
        }  else if((M1_chain_stacking[i] == -2 || M1_chain_stacking[i] == -3) && 
          (M1_chain_stacking[j] == -2 || M1_chain_stacking[j] == -3))  {
          M1_effective_stacking[i][j] = -2;
        }  else if((M1_chain_stacking[i] == -2 || M1_chain_stacking[i] == -3) && 
          (M1_chain_stacking[j] == -1 || M1_chain_stacking[j] == -4))  {
          M1_effective_stacking[i][j] = -3;
        }  else if((M1_chain_stacking[i] == -1 || M1_chain_stacking[i] == -4) && 
          (M1_chain_stacking[j] == -2 || M1_chain_stacking[j] == -3))  {
          M1_effective_stacking[i][j] = -4;
        }
        M1_effective_stacking[j][i] = M1_effective_stacking[i][j];
      }
    }
  }
  
  
  //  initialize space
  M2_nuc_del_penalty = vector< vector<int> > (M2_seq.length());
  M2_effective_stacking = vector< vector<int> > (M2_seq.length());
  for(i = 0; i < M2_seq.length(); ++ i)  {
    M2_nuc_del_penalty[i] = vector<int> (M2_seq.length(), 0);
    M2_effective_stacking[i] = vector<int> (M2_seq.length(), 0);
  }
  for(i = 0; i < M2_seq.length() - 1; ++ i)  {
    M2_effective_stacking[i][i + 1] = M2_chain_stacking[i];
    M2_effective_stacking[i + 1][i] = M2_chain_stacking[i];
  }
  for(i = 0; i < M2_seq.length() - 2; ++ i)  {
    for(j = i + 2; j < M2_seq.length(); ++ j)  {
      bool between_break_point = false;
      //  check for breakpoints
      for(k = 0; k < M2_break_points.size(); ++ k)  {
        if(i < (unsigned int) M2_break_points[k] && j >= (unsigned int) M2_break_points[k])  {
          //  penalty equals to 0, no need to update the score
          between_break_point = true;
          break;
        }
      }
      if(between_break_point)  {
        continue;
      } else  {
        //  evaluate the score for deleting nucleotides between i and j
        //  need to apply more sophisticated GAP opening AND EXTENSION scoring function
        M2_nuc_del_penalty[i][j] += scoring_function.gap_opening + scoring_function.gap_extension * (j - i - 1);
        M2_nuc_del_penalty[j][i] = M2_nuc_del_penalty[i][j];
        //  evaluate effective stacking. see structural_motif::hash_stacking for indication of numbers
        if(  (M2_chain_stacking[i] == -1 || M2_chain_stacking[i] == -4) &&
          (M2_chain_stacking[j] == -1 || M2_chain_stacking[j] == -4))  {
          M2_effective_stacking[i][j] = -1;
        }  else if((M2_chain_stacking[i] == -2 || M2_chain_stacking[i] == -3) && 
          (M2_chain_stacking[j] == -2 || M2_chain_stacking[j] == -3))  {
          M2_effective_stacking[i][j] = -2;
        }  else if((M2_chain_stacking[i] == -2 || M2_chain_stacking[i] == -3) && 
          (M2_chain_stacking[j] == -1 || M2_chain_stacking[j] == -4))  {
          M2_effective_stacking[i][j] = -3;
        }  else if((M2_chain_stacking[i] == -1 || M2_chain_stacking[i] == -4) && 
          (M2_chain_stacking[j] == -2 || M2_chain_stacking[j] == -3))  {
          M2_effective_stacking[i][j] = -4;
        }
        M2_effective_stacking[j][i] = M2_effective_stacking[i][j];
      }
    }
  }
  
  return;
}  

int MotifGraphMatching::loop_align(void)  {
  //  compute the semi-local alignment of the two sequences
  int i, j, k, l;
  //  if either or both segment is empty
  if(M1_seq.length() <= 0)  {
    if(M2_seq.length() <= 0)  {
      return 0;
    }  else  {
      return scoring_function.weight_sequence * (scoring_function.gap_opening + M2_seq.length() * scoring_function.gap_extension);
    }
  }  else if(M2_seq.length() <= 0)  {
    return scoring_function.weight_sequence * (scoring_function.gap_opening + M1_seq.length() * scoring_function.gap_extension);
  }
  
  //  allocate space
  vector< vector<int> > DP_table(M1_seq.length());
  for(i = 0; i < (int) M1_seq.length(); ++ i)  {
    DP_table[i] = vector<int> (M2_seq.length());
  }
  
  //  boundary case
  for(i = 0; i < (int) M1_seq.length(); ++ i)  {
    DP_table[i][0] = scoring_function.weight_sequence * scoring_function.match_nucleotide(M1_seq[i], M2_seq[0]);
  }
  for(j = 1; j < (int) M2_seq.length(); ++ j)  {
    DP_table[0][j] = scoring_function.weight_sequence * scoring_function.match_nucleotide(M1_seq[0], M2_seq[j]);
  }
  
  //  general cases
  for(i = 1; i < (int) M1_seq.length(); ++ i)  {
    for(j = 1; j < (int) M2_seq.length(); ++ j)  {
      int max_score = -MAX, temp_score;
      //int max_i, max_j;
      //  go over all entries in 0...i - 1, j - 1
      for(k = 0; k < i; ++ k)  {
        //  the previous alignment ends at k, j - 1
        temp_score = DP_table[k][j - 1];
        temp_score += scoring_function.weight_adjacent_stacking 
          * scoring_function.match_stacking_basepair(M1_effective_stacking[k][i], M2_effective_stacking[j - 1][j]);
        temp_score += scoring_function.weight_sequence * M1_nuc_del_penalty[k][i];
        temp_score += scoring_function.weight_sequence * scoring_function.match_nucleotide(M1_seq[i], M2_seq[j]);
        if(temp_score > max_score)  {
          max_score = temp_score;
          //max_i = k;
          //max_j = j - 1;
        }
      }
      //  go over all entries in i - 1, 0...j - 1
      for(k = 0; k < j; ++ k)  {
        //  the previous alignment ends at k, j - 1
        temp_score = DP_table[i - 1][k];
        temp_score += scoring_function.weight_adjacent_stacking 
          * scoring_function.match_stacking_basepair(M1_effective_stacking[i - 1][i], M2_effective_stacking[k][j]);
        temp_score += scoring_function.weight_sequence * M2_nuc_del_penalty[k][j];
        temp_score += scoring_function.weight_sequence * scoring_function.match_nucleotide(M1_seq[i], M2_seq[j]);
        if(temp_score > max_score)  {
          max_score = temp_score;
          //max_i = i - 1;
          //max_j = k;
        }
      }
      //  go over all 0...breakpoint_i, 0...breakpoint_j
      //  identify the last breakpoints for M1 and M2, respectively
      int M1_last = -1, M2_last = -1;
      for(k = 0; k < (int) M1_break_points.size(); ++ k)  {
        if(M1_break_points[k] > M1_last && M1_break_points[k] <= i)  {
          M1_last = M1_break_points[k];
        }
      }
      for(k = 0; k < (int) M2_break_points.size(); ++ k)  {
        if(M2_break_points[k] > M2_last && M2_break_points[k] <= j)  {
          M2_last = M2_break_points[k];
        }
      }
      for(k = 0; k < M1_last; ++ k)  {
        for(l = 0; l < M2_last; ++ l)  {
          //  the previous alignment ends at k, j - 1
          temp_score = DP_table[k][l];
          temp_score += scoring_function.weight_adjacent_stacking 
            * scoring_function.match_stacking_basepair(M1_effective_stacking[k][i], M2_effective_stacking[l][j]);
          temp_score += scoring_function.weight_sequence * scoring_function.match_nucleotide(M1_seq[i], M2_seq[j]);
          if(temp_score > max_score)  {
            max_score = temp_score;
            //max_i = k;
            //max_j = l;
          }
        }
      }
      //  assign score
      //  local alignment, taking the matching of i, j as start, need more sophisticated NUC SUBSTITUTION function
      int initial_match_score = scoring_function.weight_sequence * scoring_function.match_nucleotide(M1_seq[i], M2_seq[j]);
      if(max_score >  initial_match_score)  {
        DP_table[i][j] = max_score;
      }  else  {
        DP_table[i][j] = initial_match_score;
      }
    }
  }
  
  //  finds the maximum alignment score from the entire matrix
  int maximum_align_score = 0;
  for(i = 0; i < (int) M1_seq.length(); ++ i)  {
    for(j = 0; j < (int) M2_seq.length(); ++ j)  {
      if(DP_table[i][j] > maximum_align_score)  {
        maximum_align_score = DP_table[i][j];
      }
    }
  }
  /*
  for(i = 0; i < M1_seq.length(); ++ i)  {
    for(j = 0; j < M2_seq.length(); ++ j)  {
      cout << DP_table[i][j] << "  ";
    }
    cout << endl;
  }
  */
  return maximum_align_score;
}

int MotifGraphMatching::loop_align(int b1, int b2, string& M1_alignment, string& M2_alignment, int e1, int e2)  {
  //cout << "****************" << endl;
  //cout << b1 << " " << e1 << "  " << b2 << "  " << e2 << endl;
  //cout << "sequences" << endl;
  //cout << M1_seq << endl;
  //cout << M2_seq << endl;
  // check if a break point is present here, if yes then return 0
  for(auto it = M1_break_points.begin(); it != M1_break_points.end(); ++ it) {
    if(*it > b1 && *it <= e1)  {
      //cout << "break 1: " << *it << " " << b1 << "  " << e1 << endl;
      M1_alignment = M2_alignment = "~";
      return 0;
    }
  }
  for(auto it = M2_break_points.begin(); it != M2_break_points.end(); ++ it) {
    if(*it > b2 && *it <= e2)  {
      //cout << "break 2: " << *it << " " << b2 << "  " << e2 << endl;
      M1_alignment = M2_alignment = "~";
      return 0;
    }
  }
  
  //  compute the global alignment of the two sequences
  int i, j, k;
  //  clear the alignment holder
  M1_alignment = string("");
  M2_alignment = string("");
  //  compute effective length
  int len_seg_1 = e1 - b1 - 1;
  int len_seg_2 = e2 - b2 - 1;
  //cout << "lengths: " << len_seg_1 << " " << len_seg_2 << endl;
  
  //  if either or both segment is empty
  
  
  if(len_seg_1 <= 0)  {
    if(len_seg_2 <= 0)  {
      M1_alignment = "";
      M2_alignment = "";
      int return_score =  scoring_function.weight_adjacent_stacking 
        * scoring_function.match_stacking_basepair(M1_effective_stacking[b1][e1], M2_effective_stacking[b2][e2]);
      //cout << "in loop align: " << return_score << endl;
      return return_score;
    }  else  {
      //  construct alignment
      M1_alignment += string(len_seg_2, '-');
      M2_alignment += M2_seq.substr(b2 + 1, len_seg_2);
      int return_score = scoring_function.weight_adjacent_stacking 
        * scoring_function.match_stacking_basepair(M1_effective_stacking[b1][e1], M2_effective_stacking[b2][e2]);
      return_score += scoring_function.weight_sequence * M2_nuc_del_penalty[b2][e2];
      
      //cout << "in loop align 2: " << return_score << "  " << M2_nuc_del_penalty[b2][e2] << endl;
      return return_score;
    }
  }  else if(len_seg_2 <= 0)  {
    //  construct alignment
    M1_alignment += M1_seq.substr(b1 + 1, len_seg_1);
    M2_alignment += string(len_seg_1, '-');
    int return_score = scoring_function.weight_adjacent_stacking 
      * scoring_function.match_stacking_basepair(M1_effective_stacking[b1][e1], M2_effective_stacking[b2][e2]);
    return_score += scoring_function.weight_sequence * M1_nuc_del_penalty[b1][e1];
    //cout << "in loop align 1: " << return_score << "  " << M1_nuc_del_penalty[b1][e1] << endl;
    return return_score;
  }
  
  //  allocate space. note that we need to include e1-e2 matching into the matrix, but not b1-b2 matching
  vector< vector<int> > DP_table(len_seg_1 + 2);
  vector< vector<int> > record_i(len_seg_1 + 2);
  vector< vector<int> > record_j(len_seg_1 + 2);
  for(i = 0; i < len_seg_1 + 2; ++ i)  {
    DP_table[i] = vector<int> (len_seg_2 + 2);
    record_i[i] = vector<int> (len_seg_2 + 2);
    record_j[i] = vector<int> (len_seg_2 + 2);
  }
  //cout << "Here!" << endl;
  //  boundary case
  DP_table[0][0] = 0;
  record_i[0][0] = record_j[0][0] = -1;
  for(i = 1; i < len_seg_1 + 2; ++ i)  {
    DP_table[i][0] = scoring_function.weight_sequence * (scoring_function.gap_opening + i * scoring_function.gap_extension);
    record_i[i][0] = record_j[i][0] = 0;
  }
  for(j = 1; j < len_seg_2 + 2; ++ j)  {
    DP_table[0][j] = scoring_function.weight_sequence * (scoring_function.gap_opening + j * scoring_function.gap_extension);
    record_i[0][j] = record_j[0][j] = 0;
  }
  
  
  //  general cases
  for(i = 1; i < len_seg_1 + 2; ++ i)  {
    for(j = 1; j < len_seg_2 + 2; ++ j)  {
      //cout << "index: " << i << " " << j << endl;
      int max_score = -MAX, temp_score;
      int max_i = i - 1, max_j = j - 1;
      //  go over all entries in 0...i - 1, j - 1
      for(k = 0; k < i; ++ k)  {
        //  the previous alignment ends at k, j - 1
        temp_score = DP_table[k][j - 1];
        temp_score += scoring_function.weight_adjacent_stacking
          * scoring_function.match_stacking_basepair(M1_effective_stacking[b1 + k][b1 + i], M2_effective_stacking[b2 + j - 1][b2 + j]);
        temp_score += scoring_function.weight_sequence * M1_nuc_del_penalty[b1 + k][b1 + i];    
        temp_score += scoring_function.weight_sequence * scoring_function.match_nucleotide(M1_seq[b1 + i], M2_seq[b2 + j]);
        //cout << "k: " << k << " " << temp_score << endl;
        if(temp_score > max_score)  {
          max_score = temp_score;
          max_i = k;
          max_j = j - 1;
        }
      }
      //  go over all entries in i - 1, 0...j - 1
      for(k = 0; k < j; ++ k)  {
        //  the previous alignment ends at k, j - 1
        
        temp_score = DP_table[i - 1][k];
        temp_score += scoring_function.weight_adjacent_stacking 
          * scoring_function.match_stacking_basepair(M1_effective_stacking[b1 + i - 1][b1 + i], M2_effective_stacking[b2 + k][b2 + j]);
        temp_score += scoring_function.weight_sequence * M2_nuc_del_penalty[b2 + k][b2 + j];
        temp_score += scoring_function.weight_sequence * scoring_function.match_nucleotide(M1_seq[b1 + i], M2_seq[b2 + j]);
        
        if(temp_score > max_score)  {
          max_score = temp_score;
          max_i = i - 1;
          max_j = k;
        }
      }
      //  assign score
      DP_table[i][j] = max_score;
      record_i[i][j] = max_i;
      record_j[i][j] = max_j;
      //cout << "max index: " << max_i << " " << max_j << endl;
    }
  }
  
  /*
  for(i = 0; i < len_seg_1 + 2; ++ i)  {
    for(j = 0; j < len_seg_2 + 2; ++ j)  {
      cout << DP_table[i][j] << "\t";
    }
    cout << endl;
  }
  cout << endl;
  
  for(i = 0; i < len_seg_1 + 2; ++ i)  {
    for(j = 0; j < len_seg_2 + 2; ++ j)  {
      cout << record_i[i][j] << "\t";
    }
    cout << endl;
  }
  cout << endl;
  
  for(i = 0; i < len_seg_1 + 2; ++ i)  {
    for(j = 0; j < len_seg_2 + 2; ++ j)  {
      cout << record_j[i][j] << "\t";
    }
    cout << endl;
  }
  cout << endl;
  */
  //cout << "before trace back" << endl;
  // TODO: examine the last row and column and apply stacking matching and gap penalty to e1 and e2
  
  //  traceback and recover alignment
  
  int tb_i = e1 - b1, tb_j = e2 - b2;
  int hold_tb_i, hold_tb_j;
  while(tb_i > 0 || tb_j > 0)  {
    //  add tb_i and tb_j to the existing alignment if neither of them is on boundary
    if(tb_i > 0 && tb_j > 0)  {
      M1_alignment += M1_seq.substr(b1 + tb_i, 1);
      M2_alignment += M2_seq.substr(b2 + tb_j, 1);
    }
    // jump back
    hold_tb_i = record_i[tb_i][tb_j];
    hold_tb_j = record_j[tb_i][tb_j];
    //cout << tb_i << " " << tb_j << "  " << hold_tb_i << " " << hold_tb_j << endl;
    if(tb_i - hold_tb_i == 1 && tb_j - hold_tb_j == 1)  {
      // do nothing
      ;
    } else if((tb_i - hold_tb_i > 1 && tb_j - hold_tb_j == 1) || (tb_i - hold_tb_i > 0 && hold_tb_j == 0))  {
      string seg = M1_seq.substr(b1 + hold_tb_i + 1, tb_i - hold_tb_i);
      //cout << "segA: " << seg << endl;
      M1_alignment += string(seg.rbegin(), seg.rend());
      M2_alignment += string(tb_i - hold_tb_i, '-');
    } else if((tb_i - hold_tb_i == 1 && tb_j - hold_tb_j > 1) || (hold_tb_i == 0 && tb_j - hold_tb_j > 0)) {
      string seg = M2_seq.substr(b2 + hold_tb_j + 1, tb_j - hold_tb_j);
      //cout << "segB: " << seg << endl;
      M1_alignment += string(tb_j - hold_tb_j, '-');
      M2_alignment += string(seg.rbegin(), seg.rend());
    } else if(hold_tb_i > 0 && hold_tb_j > 0)  {
      //cerr << "Error in alignning loop regions. Invalid alignment path detected." << endl;
    }
    tb_i = hold_tb_i;
    tb_j = hold_tb_j;
    //cout << M1_alignment << "\t" << M2_alignment << endl;
  }
  
  //cout << M1_alignment << endl;
  //cout << M2_alignment << endl;
  //  reverse the alignment
  M1_alignment = string(M1_alignment.rbegin(), M1_alignment.rend());
  M2_alignment = string(M2_alignment.rbegin(), M2_alignment.rend());
  // chomp the last character because it corresponds to the matching of e1-e2
  M1_alignment.resize(M1_alignment.length() - 1);
  M2_alignment.resize(M2_alignment.length() - 1);
  
  //  evaluate alignment score, remove the sequence alignment scores for beginning and ending nucleotides
  int alignment_score = DP_table[len_seg_1 + 1][len_seg_2 + 1];
  alignment_score -= scoring_function.weight_sequence * scoring_function.match_nucleotide(M1_seq[e1], M2_seq[e2]);
  //cout << "entered here:  " << b1 << "  " << b2 << "  " << e1 << "  " << e2 << endl;
  //cout << "alignment score: " << alignment_score << endl;
  //cout << M1_alignment << endl;
  //cout << M2_alignment << endl;
  //cout << "Score: " << alignment_score << endl;
  //cout << "in loop align: " << alignment_score << endl;
  return alignment_score;
}


int MotifGraphMatching::loop_align(string &M1_alignment, string &M2_alignment, int e1, int e2)  {
  
  //  compute the semi-local alignment of the two sequences
  int i, j, k, l;
  //  clear the alignment holder
  M1_alignment = string("");
  M2_alignment = string("");
  //  if either or both segment is empty
  if(e1 <= 0)  {
    if(e2 <= 0)  {
      //  need more sophisticated STACKING SUBSTITUTION function
      return 0;
    }  else  {
      //  construct alignment
      M1_alignment = string(1, '~');
      M2_alignment = string(1, '~');
      return 0;
    }
  }  else if(e2 <= 0)  {
    //  construct alignment
    M1_alignment = string(1, '~');
    M2_alignment = string(1, '~');
    return 0;
  }
  
  //  allocate space
  vector< vector<int> > DP_table(e1 + 1);
  vector< vector<int> > record_i(e1 + 1);
  vector< vector<int> > record_j(e1 + 1);
  for(i = 0; i < e1 + 1; ++ i)  {
    DP_table[i] = vector<int> (e2 + 1);
    record_i[i] = vector<int> (e2 + 1);
    record_j[i] = vector<int> (e2 + 1);
  }
  
  //  boundary case
  for(i = 0; i <= e1; ++ i)  {
    DP_table[i][0] = scoring_function.weight_sequence * scoring_function.match_nucleotide(M1_seq[i], M2_seq[0]);
    //  label as terminal
    record_i[i][0] = record_j[i][0] = -1;
  }
  for(j = 1; j <= e2; ++ j)  {
    DP_table[0][j] = scoring_function.weight_sequence * scoring_function.match_nucleotide(M1_seq[0], M2_seq[j]);
    //  label as terminal
    record_i[0][j] = record_j[0][j] = -1;
  }
  
  //  general cases
  for(i = 1; i <= e1; ++ i)  {
    for(j = 1; j <= e2; ++ j)  {
      int max_score = -MAX, temp_score;
      int max_i, max_j;
      //  go over all entries in 0...i - 1, j - 1
      for(k = 0; k < i; ++ k)  {
        //  the previous alignment ends at k, j - 1
        temp_score = DP_table[k][j - 1];
        temp_score += scoring_function.weight_adjacent_stacking 
          * scoring_function.match_stacking_basepair(M1_effective_stacking[k][i], M2_effective_stacking[j - 1][j]);
        temp_score += scoring_function.weight_sequence * M1_nuc_del_penalty[k][i];
        temp_score += scoring_function.weight_sequence * scoring_function.match_nucleotide(M1_seq[i], M2_seq[j]);
        if(temp_score > max_score)  {
          max_score = temp_score;
          max_i = k;
          max_j = j - 1;
        }
      }
      //  go over all entries in i - 1, 0...j - 1
      for(k = 0; k < j; ++ k)  {
        //  the previous alignment ends at k, j - 1
        temp_score = DP_table[i - 1][k];
        temp_score += scoring_function.weight_adjacent_stacking 
          * scoring_function.match_stacking_basepair(M1_effective_stacking[i - 1][i], M2_effective_stacking[k][j]);
        temp_score += scoring_function.weight_sequence * M2_nuc_del_penalty[k][j];
        temp_score += scoring_function.weight_sequence * scoring_function.match_nucleotide(M1_seq[i], M2_seq[j]);
        if(temp_score > max_score)  {
          max_score = temp_score;
          max_i = i - 1;
          max_j = k;
        }
      }
      //  go over all 0...breakpoint_i, 0...breakpoint_j
      //  identify the last breakpoints for M1 and M2, respectively
      int M1_last = -1, M2_last = -1;
      for(k = 0; k < (int) M1_break_points.size(); ++ k)  {
        if(M1_break_points[k] > M1_last && M1_break_points[k] <= i)  {
          M1_last = M1_break_points[k];
        }
      }
      for(k = 0; k < (int) M2_break_points.size(); ++ k)  {
        if(M2_break_points[k] > M2_last && M2_break_points[k] <= j)  {
          M2_last = M2_break_points[k];
        }
      }
      for(k = 0; k < M1_last; ++ k)  {
        for(l = 0; l < M2_last; ++ l)  {
          //  the previous alignment ends at k, j - 1
          temp_score = DP_table[k][l];
          temp_score += scoring_function.weight_adjacent_stacking 
            * scoring_function.match_stacking_basepair(M1_effective_stacking[k][i], M2_effective_stacking[l][j]);
          temp_score += scoring_function.weight_sequence * scoring_function.match_nucleotide(M1_seq[i], M2_seq[j]);
          if(temp_score > max_score)  {
            max_score = temp_score;
            max_i = k;
            max_j = l;
          }
        }
      }
      //  assign score
      //  local alignment, taking the matching of i, j as start, need more sophisticated NUC SUBSTITUTION function
      int initial_match_score = (M1_seq[i] == M2_seq[j]) ? 1 : -1;
      if(max_score >  initial_match_score)  {
        DP_table[i][j] = max_score;
        record_i[i][j] = max_i;
        record_j[i][j] = max_j;
      }  else  {
        DP_table[i][j] = initial_match_score;
        //  label as terminal
        record_i[i][j] = record_j[i][j] = -1;
      }
    }
  }
  
  //  finds the maximum alignment score from edge
  int maximum_align_score = DP_table[e1][e2];
  int tb_i = e1, tb_j = e2;  
  while(tb_i > -1 && tb_j > -1)  {
    //  add tb_i and tb_j to the existing alignment
    M1_alignment.append(1, M1_seq[tb_i]);
    M2_alignment.append(1, M2_seq[tb_j]);
    if(record_i[tb_i][tb_j] == tb_i - 1 || record_j[tb_i][tb_j] == tb_j - 1)  {
      //  no break point is involved
      if(record_i[tb_i][tb_j] == tb_i - 1 && record_j[tb_i][tb_j] < tb_j - 1)  {
        int len = tb_j - record_j[tb_i][tb_j] - 1;
        for(k = 0; k < len; ++ k)  {
          M1_alignment.append(1, '-');
          M2_alignment.append(1, M2_seq[tb_j - k - 1]);
        }
      }  else if(record_i[tb_i][tb_j] < tb_i - 1 && record_j[tb_i][tb_j] == tb_j - 1)  {
        int len = tb_i - record_i[tb_i][tb_j] - 1;
        for(k = 0; k < len; ++ k)  {
          M1_alignment.append(1, M1_seq[tb_i - k - 1]);
          M2_alignment.append(1, '-');
        }
      }
    }  else  {
      //  a local alignment across break point is taken
      M1_alignment.append(1, '~');
      M2_alignment.append(1, '~');
    }
    int hold_tb_i = record_i[tb_i][tb_j];
    int hold_tb_j = record_j[tb_i][tb_j];
    tb_i = hold_tb_i;
    tb_j = hold_tb_j;
  }
  
  //  reverse the alignment
  for(i = 0; i < (int) M1_alignment.size() / 2; ++ i)  {
    char temp = M1_alignment[i];
    M1_alignment[i] = M1_alignment[M1_alignment.size() - 1 - i];
    M1_alignment[M1_alignment.size() - 1 - i] = temp;
  }
  
  for(i = 0; i < (int) M2_alignment.size() / 2; ++ i)  {
    char temp = M2_alignment[i];
    M2_alignment[i] = M2_alignment[M2_alignment.size() - 1 - i];
    M2_alignment[M2_alignment.size() - 1 - i] = temp;
  }
  //  delete the alignment of e1 and e2 themselves
  M1_alignment.resize(M1_alignment.size() - 1);
  M2_alignment.resize(M2_alignment.size() - 1);
  /*
  for(i = 0; i < M1_alignment.size(); ++ i)  {
    cout << M1_alignment[i];
  }
  cout << endl;
  for(i = 0; i < M2_alignment.size(); ++ i)  {
    cout << M2_alignment[i];
  }
  cout << endl;
  */
  /*
  for(i = 0; i <= e1; ++ i)  {
    for(j = 0; j <= e2; ++ j)  {
      cout << DP_table[i][j] << "  ";
    }
    cout << endl;
  }
  */
  return maximum_align_score;
}

int MotifGraphMatching::loop_align(int b1, int b2, string &M1_alignment, string &M2_alignment)  {
  //  compute the global alignment of the two sequences
  int i, j, k, l;
  //  clear the alignment holder
  M1_alignment = string("");
  M2_alignment = string("");
  //  compute effective length
  int len_seg_1 = M1_seq.length() - b1;
  int len_seg_2 = M2_seq.length() - b2;
  
  //  if either or both segment is empty
  if(len_seg_1 <= 1)  {
    if(len_seg_2 <= 1)  {
      return 0;
    }  else  {
      //  construct alignment
      M1_alignment = string(1, '~');
      M2_alignment = string(1, '~');
      return 0;
    }
  }  else if(len_seg_2 <= 1)  {
    //  construct alignment
    M1_alignment = string(1, '~');
    M2_alignment = string(1, '~');
    return 0;
  }
  
  //  allocate space
  vector< vector<int> > DP_table(len_seg_1);
  vector< vector<int> > record_i(len_seg_1);
  vector< vector<int> > record_j(len_seg_1);
  for(i = 0; i < len_seg_1; ++ i)  {
    DP_table[i] = vector<int> (len_seg_2);
    record_i[i] = vector<int> (len_seg_2);
    record_j[i] = vector<int> (len_seg_2);
  }
  
  //  boundary case
  for(i = 0; i < len_seg_1; ++ i)  {
    DP_table[i][0] = -MAX;
    record_i[i][0] = record_j[i][0] = -1;
  }
  for(j = 1; j < len_seg_2; ++ j)  {
    DP_table[0][j] = -MAX;
    record_i[0][j] = record_j[0][j] = -1;
  }
  //  need more sophisticated NUC SUBSTITUTION scoring function
  DP_table[0][0] = scoring_function.weight_sequence * scoring_function.match_nucleotide(M1_seq[b1], M2_seq[b2]);
  
  //  general cases
  for(i = 1; i < len_seg_1; ++ i)  {
    for(j = 1; j < len_seg_2; ++ j)  {
      int max_score = -MAX, temp_score;
      int max_i = -1, max_j = -1;
      //  go over all entries in 0...i - 1, j - 1
      for(k = 0; k < i; ++ k)  {
        //  the previous alignment ends at k, j - 1
        temp_score = DP_table[k][j - 1];
        temp_score += scoring_function.weight_adjacent_stacking 
          * scoring_function.match_stacking_basepair(M1_effective_stacking[b1 + k][b1 + i], M2_effective_stacking[b2 + j - 1][b2 + j]);
        temp_score += scoring_function.weight_sequence * M1_nuc_del_penalty[b1 + k][b1 + i];
        temp_score += scoring_function.weight_sequence * scoring_function.match_nucleotide(M1_seq[b1 + i], M2_seq[b2 + j]);
        if(temp_score > max_score)  {
          max_score = temp_score;
          max_i = k;
          max_j = j - 1;
        }
      }
      //  go over all entries in i - 1, 0...j - 1
      for(k = 0; k < j; ++ k)  {
        //  the previous alignment ends at k, j - 1
        temp_score = DP_table[i - 1][k];
        temp_score += scoring_function.weight_adjacent_stacking 
          * scoring_function.match_stacking_basepair(M1_effective_stacking[b1 + i - 1][b1 + i], M2_effective_stacking[b2 + k][b2 + j]);
        temp_score += scoring_function.weight_sequence * M2_nuc_del_penalty[b2 + k][b2 + j];
        temp_score += scoring_function.weight_sequence * scoring_function.match_nucleotide(M1_seq[b1 + i], M2_seq[b2 + j]);
        if(temp_score > max_score)  {
          max_score = temp_score;
          max_i = i - 1;
          max_j = k;
        }
      }
      //  go over all 0...breakpoint_i, 0...breakpoint_j
      //  identify the last breakpoints for M1 and M2, respectively
      int M1_last = -1, M2_last = -1;
      for(k = 0; k < (int) M1_break_points.size(); ++ k)  {
        if(M1_break_points[k] - b1 > M1_last && M1_break_points[k] - b1 <= i)  {
          M1_last = M1_break_points[k] - b1;
        }
      }
      for(k = 0; k < (int) M2_break_points.size(); ++ k)  {
        if(M2_break_points[k] - b2 > M2_last && M2_break_points[k] - b2 <= j)  {
          M2_last = M2_break_points[k] - b2;
        }
      }
      for(k = 0; k < M1_last; ++ k)  {
        for(l = 0; l < M2_last; ++ l)  {
          //  the previous alignment ends at k, j - 1
          temp_score = DP_table[k][l];
          temp_score += scoring_function.weight_adjacent_stacking 
            * scoring_function.match_stacking_basepair(M1_effective_stacking[k + b1][i + b1], M2_effective_stacking[l + b2][j + b2]);
          temp_score += scoring_function.weight_sequence * scoring_function.match_nucleotide(M1_seq[i + b1], M2_seq[j + b2]);
          if(temp_score > max_score)  {
            max_score = temp_score;
            max_i = k;
            max_j = l;
          }
        }
      }
      //  assign score
      DP_table[i][j] = max_score;
      record_i[i][j] = max_i;
      record_j[i][j] = max_j;
    }
  }
  
  
  //  identify the maximum alignment score
  int maximal_alignment_score = -MAX;
  int tb_i = -1, tb_j = -1;
  for(i = 0; i < len_seg_1; ++ i)  {
    for(j = 0; j < len_seg_2; ++ j)  {
      if(DP_table[i][j] > maximal_alignment_score)  {
        maximal_alignment_score = DP_table[i][j];
        tb_i = i;
        tb_j = j;
      }
    }
  }
  
  
  //  traceback and recover alignment
  while(tb_i > 0 && tb_j > 0)  {
    //  add tb_i and tb_j to the existing alignment
    M1_alignment.append(1, M1_seq[b1 + tb_i]);
    M2_alignment.append(1, M2_seq[b2 + tb_j]);
    if(record_i[tb_i][tb_j] == tb_i - 1 || record_j[tb_i][tb_j] == tb_j - 1)  {
      //  no break point is involved
      if(record_i[tb_i][tb_j] == tb_i - 1 && record_j[tb_i][tb_j] < tb_j - 1)  {
        int len = tb_j - record_j[tb_i][tb_j] - 1;
        for(k = 0; k < len; ++ k)  {
          M1_alignment.append(1, '-');
          M2_alignment.append(1, M2_seq[b2 + tb_j - k - 1]);
        }
      }  else if(record_i[tb_i][tb_j] < tb_i - 1 && record_j[tb_i][tb_j] == tb_j - 1)  {
        int len = tb_i - record_i[tb_i][tb_j] - 1;
        for(k = 0; k < len; ++ k)  {
          M1_alignment.append(1, M1_seq[b1 + tb_i - k - 1]);
          M2_alignment.append(1, '-');
        }
      }
    }  else  {
      //  a local alignment across break point is taken
      M1_alignment.append(1, '~');
      M2_alignment.append(1, '~');
    }
    int hold_tb_i = record_i[tb_i][tb_j];
    int hold_tb_j = record_j[tb_i][tb_j];
    tb_i = hold_tb_i;
    tb_j = hold_tb_j;
  }
  
  
  //  reverse the alignment
  for(i = 0; i < (int) M1_alignment.size() / 2; ++ i)  {
    char temp = M1_alignment[i];
    M1_alignment[i] = M1_alignment[M1_alignment.size() - 1 - i];
    M1_alignment[M1_alignment.size() - 1 - i] = temp;
  }
  
  for(i = 0; i < (int) M2_alignment.size() / 2; ++ i)  {
    char temp = M2_alignment[i];
    M2_alignment[i] = M2_alignment[M2_alignment.size() - 1 - i];
    M2_alignment[M2_alignment.size() - 1 - i] = temp;
  }
    
  //  evaluate alignment score, remove the sequence alignment scores for beginning and ending nucleotides
  maximal_alignment_score -= scoring_function.weight_sequence * scoring_function.match_nucleotide(M1_seq[b1], M2_seq[b2]);
  return maximal_alignment_score;
}

bool MotifGraphMatching::significant_match(vector<int> interaction_1, vector<int> interaction_2)  {
  if(interaction_1[2] > 0 && interaction_2[2] > 0 && interaction_1[3] == interaction_2[3] &&
    interaction_1[4] == interaction_2[4] && interaction_1[5] == interaction_2[5] && interaction_1[4] != (int) 'W' &&
    interaction_1[5] != 1)  {
    //  similar type of non-canonical base pairs, same edges orientation
    return true;
  }  else  {
    return false;
  }
}

bool MotifGraphMatching::is_isosteric(vector<int> interaction_1, vector<int> interaction_2)  {
  if(interaction_1[2] > 0 && interaction_2[2] > 0 && interaction_1[2] != MAX && interaction_2[2] != MAX &&
    scoring_function.basepair_index[interaction_1[2]] == scoring_function.basepair_index[interaction_2[2]])  {
    if(scoring_function.basepair_index[interaction_1[2]] != 0 && scoring_function.basepair_index[interaction_1[2]] != 1)  {
      // if both are canonical base pair, then they does not counted as isosteric matching
      return true;
    } else  {
      return false;
    }
  }  else  {
    return false;
  }
}

int MotifGraphMatching::maximal_seq_score(int begin, int end, string M_seq, vector <int> M_chain_stacking)  {
  int i;
  int max_score = 0;
  
  //cout << "maximal_seq_score called:  " << begin << "  " << end << endl;
  //  evaluate nucleotide substitution
  for(i = begin; i <= end; ++ i)  {
    max_score += scoring_function.weight_sequence * scoring_function.match_nucleotide(M_seq[i], M_seq[i]);
  }
  //cout << "max_score 1:  " << max_score << endl;
  //  evaluate adjacent stacking substitution
  for(i = 0; i < (int) M_chain_stacking.size(); ++ i)  {
    if(M_chain_stacking[i] >= begin && M_chain_stacking[i] < end)  {
      max_score += scoring_function.weight_adjacent_stacking * 
        scoring_function.match_stacking_basepair(-1, -1);
    }
  }
  //cout << "max_score 2:  " << max_score << endl;
  return max_score;
}

int MotifGraphMatching::asymmetric_penalty(std::vector<std::vector<int> >& matched, int p1, int p2)  {

  int i, j, k;  
  if(matched.size() == 0 || M1_interaction[p1][2] <= 2)  {
    return 0;
  }
  //  finds the non-canonical base pairs that are directly before the current matched pairs p1 and p2
  //  sort the base pairs based on their locations
  for(i = 0; i < (int) matched.size() - 1; ++ i)  {
    for(j = i + 1; j < (int) matched.size(); ++ j)  {
      if((M1_interaction[matched[i][0]][0] > M1_interaction[matched[j][0]][0]) || 
        (M1_interaction[matched[i][0]][0] == M1_interaction[matched[j][0]][0] && 
         M1_interaction[matched[i][0]][1] > M1_interaction[matched[j][0]][1]))  {
        vector<int> temp(2);
        temp = matched[i];
        matched[i] = matched[j];
        matched[j] = temp;
      }
    }
  }
  int record_i = 0;
  //cout << (int) matched.size() << endl;
  for(i = 0; i < (int) matched.size(); ++ i)  {
    //cout << i << endl;
    if((M1_interaction[matched[i][0]][0] > M1_interaction[p1][0]) && 
        (M1_interaction[matched[i][0]][1] < M1_interaction[p1][1]) &&
        (M1_interaction[matched[i][0]][2] > 2))  {
      record_i = i;
      break;
    }
  }
  vector <int> adjacent_pairs;
  adjacent_pairs.push_back(record_i);
  
  int max_index;
  if(i >= 0)  {
    max_index = M1_interaction[matched[record_i][0]][1];
  }  else  {
    return 0;
  }
  for(i = record_i + 1; i < (int) matched.size(); ++ i)  {
    if(M1_interaction[matched[i][0]][1] > max_index && M1_interaction[matched[i][0]][1] < M1_interaction[p1][1])  {
      max_index = M1_interaction[matched[i][0]][1];
      adjacent_pairs.push_back(i);
    }
  }  
  int penalty_score = 0;
  //  compute asymmetricity between the loops
  for(i = 0; i < (int) adjacent_pairs.size(); ++ i)  {
    k = adjacent_pairs[i];
    int sp1l = (M1_interaction[matched[k][0]][0] - M1_interaction[p1][0]);
    int sp1r = (M1_interaction[p1][1] - M1_interaction[matched[k][0]][1]);
    int a1 =  sp1l - sp1r;
    int sp2l = (M2_interaction[matched[k][1]][0] - M2_interaction[p2][0]);
    int sp2r = (M2_interaction[p2][1] - M2_interaction[matched[k][1]][1]);
    int a2 =  sp2l - sp2r; 
    int asym_size = abs(a1 - a2);
    penalty_score += scoring_function.weight_isosteric * asym_size * scoring_function.asym_nuc;
    if(a1 * a2 < 0){
      penalty_score += scoring_function.weight_isosteric * scoring_function.asym_loop;
    }
      
    if(sp1l == 1 && sp1r == 1)  {
      //  in case of a stacking case
      //cout << "entered" << endl;
      if(sp2l == 1 && sp2r == 1 && is_isosteric(M1_interaction[matched[k][0]], M2_interaction[matched[k][1]]) && 
        is_isosteric(M1_interaction[p1], M2_interaction[p2]))  {
        penalty_score += scoring_function.weight_isosteric * scoring_function.cons_stacking;
      }  else if((sp2l > 1 || sp2r > 1) && sp2l != sp2r)  {
        //  asymmetric internal loop or bulge
        if(sp2l == 1 || sp2r == 1)  {
          //  a bulge
          penalty_score += scoring_function.weight_isosteric * scoring_function.stack_to_bulge;
        }  else  {
          //  an aymmetric internal loop
          penalty_score += scoring_function.weight_isosteric * scoring_function.stack_to_internal_asym;
        }
      }  else if((sp2l > 1 || sp2r > 1) && sp2l == sp2r)  {
        //  symmetric internal loop
        penalty_score += scoring_function.weight_isosteric * scoring_function.stack_to_internal_sym;
      }
    }
  }
  //cout << "asymmetric penalty score:  " << penalty_score << endl;
  return penalty_score;
}

void MotifGraphMatching::find_aligned_boundary(const vector<int>& clique, vector<int>& boundary)  {
  vector<int> M1_nuc, M2_nuc;
  get_aligned_nucs(clique, M1_nuc, M2_nuc);
  //cout << "Record_alignment:  stuck here 1" << endl;
  list<AlignedRegionType> M1_regions, M2_regions;
  get_aligned_regions(M1_nuc, M1_break_points, M1_regions);
  get_aligned_regions(M2_nuc, M2_break_points, M2_regions);
  for(auto it = M1_regions.begin(); it != M1_regions.end(); ++ it) {
    boundary.push_back(it->begin);
    boundary.push_back(it->end);
  }
  return;
}

void MotifGraphMatching::CountDeletedInteractions(
   const std::vector<int>& aligned_pairs, const std::vector<std::vector<int> >& interactions, 
   const std::list<AlignedRegionType>& aligned_region, int& num_pair_deleted, int& num_stack_deleted
) {
  //cout << "******************" << endl;
  // initialization
  num_pair_deleted = num_stack_deleted = 0;
  // for each base pair that is not aligned
  for(int i = 0; i < (int) aligned_pairs.size(); ++ i)  {
    //cout << "Interaction: " << i << endl;
    if(!aligned_pairs[i])  {
      bool local_tag = true;
      for(auto it = aligned_region.begin(); it != aligned_region.end(); ++ it)  {
        //cout << "Region:  " << it->begin << "  " << it->end << " --  " << endl;
        if((interactions[i][0] >= it->begin && interactions[i][0] <= it->end) ||
          (interactions[i][1] >= it->begin && interactions[i][1] <= it->end))  {
          //cout << it->begin << "  " << it->end << " --  " << M1_interaction[i][0] << "  " << M1_interaction[i][1] << endl;
          local_tag = false;
          break;
        }
      }
      if(!local_tag)  {
        if(interactions[i][2] > 0 && interactions[i][2] != MAX)  {
          ++ num_pair_deleted;
        }  else  if(interactions[i][2] < 0) {
          ++ num_stack_deleted;
        }
      }
    }
  }  
  return;
}

int MotifGraphMatching::compute_alignment_score(const vector<int>& estimated_clique)  {
  //  compute the alignment score given the base-interaction matching retrieved from estimated_clique
  int i;
  vector<int> M1_nuc, M2_nuc;
  vector<int> M1_matched, M2_matched;
  M1_matched = vector<int> (M1_seq.length(), 0);
  M2_matched = vector<int> (M2_seq.length(), 0);
  int evaluate_score = 0;
  vector< vector<int> > matched_pairs;
  vector<int> M1_pairs_record(M1_interaction.size(), 0);
  vector<int> M2_pairs_record(M2_interaction.size(), 0);
  //  construct the nucleotides that particitate in selected interactions
  for(i = 0; i < (int) estimated_clique.size(); ++ i)  {
    //cout << "index: " << i << endl;
    //  the two base interactions idx_i and idx_j match  
    int idx_i = estimated_clique[i] / M2_num_interaction;
    int idx_j = estimated_clique[i] % M2_num_interaction;
    M1_pairs_record[idx_i] = M2_pairs_record[idx_j] = 1;
    //  check for insertions that results in a bulge
    
    evaluate_score += asymmetric_penalty(matched_pairs, idx_i, idx_j);
    //cout << "evaluate_score asymm:  " << evaluate_score << endl;
    //  push the newly matched pairs  
    vector<int> current_matched = {idx_i, idx_j};
    matched_pairs.push_back(current_matched);
    //cout << "Processing interactions:  " << idx_i << "  " << idx_j << endl;
    if((M1_interaction[idx_i][2] > 0 && M2_interaction[idx_j][2] > 0) ||
      (M1_interaction[idx_i][2] > 0 && M2_interaction[idx_j][2] == 0) ||
      (M1_interaction[idx_i][2] == 0 && M2_interaction[idx_j][2] > 0))  {
      //  match of two base pairs, hbond may present
      //cout << "matching basepair" << endl;
      evaluate_score += scoring_function.weight_isosteric 
        * scoring_function.match_isosteric_basepair(M1_interaction[idx_i][2], M2_interaction[idx_j][2]);
      //cout << "basepair type:  " << M1_interaction[idx_i][2] << "  " << M2_interaction[idx_j][2] << endl;
      //cout << "pair score:  " << scoring_function.match_isosteric_basepair(M1_interaction[idx_i][2], M2_interaction[idx_j][2]) << endl;
    }  else if((M1_interaction[idx_i][2] < 0 && M2_interaction[idx_j][2] < 0) ||
      (M1_interaction[idx_i][2] < 0 && M2_interaction[idx_j][2] == 0) ||
      (M1_interaction[idx_i][2] == 0 && M2_interaction[idx_j][2] < 0))  {
      //  match of two non-adjancent stackings. hbond may present
      evaluate_score += scoring_function.weight_nonadjacent_stacking
        * scoring_function.match_stacking_basepair(M1_interaction[idx_i][2], M2_interaction[idx_j][2]);
      //cout << "basepair type:  " << M1_interaction[idx_i][2] << "  " << M2_interaction[idx_j][2] << endl;
      //cout << "stacking score:  " << scoring_function.match_stacking_basepair(M1_interaction[idx_i][2], M2_interaction[idx_j][2]) << endl;
    }  else if(M1_interaction[idx_i][2] == 0 && M2_interaction[idx_j][2] == 0)  {
      //  match of two hydrogen bonds
      evaluate_score += scoring_function.weight_isosteric * 
        (scoring_function.hbond_match_base + scoring_function.hbond_match_bonus);
      //cout << "3:  " << scoring_function.hbond_match_bonus << endl;
    }
    //cout << "evaluate score:  " << evaluate_score << endl;
    evaluate_score += scoring_function.weight_sequence 
      * (scoring_function.match_nucleotide(M1_seq[M1_interaction[idx_i][0]], M2_seq[M2_interaction[idx_j][0]]) 
      + scoring_function.match_nucleotide(M1_seq[M1_interaction[idx_i][1]], M2_seq[M2_interaction[idx_j][1]]));
    //cout << "nuc match score: " << scoring_function.match_nucleotide(M1_seq[M1_interaction[idx_i][0]], M2_seq[M2_interaction[idx_j][0]]) + scoring_function.match_nucleotide(M1_seq[M1_interaction[idx_i][1]], M2_seq[M2_interaction[idx_j][1]]) << endl;
    //cout << "evaluate score:  " << evaluate_score << endl;
    //  copy the matched indices
    M1_nuc.push_back(M1_interaction[idx_i][0]);
    M1_nuc.push_back(M1_interaction[idx_i][1]);
    M2_nuc.push_back(M2_interaction[idx_j][0]);
    M2_nuc.push_back(M2_interaction[idx_j][1]);
    //  check for bonus score of triple interaction
    if(M1_matched[M1_interaction[idx_i][0]] > 0 && M2_matched[M2_interaction[idx_j][0]] > 0)  {
      //cout << "Triple bonus added" << endl;
      evaluate_score += scoring_function.weight_isosteric * scoring_function.triple_interaction_bonus;
    }
    if(M1_matched[M1_interaction[idx_i][1]] > 0 && M2_matched[M2_interaction[idx_j][1]] > 0)  {
      //cout << "Triple bonus added" << endl;
      evaluate_score += scoring_function.weight_isosteric * scoring_function.triple_interaction_bonus;
    }
    if(significant_match(M1_interaction[idx_i], M2_interaction[idx_j]))  {
      M1_matched[M1_interaction[idx_i][0]] = M1_matched[M1_interaction[idx_i][1]] = 1;
      M2_matched[M2_interaction[idx_j][0]] = M2_matched[M2_interaction[idx_j][1]] = 1;      
    }
  }
  //cout << "End of pair evaluation score:  " << evaluate_score << endl;
  //  assert that the sizes of the matched nucleotides are the same
  if(M1_nuc.size() != M2_nuc.size())  {
    //cout << "The sizes of matched nucleotides do not match: MotifGraphMatching::compute_alignment_score" << endl;
    return 0;
  }  else if(M1_nuc.size() == 0)  {
    return loop_align();
  }
  //cout << "Esitmate score:  " << evaluate_score << endl;
  //  sort the matched nucleotides
  sort(M1_nuc.begin(), M1_nuc.end());
  sort(M2_nuc.begin(), M2_nuc.end());
  //  for each segment evaluate the corresponding alignment score
  string M1_hold_alignment;
  string M2_hold_alignment;
  //evaluate_score += loop_align(M1_hold_alignment, M2_hold_alignment, M1_nuc[0], M2_nuc[0]);
  /*
  for(auto it = M1_nuc.begin(); it != M1_nuc.end(); ++ it) {
    cout << *it << "  ";
  }
  cout << endl;
  for(auto it = M2_nuc.begin(); it != M2_nuc.end(); ++ it) {
    cout << *it << "  ";
  }
  cout << endl;
  */
  //cout << M1_nuc[0] << "  " << M2_nuc[0] << endl;
  //cout << loop_align(M1_hold_alignment, M2_hold_alignment, M1_nuc[0], M2_nuc[0]) << endl;
  for(i = 0; i < (int) M1_nuc.size() - 1; ++ i)  {
    //  evaluate the alignment scrore of the segment
    if(M1_nuc[i] != M1_nuc[i + 1] || M2_nuc[i] != M2_nuc[i + 1])  {
      //cout << "call loop align: " << M1_nuc[i] << " " <<  M2_nuc[i] << "  " << M1_nuc[i + 1] << " " <<  M2_nuc[i + 1] << endl;
      evaluate_score += loop_align(M1_nuc[i], M2_nuc[i], M1_hold_alignment, M2_hold_alignment, M1_nuc[i + 1], M2_nuc[i + 1]);
      //cout << "Evaluate score:  " << evaluate_score << endl;
    }
  }
  //cout << "score after loop align:  " << evaluate_score << endl;
  //evaluate_score += loop_align(M1_nuc[M1_nuc.size() - 1], M2_nuc[M2_nuc.size() - 1], M1_hold_alignment, M2_hold_alignment);
  //  apply penalty for unmatched base pairs in the motif query
  vector<int> M1_nuc_check, M2_nuc_check;
  get_aligned_nucs(estimated_clique, M1_nuc_check, M2_nuc_check);
  //cout << "Record_alignment:  stuck here 1" << endl;
  list<AlignedRegionType> M1_regions, M2_regions;
  get_aligned_regions(M1_nuc_check, M1_break_points, M1_regions);
  get_aligned_regions(M2_nuc_check, M2_break_points, M2_regions);
  
  int unaligned_pair_query = 0, unaligned_stacking_query = 0;
  int unaligned_pair_target = 0, unaligned_stacking_target = 0;
  
  CountDeletedInteractions(
      M1_pairs_record, M1_interaction, M1_regions, unaligned_pair_query, unaligned_stacking_query
  );
  CountDeletedInteractions(
      M2_pairs_record, M2_interaction, M2_regions, unaligned_pair_target, unaligned_stacking_target
  );
  
  //cout << "Num deleted query interactions:  " << unaligned_pair_query << "  " << unaligned_stacking_query << endl;
  //cout << "Num deleted target interactions:  " << unaligned_pair_target << "  " << unaligned_stacking_target << endl;
  evaluate_score += scoring_function.weight_isosteric * (-1200 * unaligned_pair_query) 
      + scoring_function.weight_nonadjacent_stacking * (-800 * unaligned_stacking_query);
  evaluate_score += scoring_function.weight_isosteric * (-120 * unaligned_pair_target) 
      + scoring_function.weight_nonadjacent_stacking * (-80 * unaligned_stacking_target);
  //cout << "score after deletion penalty:  " << evaluate_score << endl;
  return evaluate_score;
}

int MotifGraphMatching::compute_alignment_score_max(const vector<int>& estimated_clique)  {
  //  compute the alignment score given the base-interaction matching retrieved from estimated_clique
  vector<int> M1_nuc, M2_nuc;
  vector<int> M1_matched, M2_matched;
  M1_matched = vector<int> (M1_seq.length(), 0);
  M2_matched = vector<int> (M2_seq.length(), 0);
  int i, j;
  int evaluate_score = 0;
  vector< vector<int> > matched_pairs;
  vector<int> M1_pairs_record(M1_interaction.size(), 0);
  vector<int> M2_pairs_record(M2_interaction.size(), 0);
  //  construct the nucleotides that particitate in selected interactions
  for(i = 0; i < (int) estimated_clique.size(); ++ i)  {
    //  the two base interactions idx_i and idx_j match
    int idx_i = estimated_clique[i] / M2_num_interaction;
    int idx_j = estimated_clique[i] % M2_num_interaction;
    M1_pairs_record[idx_i] = M2_pairs_record[idx_j] = 1;
    //  check for insertions that results in a bulge
    evaluate_score += asymmetric_penalty(matched_pairs, idx_i, idx_j);
    //  push the newly matched pairs
    matched_pairs.resize(matched_pairs.size() + 1);
    matched_pairs[matched_pairs.size() - 1] = vector<int> (2);
    matched_pairs[matched_pairs.size() - 1][0] = idx_i;
    matched_pairs[matched_pairs.size() - 1][1] = idx_j;
    //cout << "Processing interactions:  " << idx_i << "  " << idx_j << endl;
    if((M1_interaction[idx_i][2] > 0 && M2_interaction[idx_j][2] > 0) ||
      (M1_interaction[idx_i][2] > 0 && M2_interaction[idx_j][2] == 0) ||
      (M1_interaction[idx_i][2] == 0 && M2_interaction[idx_j][2] > 0))  {
      //  match of two base pairs, hbond may present
      //cout << "matching basepair" << endl;
      evaluate_score += scoring_function.weight_isosteric 
        * scoring_function.match_isosteric_basepair(M1_interaction[idx_i][2], M2_interaction[idx_j][2]);
      //cout << "basepair type:  " << M1_interaction[idx_i][2] << "  " << M2_interaction[idx_j][2] << endl;
      //cout << "score:  " << scoring_function.match_isosteric_basepair(M1_interaction[idx_i][2], M2_interaction[idx_j][2]) << endl;
    }  else if((M1_interaction[idx_i][2] < 0 && M2_interaction[idx_j][2] < 0) ||
      (M1_interaction[idx_i][2] < 0 && M2_interaction[idx_j][2] == 0) ||
      (M1_interaction[idx_i][2] == 0 && M2_interaction[idx_j][2] < 0))  {
      //  match of two non-adjancent stackings. hbond may present
      evaluate_score += scoring_function.weight_nonadjacent_stacking
        * scoring_function.match_stacking_basepair(M1_interaction[idx_i][2], M2_interaction[idx_j][2]);
      //cout << "basepair type:  " << M1_interaction[idx_i][2] << "  " << M2_interaction[idx_j][2] << endl;
      //cout << "score:  " << scoring_function.match_stacking_basepair(M1_interaction[idx_i][2], M2_interaction[idx_j][2]) << endl;
    }  else if(M1_interaction[idx_i][2] == 0 && M2_interaction[idx_j][2] == 0)  {
      //  match of two hydrogen bonds
      evaluate_score += scoring_function.weight_isosteric * 
        (scoring_function.hbond_match_base + scoring_function.hbond_match_bonus);
      //cout << "3:  " << scoring_function.hbond_match_bonus << endl;
    }
    evaluate_score += scoring_function.weight_sequence 
      * (scoring_function.match_nucleotide(M1_seq[M1_interaction[idx_i][0]], M2_seq[M2_interaction[idx_j][0]]) 
      + scoring_function.match_nucleotide(M1_seq[M1_interaction[idx_i][1]], M2_seq[M2_interaction[idx_j][1]]));
    //  copy the matched indices
    M1_nuc.push_back(M1_interaction[idx_i][0]);
    M1_nuc.push_back(M1_interaction[idx_i][1]);
    M2_nuc.push_back(M2_interaction[idx_j][0]);
    M2_nuc.push_back(M2_interaction[idx_j][1]);
    //  check for bonus score of triple interaction
    if(M1_matched[M1_interaction[idx_i][0]] > 0 && M2_matched[M2_interaction[idx_j][0]] > 0)  {
      evaluate_score += scoring_function.weight_isosteric * scoring_function.triple_interaction_bonus;
    }
    if(M1_matched[M1_interaction[idx_i][1]] > 0 && M2_matched[M2_interaction[idx_j][1]] > 0)  {
      evaluate_score += scoring_function.weight_isosteric * scoring_function.triple_interaction_bonus;
    }
    if(significant_match(M1_interaction[idx_i], M2_interaction[idx_j]))  {
      M1_matched[M1_interaction[idx_i][0]] = M1_matched[M1_interaction[idx_i][1]] = 1;
      M2_matched[M2_interaction[idx_j][0]] = M2_matched[M2_interaction[idx_j][1]] = 1;      
    }
  }
  //  assert that the sizes of the matched nucleotides are the same
  if(M1_nuc.size() != M2_nuc.size())  {
    cout << "The sizes of matched nucleotides do not match: MotifGraphMatching::compute_alignment_score" << endl;
    return 0;
  }  else if(M1_nuc.size() == 0)  {
    return maximal_seq_score(0, M1_seq.length() - 1, M1_seq, M1_chain_stacking) < 
      maximal_seq_score(0, M2_seq.length() - 1, M2_seq, M2_chain_stacking) ?
      maximal_seq_score(0, M1_seq.length() - 1, M1_seq, M1_chain_stacking) :
      maximal_seq_score(0, M2_seq.length() - 1, M2_seq, M2_chain_stacking);
  }
  //  sort the matched nucleotides
  for(i = 0; i < (int) M1_nuc.size() - 1; ++ i)  {
    for(j = i + 1; j < (int) M1_nuc.size(); ++ j)  {
      //cout << "@  " << i << "  " << j << endl;
      if(M1_nuc[i] > M1_nuc[j])  {
        int temp = M1_nuc[i];
        M1_nuc[i] = M1_nuc[j];
        M1_nuc[j] = temp;
      }
    }
  }
  for(i = 0; i < (int) M2_nuc.size() - 1; ++ i)  {
    for(j = i + 1; j < (int) M2_nuc.size(); ++ j)  {
      //cout << "@  " << i << "  " << j << endl;
      if(M2_nuc[i] > M2_nuc[j])  {
        int temp = M2_nuc[i];
        M2_nuc[i] = M2_nuc[j];
        M2_nuc[j] = temp;
      }
    }
  }
  //  for each segment evaluate the corresponding alignment score
  string M1_hold_alignment;
  string M2_hold_alignment;
  evaluate_score += maximal_seq_score(0, M1_nuc[0], M1_seq, M1_chain_stacking) < 
    maximal_seq_score(0, M2_nuc[0], M2_seq, M2_chain_stacking) ?
    maximal_seq_score(0, M1_nuc[0], M1_seq, M1_chain_stacking) : 
    maximal_seq_score(0, M2_nuc[0], M2_seq, M2_chain_stacking);
  for(i = 0; i < (int) M1_nuc.size() - 1; ++ i)  {
    //  evaluate the alignment scrore of the segment
    if(M1_nuc[i] != M1_nuc[i + 1] || M2_nuc[i] != M2_nuc[i + 1])  {
      evaluate_score += maximal_seq_score(M1_nuc[i], M1_nuc[i + 1], M1_seq, M1_chain_stacking) < 
        maximal_seq_score(M2_nuc[i], M2_nuc[i + 1], M2_seq, M2_chain_stacking) ?
        maximal_seq_score(M1_nuc[i], M1_nuc[i + 1], M1_seq, M1_chain_stacking) : 
        maximal_seq_score(M2_nuc[i], M2_nuc[i + 1], M2_seq, M2_chain_stacking);
    }
  }
  evaluate_score += maximal_seq_score(M1_nuc[M1_nuc.size() - 1], M1_seq.length() - 1, M1_seq, M1_chain_stacking) < 
    maximal_seq_score(M2_nuc[M2_nuc.size() - 1], M2_seq.length() - 1, M2_seq, M2_chain_stacking) ?
    maximal_seq_score(M1_nuc[M1_nuc.size() - 1], M1_seq.length() - 1, M1_seq, M1_chain_stacking) : 
    maximal_seq_score(M2_nuc[M2_nuc.size() - 1], M2_seq.length() - 1, M2_seq, M2_chain_stacking);
  
  
  
  vector<int> M1_nuc_check, M2_nuc_check;
  get_aligned_nucs(estimated_clique, M1_nuc_check, M2_nuc_check);
  //cout << "Record_alignment:  stuck here 1" << endl;
  list<AlignedRegionType> M1_regions, M2_regions;
  get_aligned_regions(M1_nuc_check, M1_break_points, M1_regions);
  get_aligned_regions(M2_nuc_check, M2_break_points, M2_regions);
  
  int unaligned_pair_query = 0, unaligned_stacking_query = 0;
  int unaligned_pair_target = 0, unaligned_stacking_target = 0;
  CountDeletedInteractions(
      M1_pairs_record, M1_interaction, M1_regions, unaligned_pair_query, unaligned_stacking_query
  );
  CountDeletedInteractions(
      M2_pairs_record, M2_interaction, M2_regions, unaligned_pair_target, unaligned_stacking_target
  );
  evaluate_score += scoring_function.weight_isosteric * (-1200 * unaligned_pair_query) 
      + scoring_function.weight_nonadjacent_stacking * (-800 * unaligned_stacking_query);
  evaluate_score += scoring_function.weight_isosteric * (-120 * unaligned_pair_target) 
      + scoring_function.weight_nonadjacent_stacking * (-80 * unaligned_stacking_target);
  return evaluate_score;
}

int MotifGraphMatching::compute_alignment_score_upper_bound(const vector<int>& current_clique, const vector<int>& candidates)  {
  int i, j;
  //cout << "printing from compute_alignment_score_upper_bound" << endl;
  int evaluate_score = 0;
  //  call compute_alignment_score to estimate the score for existing matching
  evaluate_score = compute_alignment_score_max(current_clique);
  //cout << "evaluate_score in compute_alignment_score_max:  " << evaluate_score << endl;
  //  assume all candidates can match, then compute the upper bound
  vector<int> M1_matched(M1_seq.length(), 0);
  vector<int> M2_matched(M2_seq.length(), 0);
  vector<int> M1_pairs_potential(M1_interaction.size(), 0);
  vector<int> M2_pairs_potential(M2_interaction.size(), 0);
  for(i = 0; i < (int) current_clique.size(); ++ i)  {
    //  the two base interactions idx_i and idx_j match
    int idx_i = current_clique[i] / M2_num_interaction;
    int idx_j = current_clique[i] % M2_num_interaction;    
    if(significant_match(M1_interaction[idx_i], M2_interaction[idx_j]))  {
      M1_matched[M1_interaction[idx_i][0]] = M1_matched[M1_interaction[idx_i][1]] = 1;
      M2_matched[M2_interaction[idx_j][0]] = M2_matched[M2_interaction[idx_j][1]] = 1;
    }
  }  
  vector<int> M1_repeat(M1_interaction.size(), 0);
  vector<int> M2_repeat(M2_interaction.size(), 0);
  vector<int> M1_max_score(M1_interaction.size(), 0);
  vector<int> M2_max_score(M2_interaction.size(), 0);
  int num_edges_candidates = 0;
  unordered_map<int, int> degrees_candidates;
  for(i = 0; i < (int) candidates.size(); ++ i) {
    degrees_candidates[candidates[i]] = 0;
  }
  for(i = 0; i < (int) candidates.size() - 1; ++ i)  {
    for(j = i + 1; j < (int) candidates.size(); ++ j)  {
      if(interaction_compatibility[candidates[i]][candidates[j]] == 1)  {
        ++ num_edges_candidates;
        ++ degrees_candidates[candidates[i]];
        ++ degrees_candidates[candidates[j]];
      }
    }
  }
  int min_degree = M1_num_interaction * M2_num_interaction;
  for(auto it = degrees_candidates.begin(); it != degrees_candidates.end(); ++ it)  {
    if(it->second < min_degree)  {
      min_degree = it->second;
    }
  }
  int num_candidates = 0;
  for(i = 0; i < (int) candidates.size(); ++ i)  {  
    ++ num_candidates;
    //  the two base interactions idx_i and idx_j match
    int idx_i = candidates[i] / M2_num_interaction;
    int idx_j = candidates[i] % M2_num_interaction;
    M1_pairs_potential[idx_i] = M2_pairs_potential[idx_j] = 1;
    int evaluate_score = 0;    
    if((M1_interaction[idx_i][2] > 0 && M2_interaction[idx_j][2] > 0) ||
      (M1_interaction[idx_i][2] > 0 && M2_interaction[idx_j][2] == 0) ||
      (M1_interaction[idx_i][2] == 0 && M2_interaction[idx_j][2] > 0))  {
      //  match of two base pairs, hbond may present
      //cout << "matching basepair" << endl;
      evaluate_score = scoring_function.weight_isosteric 
        * scoring_function.match_isosteric_basepair(M1_interaction[idx_i][2], M2_interaction[idx_j][2]);
      //cout << "basepair type:  " << M1_interaction[idx_i][2] << "  " << M2_interaction[idx_j][2] << endl;
      //cout << "score:  " << scoring_function.match_isosteric_basepair(M1_interaction[idx_i][2], M2_interaction[idx_j][2]) << endl;
    }  else if((M1_interaction[idx_i][2] < 0 && M2_interaction[idx_j][2] < 0) ||
      (M1_interaction[idx_i][2] < 0 && M2_interaction[idx_j][2] == 0) ||
      (M1_interaction[idx_i][2] == 0 && M2_interaction[idx_j][2] < 0))  {
      //  match of two non-adjancent stackings. hbond may present
      evaluate_score = scoring_function.weight_nonadjacent_stacking
        * scoring_function.match_stacking_basepair(M1_interaction[idx_i][2], M2_interaction[idx_j][2]);
      //cout << "basepair type:  " << M1_interaction[idx_i][2] << "  " << M2_interaction[idx_j][2] << endl;
      //cout << "score:  " << scoring_function.match_stacking_basepair(M1_interaction[idx_i][2], M2_interaction[idx_j][2]) << endl;
    }  else if(M1_interaction[idx_i][2] == 0 && M2_interaction[idx_j][2] == 0)  {
      //  match of two hydrogen bonds
      evaluate_score = scoring_function.weight_isosteric * 
        (scoring_function.hbond_match_base + scoring_function.hbond_match_bonus);
      //cout << "3:  " << scoring_function.hbond_match_bonus << endl;
    }
    //  no need to add nucleotide match, as it is already taken into account by Compute_alignment_score_max
    if(evaluate_score > M1_max_score[idx_i])  {
      M1_max_score[idx_i] = evaluate_score;
    }
    if(evaluate_score > M2_max_score[idx_j])  {
      M2_max_score[idx_j] = evaluate_score;
    }
  }
  int expected_matches = (sqrt(1 + 8 * num_edges_candidates) + 1) / 2;  
  //cout << "num_edges_candidates:  " << num_edges_candidates << endl;
  //cout << "min_degree:  " << min_degree << endl;
  //cout << "expected_matches:  " << expected_matches << endl;
  
  //vector<int> boundary;
  //find_aligned_boundary(current_clique, boundary);
  vector<int> M1_nuc_check, M2_nuc_check;
  get_aligned_nucs(current_clique, M1_nuc_check, M2_nuc_check);
  //cout << "Record_alignment:  stuck here 1" << endl;
  list<AlignedRegionType> M1_regions, M2_regions;
  get_aligned_regions(M1_nuc_check, M1_break_points, M1_regions);
  get_aligned_regions(M2_nuc_check, M2_break_points, M2_regions);
  
  int potential_pair_query = 0, potential_stacking_query = 0;
  int potential_pair_target = 0, potential_stacking_target = 0;
  
  CountDeletedInteractions(
      M1_pairs_potential, M1_interaction, M1_regions, potential_pair_query, potential_stacking_query
  );
  CountDeletedInteractions(
      M2_pairs_potential, M2_interaction, M2_regions, potential_pair_target, potential_stacking_target
  );
  
  evaluate_score -= scoring_function.weight_isosteric * (-1200 * potential_pair_query)
      + scoring_function.weight_isosteric * (-800 * potential_stacking_query);
  evaluate_score -= scoring_function.weight_isosteric * (-120 * potential_pair_target)
      + scoring_function.weight_isosteric * (-80 * potential_stacking_target);
  //  sort the max_scores
  for(i = 0; i < (int) M1_max_score.size() - 1; ++ i)  {
    for(j = i + 1; j < (int) M1_max_score.size(); ++ j)  {
      if(M1_max_score[i] < M1_max_score[j]){
        int temp = M1_max_score[i];
        M1_max_score[i] = M1_max_score[j];
        M1_max_score[j] = temp;
      }
    }
  }
  for(i = 0; i < (int) M2_max_score.size() - 1; ++ i)  {
    for(j = i + 1; j < (int) M2_max_score.size(); ++ j)  {
      if(M2_max_score[i] < M2_max_score[j]){
        int temp = M2_max_score[i];
        M2_max_score[i] = M2_max_score[j];
        M2_max_score[j] = temp;
      }
    }
  }
  int M1_score_bound = (int) M1_max_score.size() < expected_matches ? (int) M1_max_score.size() : expected_matches;
  int M2_score_bound = (int) M2_max_score.size() < expected_matches ? (int) M2_max_score.size() : expected_matches;
  int M1_sum_mScore = 0, M2_sum_mScore = 0;
  for(i = 0; i < M1_score_bound; ++ i)  {
    M1_sum_mScore += M1_max_score[i];
  }
  for(i = 0; i < M2_score_bound; ++ i)  {
    M2_sum_mScore += M2_max_score[i];
  }
  //int evaluate_score_T = evaluate_score + M1_sum_mScore_T < M2_sum_mScore_T ? M1_sum_mScore_T : M2_sum_mScore_T;
  evaluate_score += M1_sum_mScore < M2_sum_mScore ? M1_sum_mScore : M2_sum_mScore;
  //cout << "evaluate_score:  " << evaluate_score_T << "  " << evaluate_score << endl;
  return evaluate_score;
}

float MotifGraphMatching::chebyshev_inequality()  {
  float estimate_pvalue = 1;
  if(consider_pvalue)  {
    if(max_score > mean)  {
      estimate_pvalue = 1 / (((max_score - mean) / std) * ((max_score - mean) / std));
    } else  {
      estimate_pvalue = 1;
    }
  }
  return estimate_pvalue;
}

void MotifGraphMatching::print_alignment_search(void)  {
  if(lower_bound > max_score)  {
    return;
  }  else  {
    float normalized_score = (float) max_score / 10000.0;
    if(consider_pvalue)  {
      cout << M1_info << "  " << M2_info << "  " << normalized_score << "  " << pvalue << endl;
    }  else  {
      cout << M1_info << "  " << M2_info << "  " << normalized_score << endl;
    }
  }
  return;
}

void MotifGraphMatching::get_aligned_nucs(
    const std::vector<int>& clique, 
    std::vector<int>& nucs_first, 
    std::vector<int>& nucs_second
) {
  unsigned int i, j;
  //cout << "get_aligned_nucs:  clique size:  " << clique.size() << endl; 
  //  construct the nucleotides that particitate in selected interactions
  for(i = 0; i < clique.size(); ++ i)  {
    //  the two base interactions idx_i and idx_j match
    int idx_i = clique[i] / M2_num_interaction;
    int idx_j = clique[i] % M2_num_interaction;
    //  copy the matched indices
    nucs_first.push_back(M1_interaction[idx_i][0]);
    nucs_first.push_back(M1_interaction[idx_i][1]);
    nucs_second.push_back(M2_interaction[idx_j][0]);
    nucs_second.push_back(M2_interaction[idx_j][1]);
  }
  //cout << "get_aligned_nucs:  stuck here 1" << endl;  
  //  assert that the sizes of the matched nucleotides are the same
  if(nucs_first.size() != nucs_second.size())  {
    cout << "The sizes of matched nucleotides do not match: MotifGraphMatching::compute_alignment_score" << endl;
    return;
  }
  if(nucs_first.size() == 0)  {
    return;
  }
  //cout << nucs_first.size() << "  " << nucs_second.size() << endl;
  //  sort the matched nucleotides
  for(i = 0; i < nucs_first.size() - 1; ++ i)  {
    for(j = i + 1; j < nucs_first.size(); ++ j)  {
      if(nucs_first[i] > nucs_first[j])  {
        int temp = nucs_first[i];
        nucs_first[i] = nucs_first[j];
        nucs_first[j] = temp;
      }
    }
  }
  //cout << "get_aligned_nucs:  stuck here 2" << endl;
  for(i = 0; i < nucs_second.size() - 1; ++ i)  {
    for(j = i + 1; j < nucs_second.size(); ++ j)  {
      if(nucs_second[i] > nucs_second[j])  {
        int temp = nucs_second[i];
        nucs_second[i] = nucs_second[j];
        nucs_second[j] = temp;
      }
    }
  }
  //cout << "get_aligned_nucs:  before return" << endl;  
  return;
}

void MotifGraphMatching::get_aligned_regions(
    const std::vector<int>& aligned_nucs,
    const std::vector<int>& break_points,
    std::list<AlignedRegionType>& aligned_regions
) {
  /*
  cout << "********get_aligned_region called********" << endl;
  for(auto it = aligned_nucs.begin(); it != aligned_nucs.end(); ++ it) {
    cout << *it << "  ";
  }
  cout << endl;
  for(auto it = break_points.begin(); it != break_points.end(); ++ it) {
    cout << *it << "  ";
  }
  cout << endl;
  */
  //cout << "output:" << endl;  
  
  int i, j, k;
  int prev_end = -1;
  for(i = 0; i < (int) break_points.size(); ++ i) {
    for(j = prev_end + 1; j < (int) aligned_nucs.size(); ++ j) {
      if(break_points[i] <= aligned_nucs[j])  {
        int new_region_start = prev_end + 1;
        int new_region_end = j - 1;
        if(new_region_end >= new_region_start)  {
          AlignedRegionType single_region;
          single_region.begin = aligned_nucs[new_region_start];
          single_region.end = aligned_nucs[new_region_end];
          //cout << single_region.begin << "  " << single_region.end << endl; 
          for(k = new_region_start; k <= new_region_end; ++ k) {
            single_region.nucs.push_back(aligned_nucs[k]);
          }
          aligned_regions.push_back(single_region);
        }
        prev_end = new_region_end;
        break;
      }
    }
  }
  
  if(prev_end < (int) aligned_nucs.size() - 1)  {
    AlignedRegionType single_region;
    single_region.begin = aligned_nucs[prev_end + 1];
    single_region.end = aligned_nucs[aligned_nucs.size() - 1];
    //cout << single_region.begin << "  " << single_region.end << endl; 
    for(k = prev_end + 1; k <= (int) aligned_nucs.size() - 1; ++ k) {
      single_region.nucs.push_back(aligned_nucs[k]);
    }
    aligned_regions.push_back(single_region);
  }
  return;
}

void MotifGraphMatching::gen_alignment(
    const std::list<AlignedRegionType>& ref_aligned_regions,
    const std::vector<int>& ref_aligned_nucs,
    const std::vector<int>& target_aligned_nucs,
    std::string& ref_alignment,
    std::string& target_alignment
) {
  assert(ref_aligned_nucs.size() == target_aligned_nucs.size());
  //assert(*(-- ref_aligned_nucs.end()) == (-- ref_aligned_regions.end())->end);
  if(ref_aligned_nucs.size() == 0)  {
    return;
  }
  unsigned int i;
  auto region_it = ref_aligned_regions.begin();
  for(i = 0; i < ref_aligned_nucs.size() - 1; ++ i) {
    //cout << i << endl;
    if(ref_aligned_nucs[i] != ref_aligned_nucs[i + 1] && target_aligned_nucs[i] != target_aligned_nucs[i + 1])  {
      //  add the alignment of the current nucleotides
      ref_alignment.append(1, M1_seq[ref_aligned_nucs[i]]);
      target_alignment.append(1, M2_seq[target_aligned_nucs[i]]);
      if(i < ref_aligned_nucs.size() - 1 && ref_aligned_nucs[i] < region_it->end)  {
        string hold_alignment_1, hold_alignment_2;
        loop_align(
            ref_aligned_nucs[i], target_aligned_nucs[i], 
            hold_alignment_1, hold_alignment_2, 
            ref_aligned_nucs[i + 1], target_aligned_nucs[i + 1]
        );
        //  concaternate alignment
        ref_alignment += hold_alignment_1;
        target_alignment += hold_alignment_2;
      } else if(i < ref_aligned_nucs.size() - 1 && ref_aligned_nucs[i] >= region_it->end) {
        ref_alignment += "...";
        target_alignment += "...";
        ++ region_it;
      }
    }
  }
  ref_alignment.append(1, M1_seq[ref_aligned_nucs[i]]);
  target_alignment.append(1, M2_seq[target_aligned_nucs[i]]);
  return;
}

std::string MotifGraphMatching::record_alignment(const vector<int>& match_clique)  {
  //  output the alignment
  //  check the size of clique
  //cout << "Record_alignment:  stuck here 0" << endl;
  ostringstream record_out;
  int i;
  // gets the aligned key nucleotides
  vector<int> M1_nuc, M2_nuc;
  get_aligned_nucs(match_clique, M1_nuc, M2_nuc);
  //cout << "Record_alignment:  stuck here 1" << endl;
  list<AlignedRegionType> M1_regions, M2_regions;
  get_aligned_regions(M1_nuc, M1_break_points, M1_regions);
  get_aligned_regions(M2_nuc, M2_break_points, M2_regions);
  //cout << "Record_alignment:  stuck here 2" << endl;  
  string M1_alignment, M2_alignment;
  gen_alignment(M1_regions, M1_nuc, M2_nuc, M1_alignment, M2_alignment);
  //cout << "Record_alignment:  stuck here 3" << endl;
  if(M1_alignment.length() != M2_alignment.length())  {
    return "Error!\n";
  }
  string M1_aligned_region = spell_aligned_region(M1_regions, M1_hash_PDBID);
  string M2_aligned_region = spell_aligned_region(M2_regions, M2_hash_PDBID);
  
  float normalized_score = (float) max_score / 10000.0;
  
  record_out << "Query:  " << M1_info << endl;
  record_out << "Target: " << M2_info << endl; 
  record_out << "Score:  " << normalized_score << endl;
  if(consider_pvalue)  {
    record_out << std::setprecision(3) << "Pvalue: " << pvalue << endl;
  }  else  {
    record_out << "Pvalue: NA" << endl;
  }
  record_out << "Aligned region(Query):  " << M1_aligned_region << endl;
  record_out << "Aligned region(Target): " << M2_aligned_region << endl << endl;
  record_out << "        " << M1_alignment << endl;
  record_out << "        " << M2_alignment << endl << endl;
  
  record_out << "Aligned base pairs: " << endl;
  for(i = 0; i < (int) match_clique.size(); ++ i)  {
    //  the two base interactions idx_i and idx_j match
    int idx_i = match_clique[i] / M2_num_interaction;
    int idx_j = match_clique[i] % M2_num_interaction;
    string pairing_out = spell_pairing(idx_i, idx_j);
    if(!pairing_out.empty())  {
      record_out << pairing_out << endl;
    }
  }
  
  record_out << "Aligned base stackings: " << endl;
  for(i = 0; i < (int) match_clique.size(); ++ i)  {
    //  the two base interactions idx_i and idx_j match
    int idx_i = match_clique[i] / M2_num_interaction;
    int idx_j = match_clique[i] % M2_num_interaction;
    string stacking_out = spell_stacking(idx_i, idx_j);
    if(!stacking_out.empty())  {
      record_out << stacking_out << endl;
    }
  }
  return record_out.str();
}

void MotifGraphMatching::print_alignment(const vector<int>& match_clique)  {
  
  //  output the alignment
  //  check the size of clique
  //cout << "#########################################" << endl;
  int i;
  // gets the aligned key nucleotides
  vector<int> M1_nuc, M2_nuc;
  get_aligned_nucs(match_clique, M1_nuc, M2_nuc);
  list<AlignedRegionType> M1_regions, M2_regions;
  get_aligned_regions(M1_nuc, M1_break_points, M1_regions);
  get_aligned_regions(M2_nuc, M2_break_points, M2_regions);
  string M1_alignment, M2_alignment;
  //cout << "End of print" << endl;
  gen_alignment(M1_regions, M1_nuc, M2_nuc, M1_alignment, M2_alignment);
  //cout << "Good of print" << endl;
  if(M1_alignment.length() != M2_alignment.length())  {
    cout << "Error: alignment lengths differ!!!" << endl;
  }
  string M1_aligned_region = spell_aligned_region(M1_regions, M1_hash_PDBID);
  string M2_aligned_region = spell_aligned_region(M2_regions, M2_hash_PDBID);
  
  float normalized_score = (float) max_score / 10000.0;
  cout << endl;
  cout << "#  Aligning " << M1_info << " and " << M2_info << ":" << endl; 
  cout << "#  Alignment score:  " << normalized_score << endl;
  if(consider_pvalue)  {
    printf("#  P-value:  %.3e\n", pvalue);
  }  else  {
    cout << "#  P-value:  NA" << endl;
  }
  cout << "#  Query aligned region: " << M1_aligned_region << endl;
  cout << "#  Target aligned region: " << M2_aligned_region << endl << endl;
  cout << "\t" << M1_alignment << endl;
  cout << "\t" << M2_alignment << endl << endl;
  
  cout << "#  Matched base-pairing interactions: " << endl;
  for(i = 0; i < (int) match_clique.size(); ++ i)  {
    //  the two base interactions idx_i and idx_j match
    int idx_i = match_clique[i] / M2_num_interaction;
    int idx_j = match_clique[i] % M2_num_interaction;
    string pairing_out = spell_pairing(idx_i, idx_j);
    if(!pairing_out.empty())  {
      cout << pairing_out << endl;
    }
  }
  
  cout << "#  Matched base-stacking interactions: " << endl;
  for(i = 0; i < (int) match_clique.size(); ++ i)  {
    //  the two base interactions idx_i and idx_j match
    int idx_i = match_clique[i] / M2_num_interaction;
    int idx_j = match_clique[i] % M2_num_interaction;
    string stacking_out = spell_stacking(idx_i, idx_j);
    if(!stacking_out.empty())  {
      cout << stacking_out << endl;
    }
  }
  
  /*
  cout << evaluate_score << endl;
  */
  return;
}

std::string MotifGraphMatching::spell_aligned_region(
    const std::list<AlignedRegionType>& regions, 
    std::unordered_map<int, PDBIDType> hash_ID
) {
  ostringstream aligned_region;
  for(auto it = regions.begin(); it != regions.end(); ++ it) {
    aligned_region << it->begin;
    auto it_hash = hash_ID.find(it->begin);
    if(it_hash != hash_ID.end())  {
      aligned_region << ":" << it_hash->second;
    }
    aligned_region << "-";
    aligned_region << it->end;
    it_hash = hash_ID.find(it->end);
    if(it_hash != hash_ID.end())  {
      aligned_region << ":" << it_hash->second;
    }
    if(it != (-- regions.end()))  {
      aligned_region << ",";
    }
  }
  return aligned_region.str();
}

std::string MotifGraphMatching::spell_pairing(const unsigned int idx_i, const unsigned int idx_j) {
  string info;
  ostringstream pairing_info_1;
  ostringstream pairing_info_2;
  if(M1_interaction[idx_i][2] >= 0 && M2_interaction[idx_j][2] >= 0)  {
    // print the first matched pair
    PDBIDType m1_begin_alias;
    PDBIDType m1_end_alias;
    auto m1_begin_it = M1_hash_PDBID.find(M1_interaction[idx_i][0]);
    auto m1_end_it = M1_hash_PDBID.find(M1_interaction[idx_i][1]);
    if(m1_begin_it != M1_hash_PDBID.end())  {
      m1_begin_alias = m1_begin_it->second;
    }
    if(m1_end_it != M1_hash_PDBID.end())  {
      m1_end_alias = m1_end_it->second;
    }
    pairing_info_1 << "<" << M1_seq[M1_interaction[idx_i][0]]; 
    pairing_info_1 << "(" << M1_interaction[idx_i][0] << ":" << m1_begin_alias << ")"; 
    pairing_info_1 << "-" << M1_seq[M1_interaction[idx_i][1]]; 
    pairing_info_1 << "(" << M1_interaction[idx_i][1] << ":" << m1_end_alias << ")";
    if(M1_interaction[idx_i][5] == 1)  {
      pairing_info_1 << " cis ";
    }  else if(M1_interaction[idx_i][5] == 0)  {
      pairing_info_1 << " trans ";
    }  else if(M1_interaction[idx_i][5] == -1)  {
      pairing_info_1 << " hbond ";
    }
    char ca = (char) M1_interaction[idx_i][3];
    char cb = (char) M1_interaction[idx_i][4];
    pairing_info_1 << ca << "/" << cb << ">";
    // print the second matched pair
    PDBIDType m2_begin_alias;
    PDBIDType m2_end_alias;
    auto m2_begin_it = M2_hash_PDBID.find(M2_interaction[idx_j][0]);
    auto m2_end_it = M2_hash_PDBID.find(M2_interaction[idx_j][1]);
    if(m2_begin_it != M2_hash_PDBID.end())  {
      m2_begin_alias = m2_begin_it->second;
    }
    if(m2_end_it != M2_hash_PDBID.end())  {
      m2_end_alias = m2_end_it->second;
    }
    pairing_info_2 << "<" << M2_seq[M2_interaction[idx_j][0]];
    pairing_info_2 << "(" << M2_interaction[idx_j][0] << ":" << m2_begin_alias << ")";
    pairing_info_2 << "-" << M2_seq[M2_interaction[idx_j][1]];
    pairing_info_2 << "(" << M2_interaction[idx_j][1] << ":" << m2_end_alias << ")" ;
    if(M2_interaction[idx_j][5] == 1)  {
      pairing_info_2 << " cis ";
    }  else if(M2_interaction[idx_j][5] == 0)  {
      pairing_info_2 << " trans ";
    }  else if(M2_interaction[idx_j][5] == -1)  {
      pairing_info_2 << " hbond ";
    }
    ca = (char) M2_interaction[idx_j][3];
    cb = (char) M2_interaction[idx_j][4];
    pairing_info_2 << ca << "/" << cb << ">";
    info = "\t" + pairing_info_1.str() + string(45 - pairing_info_1.str().size(), ' ');
    info += "  MATCHES\t";
    info += pairing_info_2.str();
  }
  return info;
}  

std::string MotifGraphMatching::spell_stacking(const unsigned int idx_i, const unsigned int idx_j) {
  string info;
  ostringstream stacking_info_1;
  ostringstream stacking_info_2;
  if((M1_interaction[idx_i][2] < 0 && M2_interaction[idx_j][2] < 0) ||
        (M1_interaction[idx_i][2] < 0 && M2_interaction[idx_j][2] == 0) ||
        (M1_interaction[idx_i][2] == 0 && M2_interaction[idx_j][2] < 0))  {
    // print stacking in the first motif
    PDBIDType m1_begin_alias;
    PDBIDType m1_end_alias;
    auto m1_begin_it = M1_hash_PDBID.find(M1_interaction[idx_i][0]);
    auto m1_end_it = M1_hash_PDBID.find(M1_interaction[idx_i][1]);
    if(m1_begin_it != M1_hash_PDBID.end())  {
      m1_begin_alias = m1_begin_it->second;
    }
    if(m1_end_it != M1_hash_PDBID.end())  {
      m1_end_alias = m1_end_it->second;
    }
    stacking_info_1 << "<" << M1_seq[M1_interaction[idx_i][0]];
    stacking_info_1 << "(" << M1_interaction[idx_i][0] << ":" << m1_begin_alias << ")"; 
    stacking_info_1 << "-" << M1_seq[M1_interaction[idx_i][1]];
    stacking_info_1<< "(" << M1_interaction[idx_i][1] << ":" << m1_end_alias << ")";
    if(M1_interaction[idx_i][2] == -1)  {
      stacking_info_1 << " upward>";
    }  else if(M1_interaction[idx_i][2] == -2)  {
      stacking_info_1 << " downward>";
    }  else if(M1_interaction[idx_i][2] == -3)  {
      stacking_info_1 << " inward>";
    }  else if(M1_interaction[idx_i][2] == -4)  {
      stacking_info_1 << " outward>";
    }  else if(M1_interaction[idx_i][2] == 0)  {
      stacking_info_1 << " hbond>";
    }
    // print stacking in the second motif
    PDBIDType m2_begin_alias;
    PDBIDType m2_end_alias;
    auto m2_begin_it = M2_hash_PDBID.find(M2_interaction[idx_j][0]);
    auto m2_end_it = M2_hash_PDBID.find(M2_interaction[idx_j][1]);
    if(m2_begin_it != M2_hash_PDBID.end())  {
      m2_begin_alias = m2_begin_it->second;
    }
    if(m2_end_it != M2_hash_PDBID.end())  {
      m2_end_alias = m2_end_it->second;
    }
    stacking_info_2 << "<" << M2_seq[M2_interaction[idx_j][0]];
    stacking_info_2 << "(" << M2_interaction[idx_j][0] << ":" << m2_begin_alias << ")";
    stacking_info_2 << "-" << M2_seq[M2_interaction[idx_j][1]];
    stacking_info_2 << "(" << M2_interaction[idx_j][1] << ":" << m2_end_alias << ")";
    if(M2_interaction[idx_j][2] == -1)  {
      stacking_info_2 << " upward>";
    }  else if(M2_interaction[idx_j][2] == -2)  {
      stacking_info_2 << " downward>";
    }  else if(M2_interaction[idx_j][2] == -3)  {
      stacking_info_2 << " inward>";
    }  else if(M2_interaction[idx_j][2] == -4)  {
      stacking_info_2 << " outward>";
    }  else if(M2_interaction[idx_j][2] == 0)  {
      stacking_info_2 << " hbond>";
    }
    info = "\t" + stacking_info_1.str() + string(45 - stacking_info_1.str().size(), ' ');
    info += "  MATCHES\t";
    info += stacking_info_2.str();
  }  
  return info;
}
