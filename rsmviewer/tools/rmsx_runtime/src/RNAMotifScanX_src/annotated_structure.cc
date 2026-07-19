#include "annotated_structure.h"

using namespace std;

AnnotatedStructure::AnnotatedStructure(void) 
: StructuralMotif::StructuralMotif()  {
  return;
}

AnnotatedStructure::AnnotatedStructure(const string& file_name) 
: StructuralMotif::StructuralMotif(file_name)  {
  // allocate space for the shortest_path matrix
  unsigned int i, j;
  shortest_path_matrix_.resize(sequence_.length());
  for(i = 0; i < sequence_.length(); ++ i)  {
    shortest_path_matrix_[i] = vector<int> (sequence_.length(), MAX);
    // set distance to itself to be 0
    shortest_path_matrix_[i][i] = 0;
    // set distance to adjacent nucleotides to 1
    if(i > 0)  {
      shortest_path_matrix_[i][i - 1] = 1;
    } 
    if(i < sequence_.length() - 1)  {
      shortest_path_matrix_[i][i + 1] = 1;
    }
    // set distance to interacting nucleotides to 1
    for(j = 0; j < sequence_.length(); ++ j)  {
      if(interaction_adjacency_[i][j] >= 0 && interaction_adjacency_[i][j] != MAX)  {
        // considering hydrogen bonding and base-pairing
        shortest_path_matrix_[i][j] = 1;
      }
    }
  }
  
  for(i = 0; i < chain_stacking_.size(); ++ i)  {
    if(chain_stacking_[i] != MAX)  {
      chain_stacking_hash_[i] = chain_stacking_[i];
    }
  }
  build_adjacency_list();
  return;
}

AnnotatedStructure::AnnotatedStructure(
    const std::string& file_name, 
    const std::unordered_map<int, PDBIDType>& in_hash_PDBID
): StructuralMotif::StructuralMotif(file_name) {
  // allocate space for the shortest_path matrix
  unsigned int i, j;
  hash_PDBID_ = in_hash_PDBID;
  shortest_path_matrix_.resize(sequence_.length());
  for(i = 0; i < sequence_.length(); ++ i)  {
    shortest_path_matrix_[i] = vector<int> (sequence_.length(), MAX);
    // set distance to itself to be 0
    shortest_path_matrix_[i][i] = 0;
    // set distance to adjacent nucleotides to 1
    if(i > 0)  {
      shortest_path_matrix_[i][i - 1] = 1;
    } 
    if(i < sequence_.length() - 1)  {
      shortest_path_matrix_[i][i + 1] = 1;
    }
    // set distance to interacting nucleotides to 1
    for(j = 0; j < sequence_.length(); ++ j)  {
      if(interaction_adjacency_[i][j] >= 0 && interaction_adjacency_[i][j] != MAX)  {
      // considering hydrogen bonding and base-pairing
        shortest_path_matrix_[i][j] = 1;
      }
    }
  }
  
  for(i = 0; i < chain_stacking_.size(); ++ i)  {
    if(chain_stacking_[i] != MAX)  {
      chain_stacking_hash_[i] = chain_stacking_[i];
    }
  }
  build_adjacency_list();
  return;
}

AnnotatedStructure::~AnnotatedStructure()  {
  return;
}

void AnnotatedStructure::RemoveStackingPairs(const int max_closing_pairs)  {
  /*
  for(auto it = record_interaction_.begin(); it != record_interaction_.end(); ++ it) {
    cout << "original:  " << (*it)[0] << "  " << (*it)[1] << endl;
  }
  */
  unsigned int i, j, k;
  for(i = 0; i < interaction_adjacency_.size() - 1; ++ i) {
    for(j = i + 1; j < interaction_adjacency_[i].size(); ++ j) {
      //cout << i << "\t" << j << " : " << interaction_adjacency_[i][j] << endl;
      if(is_canonical(interaction_adjacency_[i][j]))  {
        // detect the longest span of this stack
        //cout << "hello" << i << "\t" << j << endl;
        unsigned int m = 1;
        while(i + m < sequence_.length() && j - m < sequence_.length() 
            && is_canonical(interaction_adjacency_[i + m][j - m])) {
           //cout << "m:  " << m << endl;
           ++ m;
        }
        // if m is long enough the remove the inner nested canoncial pairs
        if((int) m > 2 * max_closing_pairs)  {
          for(k = max_closing_pairs; k < m - max_closing_pairs; ++ k) {
            interaction_adjacency_[i + k][j - k] = interaction_adjacency_[j - k][i + k] = MAX;
            //cout << "to delete: "  << i + k << "\t" << j - k << endl;  
          }
        }
      }
    }
  }
  //cout << "good here" << endl
  
  vector<vector<int> > updated_interactions;
  for(auto it = record_interaction_.begin(); it != record_interaction_.end(); ++ it) {
    //cout << "KKKKK: " << (*it)[0] << "  " << (*it)[1] << endl;
    if(interaction_adjacency_[(*it)[0]][(*it)[1]] != MAX)  {
      updated_interactions.push_back(*it);
    }
  }
  record_interaction_ = updated_interactions;
  /*
  for(auto it = record_interaction_.begin(); it != record_interaction_.end(); ++ it) {
    cout << "kept:  " << (*it)[0] << "  " << (*it)[1] << endl;
  }
  */
  return;
}

void AnnotatedStructure::get_near_nucleotides(int pivot, unsigned int dist, list<int>& adjacent_nucs)  {
  assert(pivot >= 0 && pivot < (int) sequence_.length());
  assert(dist > 0);
  int i;
  for(i = 0; i < (int) sequence_.length(); ++ i)  {
    if(shortest_path_matrix_[pivot][i] <= (int) dist)  {
      adjacent_nucs.push_back(i);
    }
  }
  return;
}

void AnnotatedStructure::get_near_nucleotides(
    const int& start, const int& start_up, const int& start_end, 
    const int& end, const int& end_up, const int& end_end,
    const int& allowed_gap, list<int>& adjacent_nucs)  {
  // define reachable nucleotides
  unordered_map<int, bool> recruited_nucs;
  //cout << "index s1: " << start - start_up << "  " << start + start_end << endl;
  //cout << "index s2: " << end - end_up << "  " << end + end_end << endl;
  
  for(int i = start - start_up - allowed_gap; i <= start + start_end + allowed_gap; ++ i) {
    recruited_nucs[i] = true;
    for(auto it = interaction_adjacency_list_[i].begin(); it != interaction_adjacency_list_[i].end(); ++ it) {
      recruited_nucs[*it] = true;
    }
  }
  for(int i = end - end_up - allowed_gap; i <= end + end_end + allowed_gap; ++ i) {
    recruited_nucs[i] = true;
    for(auto it = interaction_adjacency_list_[i].begin(); it != interaction_adjacency_list_[i].end(); ++ it) {
      recruited_nucs[*it] = true;
    }
  }
  //cout << "initial recruited: " << recruited_nucs.size() << endl;
  // record the reachable nucleotides
  for(auto it = recruited_nucs.begin(); it != recruited_nucs.end(); ++ it) {
  //  cout << "**:  " << it->first << endl;
    if(it->first >= 0 && it->first < (int) sequence_.length())  {
      adjacent_nucs.push_back(it->first);
    }
  }
  //cout << "end of function" << endl;
  return;
}

// finds the local motifs that are defined by the base pair starting and ending at
// "start" and "end", respectively
void AnnotatedStructure::get_local_motif(
    int start,  // start of the base pair 
    int end,    // end of the base pair
    int num_key_strands, 
    unsigned int dist,  // distance cutoff to define a nucleotide as adjacent
    unsigned int num_strands // number of strands allowed in the local motif
)  {
  //cout << "get_local_motif: begin:  " << start << "  " << end << endl;
  list<int> adjacent_nucs;  // the nucleotides that are adjacent to the start 
  get_near_nucleotides(start, dist + 1, adjacent_nucs); // note that the
  //for(auto it = adjacent_nucs.begin(); it != adjacent_nucs.end(); ++ it) {
  //  cout << *it << "  ";
  //}
  //cout << endl;
  list<pair<int, int> > pruned_strands;
  prune_strands(start, end, num_key_strands, num_strands, adjacent_nucs, pruned_strands);
  //cout << "***************" << endl;
  //for(auto it = pruned_strands.begin(); it != pruned_strands.end(); ++ it) {
  //  cout << it->first << "  " << it->second << endl;
  //}
  MotifBoundType single_instance = convert_strand_to_bound_array(pruned_strands);
  update_motif_instances(single_instance);
  return;
}

void AnnotatedStructure::get_local_motif(
    std::vector<int>& anchor_info, const int& num_key_strands,
    const int& allowed_gap, const int& num_strands
) {
  //cout << "********* anchors: " << anchor_info[2] << "  " << anchor_info[3] << endl;
  //cout << "begin of get local motif with " << num_key_strands << endl;
  list<int> adjacent_nucs;  // the nucleotides that are adjacent to the start 
  get_near_nucleotides(
      anchor_info[2], anchor_info[4], anchor_info[5], 
      anchor_info[3], anchor_info[6], anchor_info[7],
      allowed_gap, adjacent_nucs
  );
  /* 
  for(auto it = adjacent_nucs.begin(); it != adjacent_nucs.end(); ++ it) {
    cout << *it << "  ";
  }
  cout << endl;
  */
  list<pair<int, int> > pruned_strands;
  prune_strands(anchor_info[2], anchor_info[3], num_key_strands, num_strands, adjacent_nucs, pruned_strands);
  //cout << "***************" << endl;
  //for(auto it = pruned_strands.begin(); it != pruned_strands.end(); ++ it) {
  //  cout << it->first << "  " << it->second << endl;
  //}
  //cout << "prune strands done" << endl;
  MotifBoundType single_instance = convert_strand_to_bound_array(pruned_strands);
  //cout << "convert strands done" << endl;
  update_motif_instances(single_instance);
  //cout << "update motif done" << endl;
  return;
}

void AnnotatedStructure::PrintAllInstances(void)  {
  for(auto it = motif_instance.begin(); it != motif_instance.end(); ++ it) {
    list<MotifBoundType> permuted_motifs;
    PermuteStrandOrientation(*it, permuted_motifs);
    for(auto it_list = permuted_motifs.begin(); it_list != permuted_motifs.end(); ++ it_list) {
      for(auto it_vector = (*it_list).begin(); it_vector != (*it_list).end(); ++ it_vector) {
        cout << *it_vector << " ";
      }
      cout << endl;
    }
  }
}

bool AnnotatedStructure::complete_covered(MotifBoundType motif_a, MotifBoundType motif_b) {
  // testing if motif_a is covered by motif_b
  unsigned int i, j;
  bool is_covered = true;
  for(i = 0; i < motif_a.size(); i += 2) {
    bool single_strand_covered = false;
    for(j = 0; j < motif_b.size(); j += 2) {
      if(motif_b[j] <= motif_a[i] && motif_b[j + 1] >= motif_a[i + 1])  {
        single_strand_covered = true;
        break;
      }
    }
    if(!single_strand_covered)  {
      is_covered = false;
      break;
    }
  }
  return is_covered;
}

MotifBoundType AnnotatedStructure::convert_strand_to_bound_array(const list<pair<int, int> >& strands)  {
  MotifBoundType bound_array;
  for(auto it = strands.begin(); it != strands.end(); ++ it) {
    bound_array.push_back(it->first);
    bound_array.push_back(it->second);
  }
  return bound_array;
}

void AnnotatedStructure::update_motif_instances(const MotifBoundType single_motif) {
  bool has_included = false;
  for(auto it = motif_instance.begin(); it != motif_instance.end(); ++ it) {
    if(complete_covered(single_motif, *it))  {
      // replace the old motif
      has_included = true;
    } else if(complete_covered(*it, single_motif)) {
      *it = single_motif;
      has_included = true;
    }
  }
  if(!has_included)  {
    motif_instance.push_back(single_motif);
  }
  return;
}

void AnnotatedStructure::add_element(int size, vector<int> permuted_elements, list<vector<int> >& orders)  {
  int i, j;
  if((int) permuted_elements.size() >= size)  {
    orders.push_back(permuted_elements);
    return;
  }
  vector<int> to_include;
  for(i = 0; i < size; ++ i) {
    bool visited = false;
    for(j = 0; j < (int) permuted_elements.size(); ++ j) {
      if(permuted_elements[j] == i)  {
        visited = true;
        break;
      } 
    }
    if(!visited)  {
      to_include.push_back(i);
    }
  }
  
  for(auto it = to_include.begin(); it != to_include.end(); ++ it) {
    vector<int> updated_elements = permuted_elements;
    updated_elements.push_back(*it);
    if((int) updated_elements.size() >= size)  {
      orders.push_back(updated_elements);
    } else  {
      add_element(size, updated_elements, orders);
    }
  }  

  return;
}

void AnnotatedStructure::all_permutation(int size, list<vector<int> >& orders)  {
  assert(size > 0);
  int i;
  for(i = 0; i < size; ++ i) {
    vector<int> permuted_elements;
    permuted_elements.push_back(i);
    add_element(size, permuted_elements, orders);
  }
  return;
}

void AnnotatedStructure::print_interaction(const list<int>& included_nucs)  {
  for(auto it_i = included_nucs.begin(); it_i != included_nucs.end(); ++ it_i) {
    for(auto it_j = it_i; it_j != included_nucs.end(); ++ it_j) {
      if(interaction_adjacency_[*it_i][*it_j] != MAX)  {
        cout << *it_i << "  " << *it_j << ":  " << interaction_adjacency_[*it_i][*it_j] << endl;
      }
    }
  } 
  return;
}

void AnnotatedStructure::print_interaction(const std::list<std::pair<int, int> >& strands)  {
  list<int> included_nucs;
  for(auto it = strands.begin(); it != strands.end(); ++ it) {
    for(int i = it->first; i <= it->second; ++ i) {
      included_nucs.push_back(i);
    }
  }    
  int num_interactions = 0;
  for(auto it_i = included_nucs.begin(); it_i != included_nucs.end(); ++ it_i) {
    for(auto it_j = it_i; it_j != included_nucs.end(); ++ it_j) {
      if(interaction_adjacency_[*it_i][*it_j] != MAX)  {
        cout << *it_i << "  " << *it_j << ":  " << interaction_adjacency_[*it_i][*it_j] << endl;
        ++ num_interactions;
      }
    }
  } 
  return;
}

bool _cmp_nucs(const int& a, const int& b)  {
  if(a < b) return true;
  return false;
}

void AnnotatedStructure::define_strands(const list<int>& nucs, list<pair<int, int> >& strands) {
  if(nucs.size() <= 0)  {
    return;
  }
  list<int> nucs_copy = nucs;
  nucs_copy.sort(_cmp_nucs);
  int prev, strand_start;
  prev = strand_start = nucs_copy.front();
  auto it = ++ nucs_copy.begin();
  while(it != nucs_copy.end()) {
    if(*it > prev && *it > prev + 2)  { // allowing 2-nucleotide gap
      strands.push_back({strand_start, prev});
      strand_start = *it; 
    }
    prev = *it;
    ++ it;
  }
  strands.push_back({strand_start, nucs_copy.back()});
  return;
}

int AnnotatedStructure::count_interactions_betwee_strands(
    const std::pair<int, int>& a, const std::pair<int, int>& b
) {
  int num_interactions = 0;
  // if the two strands have overlap, return 0
  if(!(a.first > b.second || a.second < b.first))  {
    return 0;
  }
  for(int i = a.first; i <= a.second; ++ i) {
    for(int j = b.first; j <= b.second; ++ j) {
      // interaction_adjacency_[i][j] != MAX indicates that it is an interaction
      num_interactions += (int) (interaction_adjacency_[i][j] != MAX);
    }
  }
  return num_interactions;
}

bool _cmp_num_interactions(
    const pair<list<pair<int, int> >::iterator, int >& a, 
    const pair<list<pair<int, int> >::iterator, int >& b
)  {
  if((a.second > b.second) || 
      (a.second == b.second && a.first->second - a.first->first > b.first->second - b.first->first))  {
    return true;
  }  
  return false;
}


void AnnotatedStructure::remove_excess_strands(
    const int& start, const int& end, const unsigned int& num_strands, 
    std::list<std::pair<int, int> >& strands
)  {
  // identify the key strands
  if(strands.size() <= num_strands)  {
    return;
  }
  list<pair<int, int> > key_strands;
  list<list<pair<int, int> >::iterator> to_delete;
  for(auto it = strands.begin(); it != strands.end(); ++ it) {
    if((it->first <= start && it->second >= start) || (it->first <= end && it->second >= end))  {
      key_strands.push_back(*it);
      to_delete.push_back(it);
    }
  }
  for(auto it = to_delete.begin(); it != to_delete.end(); ++ it) {
    strands.erase(*it);
  }
  // identify strands with most interactions
  list<pair<list<pair<int, int> >::iterator, int > > interaction_map;
  for(auto it = strands.begin(); it != strands.end(); ++ it) {
    int num_interactions = 0;
    for(auto it_k = key_strands.begin(); it_k != key_strands.end(); ++ it_k) {
      num_interactions += count_interactions_betwee_strands(*it, *it_k);
    }
    interaction_map.push_back({it, num_interactions});
  }
  interaction_map.sort(_cmp_num_interactions);
  // include strands 
  while(key_strands.size() < num_strands && !interaction_map.empty()) {
    key_strands.push_back(*((interaction_map.front()).first));
    interaction_map.pop_front();
  }
  strands = key_strands;
  //exit(0);
  return;
}

void AnnotatedStructure::remove_dangling(std::list<std::pair<int, int> >& strands)  {
  // remove dangling nucleotides that do not participate in any interaction
  unordered_map<int, bool> all_nucs;
  unordered_map<int, bool> with_interaction;
  for(auto it = strands.begin(); it != strands.end(); ++ it) {
    for(int i = it->first; i <= it->second; ++ i) {
      all_nucs[i] = true;
    }
  }
  for(auto it_i = all_nucs.begin(); it_i != all_nucs.end(); ++ it_i) {
    for(auto it_j = all_nucs.begin(); it_j != all_nucs.end(); ++ it_j) {
      if(interaction_adjacency_[it_i->first][it_j->first] != MAX)  {
        with_interaction[it_i->first] = true;
        with_interaction[it_j->first] = true;
      }
    }
  }
  auto it = strands.begin();
  list<list<pair<int, int> >::iterator> to_delete;
  while(it != strands.end()) {
    //cout << it->first << "  " << it->second << endl;
    while(with_interaction.find(it->first) == with_interaction.end() && it->first < it->second) {
      ++ it->first;
    }
    while(with_interaction.find(it->second) == with_interaction.end() && it->first < it->second) {
      -- it->second;
    }
    if(it->second - it->first <= 0)  {
      to_delete.push_back(it);
    }
    ++ it;
  }
  for(auto it = to_delete.begin(); it != to_delete.end(); ++ it) {
    strands.erase(*it);
  }
  return;
}

void AnnotatedStructure::split_sticky_stands(
    const int& start, const int& end, std::list<std::pair<int, int> >& strands
) {
  //cout << "BEGIN of split strands:  " << start << " " << end << endl;
  //cout << "!!!!! Input" << endl;
  //cout << "num strands: " << strands.size() << endl;
  //for(auto it = strands.begin(); it != strands.end(); ++ it) {
  //  cout << it->first << "  " << it->second << endl;
  //}
  list<pair<int, int> > key_strands;
  list<list<pair<int, int> >::iterator> to_delete;
  for(auto it = strands.begin(); it != strands.end(); ++ it) {
    if((it->first <= start && it->second >= start) || (it->first <= end && it->second >= end))  {
      key_strands.push_back(*it);
      to_delete.push_back(it);
    }
  }
  // if the two nucleotides are contained in two separate strands, then return
  if(key_strands.size() > 1)  {
    return;
  }
  // otherwise we need to cut the strands in half
  
  unordered_map<int, int> cross_count;
  pair<int, int> s_range = key_strands.front();
  //cout << "SPLIT needed: " << s_range.first << "  " << s_range.second << endl;
  for(int i = s_range.first + 1; i <= s_range.second; ++ i) {
    cross_count[i] = 0;
  }
  //cout << "Init set:  " << cross_count.size() << endl;
  for(int i = s_range.first; i < s_range.second; ++ i) {
    for(int j = i + 1; j <= s_range.second; ++ j) {
      //cout << i << "  " << j << endl;
      if(interaction_adjacency_[i][j] != MAX)  {
        for(int k = i + 1; k <= j; ++ k) {
          ++ cross_count[k];
        }
      }
    }
  }
  //cout << "Cross count resolved" << endl;
  // get the breakpoint and cut the strand
  int max_break = (cross_count.begin())->first, max_cross = (cross_count.begin())->second;
  for(auto it = cross_count.begin(); it != cross_count.end(); ++ it) {
    if(it->second > max_cross && ((it->first > start && it->first <= end) || (it->first > end && it->first <= start)))  {
      max_cross = it->second;
      max_break = it->first;
    }
  }
  key_strands.front().second = max_break - 1;
  key_strands.push_back({max_break, s_range.second});
  for(auto it = to_delete.begin(); it != to_delete.end(); ++ it) {
    strands.erase(*it);
  }
  strands.splice(strands.end(), key_strands);
  //cout << "!!!!!" << endl;
  //for(auto it = strands.begin(); it != strands.end(); ++ it) {
  //  cout << it->first << "  " << it->second << endl;
  //}
  //cout << "END of split strands" << endl;
  return;
}

void AnnotatedStructure::prune_strands(
    int start, int end,            // the start and end location of the key base pair
    int num_key_strands,
    unsigned int num_strands,      // number of strands allowed for each motif
    const list<int>& adjacent_nucs,  // list of nucleotides who may constitute the motif 
    list<pair<int, int> >& pruned_strands
) {
  /*
  cout << "**********" << endl;
  for(auto it = adjacent_nucs.begin(); it != adjacent_nucs.end(); ++ it)  {
    cout << *it << ",";
  }
  cout << endl;
  */
  // identify the key strands
  list<int> original_nucs = adjacent_nucs;
  list<pair<int, int> > strands;  // each strand consists of its starting and ending point
  define_strands(original_nucs, strands);
  remove_dangling(strands);
  //cout << "*****************" << endl; 
  //for(auto it = strands.begin(); it != strands.end(); ++ it) {
  //  cout << "strands0: " << it->first << "  " << it->second << endl;
  //}
  if(strands.size() > 1)  {
    split_sticky_stands(start, end, strands);
  }
  //cout << "*****************" << endl; 
  //for(auto it = strands.begin(); it != strands.end(); ++ it) {
  //  cout << "strands1: " << it->first << "  " << it->second << endl;
  //}
  // check if the number of strands included exceeds the number of strands allowed
  if(strands.size() > num_strands)  {
    // remove the excessive strands which has little interaction with the core structure
    remove_excess_strands(start, end, num_strands, strands);
  }
  //cout << "*****************" << endl; 
  //for(auto it = strands.begin(); it != strands.end(); ++ it) {
  //  cout << "strands2: " << it->first << "  " << it->second << endl;
  //}
  remove_dangling(strands);
  //cout << "*****************" << endl; 
  //for(auto it = strands.begin(); it != strands.end(); ++ it) {
  //  cout << "strands3: " << it->first << "  " << it->second << endl;
  //}
  //remove_dangling(strands);
  //cout << "*****************" << endl; 
  //for(auto it = strands.begin(); it != strands.end(); ++ it) {
  //  cout << "strands4: " << it->first << "  " << it->second << endl;
  //}
  pruned_strands = strands;
  return;
}


// check whether node_1 and node_2 form any interaction in structure
// if yes, return 1, else return -MAX
inline int AnnotatedStructure::adjacency_graph_dist(const unsigned int node_1, const unsigned int node_2)  {
  assert(node_1 >= 0 && node_1 < sequence_.length());
  assert(node_2 >= 0 && node_2 < sequence_.length());  
  return shortest_path_matrix_[node_1][node_2];
}

// computes the shortest path between "node (as input argument)" and all other nodes
inline void AnnotatedStructure::shortest_path_single(const unsigned int dist_cutoff, const unsigned int node)  {
  
  assert(node >= 0 && node < sequence_.length());
  
  unsigned int i;
  unordered_map<int, int> dist_to_nodes; // first field is the node ID, 
                                         // second field is distance from the current node to the pivot node
  
  pair<int, int> min_dist = {MAX, 0};
  // set initial distance
  for(i = 0; i < sequence_.length(); ++ i)  {  
    dist_to_nodes[i] = adjacency_graph_dist(node, i);
    if(i != node && dist_to_nodes[i] < min_dist.first)  {
      min_dist.first = dist_to_nodes[i];
      min_dist.second = i;
    }
  }
  
  // apply Dijkstra's algorithm to compute the shortest path
  do  {    
    if(min_dist.first > (int) dist_cutoff || dist_to_nodes.size() <= 1)  {
      break;
    }

    pair<int, int> selected = min_dist;  // first field is the distance, second field is the node ID
    shortest_path_matrix_[node][selected.second] = selected.first;
    shortest_path_matrix_[selected.second][node] = selected.first;
    // delete the selected node
    dist_to_nodes.erase(selected.second);
    
    list<int>* connections = &interaction_adjacency_list_[selected.second];
    list<int> to_traverse(connections->begin(), connections->end());
    if(selected.second > 0)  {
      to_traverse.push_back(selected.second - 1);
    }
    if(selected.second < (int) sequence_.length() - 1)  {
       to_traverse.push_back(selected.second + 1);
    }
    for(auto it = to_traverse.begin(); it != to_traverse.end(); ++ it)  {
      // see if any path is shorter if we go through the selected node
      int t_node = *it;
      if(dist_to_nodes.find(t_node) != dist_to_nodes.end() &&
          dist_to_nodes[t_node] > selected.first + adjacency_graph_dist(selected.second, t_node))  {
        dist_to_nodes[t_node] = selected.first + adjacency_graph_dist(selected.second, t_node);
      }
    } 
    
    min_dist.first = dist_to_nodes.begin()->second;
    min_dist.second = dist_to_nodes.begin()->first;
    for(auto it = dist_to_nodes.begin(); it != dist_to_nodes.end(); ++ it)  {
      if(it->second < min_dist.first)  {
        min_dist.first = it->second;
        min_dist.second = it->first;
      }
    }    
  }  while(true);  
  return;
}

// gets a list of anchoring positions that matches at least one non-canoincal base pair in the query motif
// each element in the anchor is a vector of size 4, which correspond to:
// pair start in query, pair end in query, pair start in target, and pair end in target
void AnnotatedStructure::get_matched_basepairs(StructuralMotif& query, list< vector<int> >& anchors)  {  
  for(auto it_query = query.record_interaction_.begin(); it_query != query.record_interaction_.end(); ++ it_query)  {
    for(auto it_this = this->record_interaction_.begin(); it_this != this->record_interaction_.end(); ++ it_this)  {
      if(is_noncanonical(*it_query) && is_noncanonical(*it_this))  {
        
        if((*it_query)[5] == (*it_this)[5] && (*it_query)[3] == (*it_this)[3] && (*it_query)[4] == (*it_this)[4])  {
          vector<int> matched_anchor(4);
          matched_anchor = {
              (*it_query)[0], (*it_query)[1], (*it_this)[0], (*it_this)[1]
          }; 
          matched_anchor.push_back(query.GetUpDistance((*it_query)[0])); 
          matched_anchor.push_back(query.GetDownDistance((*it_query)[0]));
          matched_anchor.push_back(query.GetUpDistance((*it_query)[1])); 
          matched_anchor.push_back(query.GetDownDistance((*it_query)[1]));
          anchors.push_back(matched_anchor);
        }  else if((*it_query)[5] == (*it_this)[5] && (*it_query)[3] == (*it_this)[4] && (*it_query)[4] == (*it_this)[3])  {
          vector<int> matched_anchor(4);
          matched_anchor = {
              (*it_query)[0], (*it_query)[1], (*it_this)[1], (*it_this)[0]
          };
          matched_anchor.push_back(query.GetUpDistance((*it_query)[1])); 
          matched_anchor.push_back(query.GetDownDistance((*it_query)[1]));
          matched_anchor.push_back(query.GetUpDistance((*it_query)[0])); 
          matched_anchor.push_back(query.GetDownDistance((*it_query)[0]));
          anchors.push_back(matched_anchor);
        }  
      }
    }
  }  
  return;
}

// interface for building the shortest path
void AnnotatedStructure::BuildShortestPath(const unsigned int dist_cutoff)  {
  unsigned int i;
  /*
  for(i = 0; i < shortest_path_matrix_.size(); ++ i) {
    cout << "#######: " << i << endl;
    for(j = 0; j < shortest_path_matrix_[i].size(); ++ j) {
      cout << shortest_path_matrix_[i][j] << endl;
    }
    cout << endl;
  }
  */
  for(i = 0; i < sequence_.length(); ++ i)  {
    shortest_path_single(dist_cutoff, i);
  }
  return;  
}

// this function builds the adjacency list from the adjacency matrix
// the purpose is for fast identification of interactions
void AnnotatedStructure::build_adjacency_list(void) {
  unsigned int i, j;
  for(i = 0; i < sequence_.length() - 1; ++ i) {
    for(j = i + 1; j < sequence_.length(); ++ j) {
      if(interaction_adjacency_[i][j] != MAX)  {
        interaction_adjacency_list_[i].push_back(j);
        interaction_adjacency_list_[j].push_back(i);
      }
    }
  }
  return;
}

void AnnotatedStructure::CopyMotifInstances(list<MotifBoundType>& motif_copy) {
  for(auto it = motif_instance.begin(); it != motif_instance.end(); ++ it) {
    motif_copy.push_back(*it);
  }
  return;
}

void AnnotatedStructure::GenerateMotif(MotifBoundType included_strands, StructuralMotif& motif) {
  
  //cout << "BEGIN OF GENERATEMOTIF" << endl;
  int i, j;
  motif.file_info_ = this->file_info_ + ":";
  motif.sequence_ = "";
  map<int, int> index_hash;
  unordered_map<int, PDBIDType> motif_hash_PDB;
  unsigned int n = 0;
  char* num_holder_1 = new char[100];
  char* num_holder_2 = new char[100];
  //cout << "size of strands: " << included_strands.size() << endl;
  for(i = 0; i < floor(included_strands.size() / 2); ++ i) {
    //cout << "index: " << included_strands[2 * i] << " " << included_strands[2 * i + 1] << endl;
    sprintf(num_holder_1, "%d", included_strands[2 * i]);
    sprintf(num_holder_2, "%d", included_strands[2 * i + 1]);
    motif.file_info_ += string(num_holder_1) + "-" + string(num_holder_2);
    if(i < floor(included_strands.size() / 2) - 1)  {
      motif.file_info_ += "_";
    }
    motif.sequence_ += 
        this->sequence_.substr(
            included_strands[2 * i], 
            included_strands[2 * i + 1] - included_strands[2 * i] + 1
        );
    if(i < floor(included_strands.size() / 2) - 1)  {
      motif.break_points_.push_back(motif.sequence_.length());
    }
    for(j = included_strands[2 * i]; j <= included_strands[2 * i + 1]; ++ j) {
      index_hash[n] = j;
      auto it_map = hash_PDBID_.find(j);
      if(it_map != hash_PDBID_.end())  {
        motif_hash_PDB[n] = it_map->second;
      }
      ++ n;
    }
  }
  motif.hash_PDBID_ = motif_hash_PDB;
  delete [] num_holder_1;
  delete [] num_holder_2;
  //cout << "good here" << endl;
  //cout << "n: " << n << endl;
  //cout << motif.file_info_ << endl;
  //cout << motif.sequence_ << endl;
  //motif.break_points_.resize(motif.break_points_.size() - 1); // remove the last break point cause
      // it is the length of the entire sequence
  motif.interaction_adjacency_.resize(motif.sequence_.length());
  motif.chain_stacking_.resize(motif.sequence_.length(), MAX);
  for(i = 0; i < (int) motif.sequence_.length(); ++ i) {
    motif.interaction_adjacency_[i].resize(motif.sequence_.length(), MAX);
  }
  for(i = 0; i < (int) n - 1; ++ i) {
    auto it = chain_stacking_hash_.find(index_hash[i]);
    if(it != chain_stacking_hash_.end() && index_hash[i + 1] - index_hash[i] == 1)  {
      motif.chain_stacking_[i] = it->second;
    }
    for(j = i + 1; j < (int) n; ++ j) {
      motif.interaction_adjacency_[i][j] = this->interaction_adjacency_[index_hash[i]][index_hash[j]];
      motif.interaction_adjacency_[j][i] = this->interaction_adjacency_[index_hash[j]][index_hash[i]];
      if(this->interaction_adjacency_[index_hash[i]][index_hash[j]] >= 0 
          && this->interaction_adjacency_[index_hash[i]][index_hash[j]] != MAX)  {
        char edge_a, edge_b;
        bool is_cis;
        interpret_base_pair(this->interaction_adjacency_[index_hash[i]][index_hash[j]],
            edge_a, edge_b, is_cis);
        vector<int> single_interaction(6);
        single_interaction[0] = i;
        single_interaction[1] = j;
        single_interaction[2] = this->interaction_adjacency_[index_hash[i]][index_hash[j]];
        single_interaction[3] = (int) edge_a;
        single_interaction[4] = (int) edge_b;
        single_interaction[5] = (int) is_cis;
        motif.record_interaction_.push_back(single_interaction);
      } else if(this->interaction_adjacency_[index_hash[i]][index_hash[j]] < 0) {
        vector<int> single_interaction(6);
        single_interaction[0] = i;
        single_interaction[1] = j;
        single_interaction[2] = this->interaction_adjacency_[index_hash[i]][index_hash[j]];
        single_interaction[3] = reverse_hash_stacking(this->interaction_adjacency_[index_hash[i]][index_hash[j]]);
        motif.record_interaction_.push_back(single_interaction);
      }
    }
  }
  
  return;
}

void AnnotatedStructure::PermuteStrandOrientation(
    const MotifBoundType& in_motif, 
    list<MotifBoundType>& permuted_motifs
) {
  int num_strands = floor(in_motif.size() / 2);
  // generate strand permutations orders
  list<vector<int> > orders;
  all_permutation(num_strands, orders);
  // generate bound array for each permutation
  for(auto it = orders.begin(); it != orders.end(); ++ it) {
    MotifBoundType shuffled_motif;
    for(auto it_vector = it->begin(); it_vector != it->end(); ++ it_vector) {
      shuffled_motif.push_back(in_motif[2 * (*it_vector)]);
      shuffled_motif.push_back(in_motif[2 * (*it_vector) + 1]);
    }
    permuted_motifs.push_back(shuffled_motif);
  }
  return;
}

void AnnotatedStructure::PredictLocalMotifs(
    StructuralMotif& query,
    const int dist,
    const int num_strands
)  {
  list< vector<int> > anchors;
  get_matched_basepairs(query, anchors);
  for(auto it = anchors.begin(); it != anchors.end(); ++ it) {
    //get_local_motif((*it)[2], (*it)[3], dist, num_strands);
    // only if the base pair in the query are from the same strand and has the same 
    // mapping order in the target, we do not need to flip the query
    if(query.IsSameStrand((*it)[0], (*it)[1]) && 
        (((*it)[0] > (*it)[1] && (*it)[2] > (*it)[3]) || ((*it)[0] < (*it)[1] && (*it)[2] < (*it)[3])))  {
      get_local_motif(*it, 1, 0, num_strands);
    } else  {
      get_local_motif(*it, 2, 0, num_strands);
    }
  }
  //cout << "Predict local motifs done" << endl;
  return;
}

void AnnotatedStructure::PredictLocalMotifs(const int dist, const int num_strands) {
  unsigned int i, j;
  for(i = 0; i < sequence_.length() - 1; ++ i) {
    for(j = i + 1; j < sequence_.length(); ++ j) {
      if(is_noncanonical(interaction_adjacency_[i][j]))  {
        get_local_motif(i, j, 2, dist, num_strands);
      }
    }
  }
  return;
}


