#include "structural_motif.h"
#include "parameters.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>

#ifndef _MOTIF_GRAPH_MATCHING_
#define _MOTIF_GRAPH_MATCHING_

#ifndef MAX
#define MAX 30000
#endif

struct AlignedRegionType  {
  int begin, end;
  std::vector<int> nucs;
};
struct MatchType  {
  int score;
  int num_matching;
  std::vector<int> motif_1;
  std::vector<int> motif_2;
};

class MotifGraphMatching  {
 public:
  std::string M1_seq, M2_seq;
  std::string M1_info, M2_info;
  int M1_num_interaction, M2_num_interaction;
  int M1_num_pairing, M2_num_pairing;
  int M1_num_stacking, M2_num_stacking;
  //  store breakpoints
  std::vector<int> M1_break_points;
  std::vector<int> M2_break_points;
  //  used to store adjacent stackings
  std::vector<int> M1_chain_stacking;
  std::vector<int> M2_chain_stacking;
  //  used to store non-adjacent interactions
  std::vector<std::vector<int> > M1_interaction;
  std::vector<std::vector<int> > M2_interaction;
  //  adjacency matrix for non-adjacent interactions
  std::vector<std::vector<int> > M1_interaction_adjacency;
  std::vector<std::vector<int> > M2_interaction_adjacency;
  //  a graph indicates whether matching of two pairs of interactions are compatible
  std::vector<std::vector <int> > interaction_compatibility;  
  //  penalty for deletion of nucleotides
  std::vector<std::vector<int> > M1_nuc_del_penalty;
  std::vector<std::vector<int> > M2_nuc_del_penalty;
  //  effective stacking orientation
  std::vector<std::vector<int> > M1_effective_stacking;
  std::vector<std::vector<int> > M2_effective_stacking;
  std::unordered_map<int, PDBIDType> M1_hash_PDBID;
  std::unordered_map<int, PDBIDType> M2_hash_PDBID;
  int lower_bound;
  int max_score;
  //  store a maximal clique found by a heuristic algorithm
  std::vector<int> maximal_clique;
  Parameters scoring_function;
  //  indicate whether output p-value and use it as a lower bound
  bool consider_pvalue;
  bool consider_cutoff;
  float mean, std;
  float pvalue_cutoff;
  float pvalue;
  /////////////////////////////////////////////////////////
  //  function definitions begin               //
  /////////////////////////////////////////////////////////
  MotifGraphMatching(
      const StructuralMotif& M1, 
      const StructuralMotif& M2, 
      const Parameters& input_scoring_function
  );
  MotifGraphMatching(
      const StructuralMotif& M1, 
      const StructuralMotif& M2, 
      const Parameters& input_scoring_function, 
      const float input_mean, const float input_std
  );
  MotifGraphMatching(
      const StructuralMotif& M1, 
      const StructuralMotif& M2, 
      const Parameters& input_scoring_function, 
      const float input_mean, const float input_std, 
      const float input_cutoff
  );
  //MotifGraphMatching(structural_motif M1, structural_motif M2, parameters input_scoring_function, 
  //  float input_mean, float input_std, float pvalue_cutoff);
  ~MotifGraphMatching();
  double GetAlignmentScore(void)  {
    return (static_cast<double>(max_score) / 10000.0);
  }
  double GetPvalue(void)  {
    return pvalue;
  }
  //  a wrapper function that computes the alignment and prints the results
  void align(void);
  void align_heuristic();
  void align_simulation();
  
  //  compute the maximal possible sequence alignment
  int maximal_seq_score(int begin, int end, std::string M_seq, std::vector <int> M_chain_stacking);
  
  //  partially order the interactions
  void sort_interaction(void);
  //  formulate the graph matching problem into a clique finding problem
  void compute_interaction_compatibility(void);
  
  //  estimating the lower bound
  void find_maximal_clique_heuristic(void);
  
  //  estimate the significance of matching
  bool significant_match(std::vector<int> interaction_1, std::vector<int> interaction_2);
  bool is_isosteric(std::vector<int> interaction_1, std::vector<int> interaction_2);
  void find_aligned_boundary(const std::vector<int>& clique, std::vector<int>& boundary);
  
  //  finding the optimal interaction matching through clique formulation
  void find_maximal_clique_optimal(void);
  void extend_clique_recursive(const std::vector<int>& current_clique, const std::vector<int>& candidates);
  int asymmetric_penalty(std::vector<std::vector<int> >& matched, int p1, int p2);
  //  compute the nucleotide deletion penalty, filling up matrices M1_nuc_del_penalty and M2_nuc_del_penalty
  //  and evaluate adjacent stacking effect
  void compute_nuc_del(void);
  
  //  compute the optimal alignment of loop regions with the consideration of adjacent base stacking effects
  //  the return value is the optimal alignment score
  int loop_align(void);  //  used to estimate the upper_bound, computes the optimal alignment disregarding base interactions
  int loop_align(int b1, int b2, std::string& M1_alignment, std::string& M2_alignment, int e1, int e2);  // alignment of two segments
  int loop_align(int b1, int b2, std::string& M1_alignment, std::string& M2_alignemnt);  // alignment with no ending-bound
  int loop_align(std::string &M1_alignment, std::string &M2_alignment, int e1, int e2);  // alignment with no beginning-bound
  
  //  compute the alignment score given the base-interaction match 
  int compute_alignment_score(const std::vector<int>& estimated_clique);
  //int compute_alignment_score_final(vector<bool> estimated_clique);
  int compute_alignment_score_max(const std::vector<int>& estimated_clique);
  int compute_alignment_score_upper_bound(const std::vector<int>& current_clique, const std::vector<int>& candidates);
  
  //  output alignment
  void print_alignment(const std::vector<int>& match_clique);
  std::string record_alignment(const std::vector<int>& match_clique);
  void print_alignment_search(void);
  //  compute the p-value of lower_bound using chebyshev inequality
  float chebyshev_inequality();
  int pvalue_lowerbound(void);
  
  //  some printing functions for intermediate results
  void interpret_interaction_compatibility(void);
  void interpret_maximal_clique(void);
  void interpret_candidates(const std::vector<bool>& candidates);
  std::string spell_stacking(const unsigned int idx_i, const unsigned int idx_j);
  std::string spell_pairing(const unsigned int idx_i, const unsigned int idx_j);
  std::string spell_aligned_region(
      const std::list<AlignedRegionType>& regions, 
      std::unordered_map<int, PDBIDType> hash_ID
  );
  void get_aligned_nucs(
      const std::vector<int>& clique, 
      std::vector<int>& nucs_first, 
      std::vector<int>& nucs_second
  );
  void get_aligned_regions(
      const std::vector<int>& aligned_nucs,
      const std::vector<int>& break_points,
      std::list<AlignedRegionType>& aligned_regions
  );
  void gen_alignment(
      const std::list<AlignedRegionType>& ref_aligned_regions,
      const std::vector<int>& ref_aligned_nucs,
      const std::vector<int>& target_aligned_nucs,
      std::string& ref_alignment,
      std::string& target_alignment
  );
  template <typename T>void ReduceSizeUnsorted(std::vector<T>& array_to_reduce, const int& index_to_delete);
  void CountDeletedInteractions(
      const std::vector<int>& aligned_pairs, const std::vector<std::vector<int> >& interactions, 
      const std::list<AlignedRegionType>& aligned_region, int& num_pair_deleted, int& num_stack_deleted
  );
};

#endif
