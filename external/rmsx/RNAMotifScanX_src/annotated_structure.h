#include "structural_motif.h"
#include "map_PDB_index.h"

#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <map>
#include <cmath>
#include <unordered_map>
#include <tuple>
#include <list>

#ifndef _Annotated_Structure_H_
#define _Annotated_Structure_H_

typedef std::pair<int, int> StrandType;
typedef std::vector<int> MotifBoundType;

class AnnotatedStructure: public StructuralMotif	{
 
 public:
	
	AnnotatedStructure(void);
	AnnotatedStructure(const std::string& file_name);
	AnnotatedStructure(const std::string& file_name, const std::unordered_map<int, PDBIDType>& in_hash_PDBID);
	~AnnotatedStructure(void);

	void BuildShortestPath(const unsigned int dist_cutoff);
	void RemoveStackingPairs(const int max_closing_pairs); // this function removes nested canonical
	    // base pairs, leaving only the boundary base pairs (at most "min_boundary_pairs" at each end)
	// when given the query motif and applies anchoring
	void PredictLocalMotifs(
	    StructuralMotif& query, 
	    const int dist, 
	    const int num_strands
	);  
	// when no query motif is given and extract all non-canoincal pairs
	void PredictLocalMotifs(const int dist, const int num_strands);  
  void PrintAllInstances(void);
  void PermuteStrandOrientation(
    const MotifBoundType& in_motif, 
    std::list<MotifBoundType>& permuted_motifs
  );
  void GenerateMotif(MotifBoundType included_strands, StructuralMotif& motif);
  void CopyMotifInstances(std::list<MotifBoundType>& motif_copy);
 protected:
	std::vector< std::vector<int> > shortest_path_matrix_;
	std::unordered_map<int, std::list<int> > interaction_adjacency_list_;
	using StructuralMotif::is_noncanonical;
	std::unordered_map<int, int> chain_stacking_hash_;
	std::list<MotifBoundType> motif_instance;
	void build_adjacency_list(void);
	std::unordered_map<int, PDBIDType> hash_PDBID_;
	
	inline void shortest_path_single(const unsigned int dist_cutoff, const unsigned int node);
	inline int adjacency_graph_dist(const unsigned int node_1, const unsigned int node_2);
	// outputting anchored base pairs in "anchors", each entry contsin 4 integer values
	// which are the start and end of base pair in query motif, and the start and end of base pair in target structure
	void get_matched_basepairs(StructuralMotif& query, std::list< std::vector<int> >& anchors);
	void get_local_motif(
	    int start, int end, int num_key_strands,
	    unsigned int dist, unsigned int num_strands
	);
	void get_local_motif(
	    std::vector<int>& anchor_info, const int& num_key_strands, 
	    const int& allowed_gap, const int& num_strands
	);
	// returns a list of adjacent nucleotide within dist
	void get_near_nucleotides(int pivot, unsigned int dist, std::list<int>& adjacent_nucs);
	void get_near_nucleotides(
    const int& start, const int& start_up, const int& start_end, 
    const int& end, const int& end_up, const int& end_end,
    const int& allowed_gap, std::list<int>& adjacent_nucs
  );
	// takes in a list of nucleotides, and output the pruned strands
	// each pruned strands should be closed by a base pair, and should not have more than two
	// consecutive WC pairs close them.
	// TODO: write this function
	void prune_strands(
	    int start, int end,            // the start and end location of the key base pair
	    int num_key_strands,
	    unsigned int num_strands,      // number of strands allowed for each motif
	    const std::list<int>& adjacent_nucs,  // list of nucleotides who may constitute the motif 
	    std::list<std::pair<int, int> >& pruned_strands
	);
	void define_strands(const std::list<int>& nucs, std::list<std::pair<int, int> >& strands);
	void print_interaction(const std::list<int>& included_nucs);
	void print_interaction(const std::list<std::pair<int, int> >& strands);
	void add_element(int size, std::vector<int> permuted_elements, std::list<std::vector<int> >& orders);
	void update_motif_instances(const MotifBoundType single_motif);
	bool complete_covered(MotifBoundType motif_a, MotifBoundType motif_b);
	MotifBoundType convert_strand_to_bound_array(const std::list<std::pair<int, int> >& strands);
  void all_permutation(int size, std::list<std::vector<int> >& orders);
  void remove_dangling(std::list<std::pair<int, int> >& strands);
  void remove_excess_strands(
      const int& start, const int& end, const unsigned int& num_strands, 
      std::list<std::pair<int, int> >& strands
  );
  int count_interactions_betwee_strands(const std::pair<int, int>& a, const std::pair<int, int>& b);
  void split_sticky_stands(
      const int& start, const int& end, std::list<std::pair<int, int> >& strands
  );
};

#endif

