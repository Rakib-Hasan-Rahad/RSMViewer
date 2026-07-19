#include "structural_motif.h"
#include "annotated_structure.h"
#include "motif_graph_matching.h"
#include "simulation.h"
#include "parameters.h"
#include "map_PDB_index.h"

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/thread.hpp>
#include <boost/threadpool.hpp>
#include <boost/bind.hpp>
#include <boost/smart_ptr.hpp>
#include <iostream>
#include <iomanip>
#include <list>
#include <string>
#include <vector>
#include <cstdlib>

#ifndef _SCAN_STRUCTURE_H_
#define _SCAN_STRUCTURE_H_

struct AlignmentInfoType  {
  std::string motif_info;
  int list_index;
  double align_score;
  double p_value;
  std::string print_info;
  std::string aligned_region_info;
  std::vector<int> aligned_nucs;
  std::vector<int> maximal_clique;
  StructuralMotif *target;
};

class ScanStructure {
 public:
  ScanStructure();
  ScanStructure(
      const std::string& in_query_file, const std::string& in_structure_file,
      const std::string& in_hash_PDBID_file,
      const Parameters& in_scoring_function,
      const unsigned int& in_num_strands,
      const unsigned int& in_extend_size, const unsigned int& in_closing_pairs,
      const double& in_pvalue_cutoff,  
      const std::string& in_out_directory,
      const int in_num_threads,
      const bool in_write_alignment,
      const bool in_heuristic
  );
  ~ScanStructure();
  void Scan(void);
 protected:
  std::string query_file_;
  StructuralMotif query_motif_;
  std::string structure_file_;
  std::string hash_PDBID_file_;
  Parameters scoring_function_;
  unsigned int num_strands_;
  unsigned int extend_size_;
  unsigned int num_closing_pairs_;
  double pvalue_cutoff_;
  int num_threads_;
  bool write_alignment_;
  bool heuristic_;
  std::string out_directory_;
  std::list<AlignmentInfoType> high_score_hits_;
  
  void AlignAll(void);
  void CleanUp(void);
  void WriteResultsTabulated(void);
  void WriteResultsComplete(void);
  void RemoveDuplicate(void);
  void MotifInfoToHash(std::string motif_info, std::unordered_map<int, int>& included_nucs);
  bool SignificantOverlap(const AlignmentInfoType& info_1, const AlignmentInfoType& info_2);
};

#endif
