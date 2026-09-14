#include "scan_structure.h"
#include "align_functor.h"

using namespace std;

ScanStructure::ScanStructure()  {
  return;
}

ScanStructure::~ScanStructure() {
  return;
}

ScanStructure::ScanStructure(
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
) {
  query_file_ = in_query_file;
  structure_file_ = in_structure_file;
  hash_PDBID_file_ = in_hash_PDBID_file;
  scoring_function_ = in_scoring_function;
  num_strands_ = in_num_strands;
  extend_size_ = in_extend_size;
  num_closing_pairs_ = in_closing_pairs;
  pvalue_cutoff_ = in_pvalue_cutoff;
  out_directory_ = in_out_directory;
  num_threads_ = in_num_threads;
  write_alignment_ = in_write_alignment;
  heuristic_ = in_heuristic;
  return;
}

void ScanStructure::Scan(void)  {
  AlignAll();
  //cout << "Finish alignment" << endl;
  CleanUp();
  //cout << "Finish cleanup" << endl;
  if(!write_alignment_)  {
    WriteResultsTabulated();
  }  else  {
    WriteResultsComplete();
  }
  //cout << "Finish writing" << endl;
  return;
}

void ScanStructure::AlignAll(void) {
  // construct the query motif
  StructuralMotif query(query_file_);
  query_motif_ = query;
  Simulation shuffle(query_file_, scoring_function_);
  shuffle.estimate_mean_std();
  // construct the large RNA structure
  std::unordered_map<int, PDBIDType> hash_PDBID;
  if(!hash_PDBID_file_.empty())  {
    MapPDBIndex index_ID(hash_PDBID_file_);
    index_ID.GetPDBMapping(hash_PDBID);
  }
  AnnotatedStructure whole_structure(structure_file_, hash_PDBID);
  whole_structure.RemoveStackingPairs(num_closing_pairs_);
  //whole_structure.BuildShortestPath(extend_size_);
  // predict local motifs in the large structure by taking non-canonical base pairs
  // in the query motif as anchors
  whole_structure.PredictLocalMotifs(query_motif_, extend_size_, num_strands_);
  list<MotifBoundType> target_list;
  whole_structure.CopyMotifInstances(target_list);
  list<AlignmentInfoType*> result_pool;
  boost::threadpool::pool alignment_pool(num_threads_);
  for(auto it = target_list.begin(); it != target_list.end(); ++ it) {
    std::list<MotifBoundType> permuted_targets;
    // permute all strand orientations
    whole_structure.PermuteStrandOrientation(*it, permuted_targets);
    for(auto it_permuted = permuted_targets.begin(); it_permuted != permuted_targets.end(); ++ it_permuted) {
      // do the alignment
      //cout << "New structure" << endl;
      //StructuralMotif *instance = new StructuralMotif;
      //whole_structure.GenerateMotif(*it_permuted, *instance);
      //cout << "instance:  " << instance->file_info_ << endl;
      //instance->PrintMotifInfo();
      AlignmentInfoType* aln_info = new AlignmentInfoType;
      aln_info->target = new StructuralMotif;
      whole_structure.GenerateMotif(*it_permuted, *(aln_info->target));
      result_pool.push_back(aln_info);
      //(result_pool.back()).target = new StructuralMotif;
      //whole_structure.GenerateMotif(*it_permuted, *(result_pool.back().target));
      // schedule multi-threaded run    
      boost::shared_ptr<AlignFunctor> job(
          new AlignFunctor(
              heuristic_, &query, aln_info->target, &scoring_function_, 
              shuffle.mean, shuffle.std, pvalue_cutoff_, result_pool.back()
          )
      );    
      boost::threadpool::schedule(alignment_pool, boost::bind(&AlignFunctor::run, job));
    }
  }
  alignment_pool.wait();
  for(auto it = result_pool.begin(); it != result_pool.end(); ++ it) {
    //cout << (*it)->motif_info << "  " << (*it)->align_score << endl;
    if((*it)->p_value <= pvalue_cutoff_)  {
      high_score_hits_.push_back(**it);
    }
  }
  return;
}

bool _cmp_motif_score(AlignmentInfoType s1, AlignmentInfoType s2)  {
  if(s1.align_score >= s2.align_score)  {
    return true;
  }  else  {
    return false;
  }
}

void ScanStructure::MotifInfoToHash(
  std::string motif_info, std::unordered_map<int, int>& included_nucs
) {
  vector<string> holder_1, holder_2;
  boost::split(holder_1, motif_info, boost::is_any_of(":"));
  boost::split(holder_2, holder_1[1], boost::is_any_of("_"));
  unsigned int i, j;
  int index = 0;
  for(i = 0; i < holder_2.size(); ++ i) {
    vector<string> holder;
    boost::split(holder, holder_2[i], boost::is_any_of("-"));
    unsigned int begin = atoi(holder[0].c_str());
    unsigned int end = atoi(holder[1].c_str());
    for(j = begin; j <= end; ++ j) {
      included_nucs[index] = j;
      ++ index;
    }
  }
  return;
}

bool ScanStructure::SignificantOverlap(
    const AlignmentInfoType& info_1, 
    const AlignmentInfoType& info_2
)  {
  unordered_map<int, int> nucs_1, nucs_2;
  MotifInfoToHash(info_1.motif_info, nucs_1);
  MotifInfoToHash(info_2.motif_info, nucs_2);
  //cout << "******************************************************" << endl;
  //cout << "IDs: " << info_1.motif_info << " " << info_2.motif_info << endl;
  unordered_map<int, bool> nucs_original_1, nucs_original_2;
  for(auto it = info_1.aligned_nucs.begin(); it != info_1.aligned_nucs.end(); ++ it) {
    if(nucs_1.find(*it) != nucs_1.end())  {
      //cout << "Get: " << *it << " " << nucs_1[*it] << endl;
      nucs_original_1[nucs_1[*it]] = true;
    }
  }
  for(auto it = info_2.aligned_nucs.begin(); it != info_2.aligned_nucs.end(); ++ it) {
    if(nucs_2.find(*it) != nucs_2.end())  {
      //cout << "Get: " << *it << " " << nucs_2[*it] << endl;
      nucs_original_2[nucs_2[*it]] = true;
    }
  }
  int num_overlap = 0;
  for(auto it = nucs_original_1.begin(); it != nucs_original_1.end(); ++ it) {
    if(nucs_original_2.find(it->first) != nucs_original_2.end())  {
      //cout << "Overlap nucs:  " << it->first << endl;
      ++ num_overlap;
    }
  }
  // compute the Jaccard similarity between sets
  if(((double) num_overlap) / ((double) (nucs_original_1.size() + nucs_original_2.size())) > 0.20)  {
    return true;
  } else  {
    return false;
  }
}

void ScanStructure::RemoveDuplicate(void) {
  list<list<AlignmentInfoType>::iterator> alignment_to_delete;
  unordered_map<int, bool> deleted_items;
  // assignning list index for each high-score motif
  int index = 0;
  for(auto it = high_score_hits_.begin(); it != high_score_hits_.end(); ++ it) {
    it->list_index = index ++;
  }
  for(auto it_i = high_score_hits_.begin(); it_i != prev(high_score_hits_.end()); ++ it_i) {
    //cout << "finish duplication checking: " << cout << it_i->motif_info << endl;
    if(deleted_items.find(it_i->list_index) != deleted_items.end())  {
      continue;
    }
    
    for(auto it_j = next(it_i); it_j != high_score_hits_.end(); ++ it_j) {
      if(SignificantOverlap(*it_i, *it_j))  {
       //cout << "item to delete:  " << it_j->motif_info << endl;
        //cout << it_i->motif_info << " " << it_j->motif_info << endl;
        deleted_items[it_j->list_index] = true;
      }
    }
  }
  
  for(auto it_i = high_score_hits_.begin(); it_i != high_score_hits_.end(); ++ it_i) {
    if(deleted_items.find(it_i->list_index) != deleted_items.end())  {
      alignment_to_delete.push_back(it_i);
    }
  }
  for(auto it_delete = alignment_to_delete.begin(); it_delete != alignment_to_delete.end(); ++ it_delete) {
    //cout << "item being deleted:  " << (*it_delete)->motif_info << endl;
    delete (*it_delete)->target;
    high_score_hits_.erase(*it_delete);
  }
  //for(auto it = high_score_hits_.begin(); it != high_score_hits_.end(); ++ it) {
  //  cout << "******:  " << it->motif_info << endl;
  //}
  return;
}

void ScanStructure::CleanUp(void) {
  high_score_hits_.sort(_cmp_motif_score);
  //cout << "finish sorting high-score hits" << endl;
  RemoveDuplicate();
  return;
}

void ScanStructure::WriteResultsTabulated(void)	{
  cout << "#fragment_ID\taligned_regions\talignment_score\tP-value" << endl;
  for(auto it_i = high_score_hits_.begin(); it_i != high_score_hits_.end(); ++ it_i) {
    cout << it_i->motif_info << "\t" << it_i->aligned_region_info << "\t" << it_i->align_score << "\t" << it_i->p_value << endl;
  }
  return;
}

void ScanStructure::WriteResultsComplete(void)  {
  for(auto it_i = high_score_hits_.begin(); it_i != high_score_hits_.end(); ++ it_i) {
    cout << "<<<<<<<<<<  ALIGNMENT BEGINS  <<<<<<<<<<" << endl;
    //cout << it_i->motif_info << " " << it_i->align_score << " " << it_i->p_value << endl;    
    MotifGraphMatching out_alignment(query_motif_, *(it_i->target), scoring_function_, 0, 0);
    out_alignment.align_simulation();
    out_alignment.max_score = it_i->align_score * 10000;
    out_alignment.pvalue = it_i->p_value;
    out_alignment.print_alignment(it_i->maximal_clique);
    cout << endl;
    cout << ">>>>>>>>>>   ALIGNMENT ENDS   >>>>>>>>>>" << endl << endl;;
  }
  return;
}
