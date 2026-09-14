#include "motif_graph_matching.h"

#ifndef _ALIGN_FUNCTOR
#define _ALIGN_FUNCTOR

class AlignFunctor  {
 private:
  bool heuristic_;
  StructuralMotif *query_;
  StructuralMotif *target_;
  Parameters *scoring_function_;
  float mean_;
  float std_;
  float pvalue_;
  AlignmentInfoType* results_;
 public:
  AlignFunctor(
      bool in_heuristic, StructuralMotif *in_query, StructuralMotif *in_target,
      Parameters *in_scoring_function, float in_mean, float in_std, float in_pvalue,
      AlignmentInfoType *in_results
  ) : heuristic_(in_heuristic), 
      query_(in_query),
      target_(in_target),
      scoring_function_(in_scoring_function),  
      mean_(in_mean),
      std_(in_std),
      pvalue_(in_pvalue),
      results_(in_results)
  {
    return;
  }
  ~AlignFunctor() {
    return;
  }
  void run()  {
    MotifGraphMatching align_motifs(
        *query_, *target_, *scoring_function_, mean_, std_, pvalue_
    );
    if(!heuristic_)  {
      align_motifs.align();
    }  else  {
      align_motifs.align_heuristic();
    }
    results_->motif_info = align_motifs.M2_info;
    results_->align_score = align_motifs.GetAlignmentScore();
    results_->p_value = align_motifs.GetPvalue();
    results_->maximal_clique = align_motifs.maximal_clique;
    std::vector<int> foo;
    align_motifs.get_aligned_nucs(align_motifs.maximal_clique, foo, results_->aligned_nucs);
    std::list<AlignedRegionType> idv_nucs;
    align_motifs.get_aligned_regions(results_->aligned_nucs, align_motifs.M2_break_points, idv_nucs);
    results_->aligned_region_info = align_motifs.spell_aligned_region(idv_nucs, align_motifs.M2_hash_PDBID);
    return;
  }
};

#endif
