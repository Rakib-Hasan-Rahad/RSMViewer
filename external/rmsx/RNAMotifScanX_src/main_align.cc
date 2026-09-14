#include "structural_motif.h"
#include "annotated_structure.h"
#include "motif_graph_matching.h"
#include "simulation.h"
#include "parameters.h"

#include <boost/program_options.hpp>
#include <iostream>
#include <list>
#include <string>
#include <vector>
#include <cstdlib>

static ParameterType default_parameters;

void InitParameters(void)  {
  char* home_dir_char = getenv("RNAMOTIFSCANX_PATH");
  if(home_dir_char == NULL)  {
    std::cerr << "Environmental variable $RNAMOTIFSCANX_PATH is not set!" << std::endl;
    std::cerr << "Please follow instructions in README." << std::endl;
    exit(0);
  }
  std::string home_dir(home_dir_char);
  default_parameters.isosteric_file_name = home_dir + "/mat/iso.mat";
  default_parameters.stacking_file_name = home_dir + "/mat/stk.mat";
  default_parameters.nucleotide_file_name = home_dir + "/mat/nuc.mat";
  default_parameters.gap_opening = -5.0;
  default_parameters.gap_extension = -30.0;
  default_parameters.triple_interaction_bonus = 3.0;
  default_parameters.hbond_match_base = 2.0;
  default_parameters.hbond_match_bonus = 2.0;
  default_parameters.weight_isosteric = 2.0;
  default_parameters.weight_nonadjacent_stacking = 1.0; 
  default_parameters.weight_adjacent_stacking = 0.2;
  default_parameters.weight_sequence = 0.1;
  default_parameters.pair_deletion = -12;
  default_parameters.stacking_deletion = -2;
  default_parameters.asym_nuc = -0.2; 
  default_parameters.asym_loop = -2; 
  default_parameters.cons_stacking = 0.5; 
  default_parameters.stack_to_bulge = -1.5;
  default_parameters.stack_to_internal_asym = -1;
  default_parameters.stack_to_internal_sym = -0.5;
  return;  
}

namespace boost_po = boost::program_options;

static std::string query_file;
static std::string target_file;

void print_help_message() {
  std::cout << "Usage: align query target" << std::endl;
  //std::cout << "  [--out_dir --query_motif --extend_size --max_num_strands --max_closing_pairs]" << std::endl;
  std::cout << std::endl;
  return;
}

/*
int main()  {
  InitParameters();
  Parameters scoring_function(default_parameters);

  std::string query_file = "/home/cczhong/Codes/RNAMotifScanX/Queries/k-turn_consensus.struct";
  std::string target_file = "/home/cczhong/Codes/RNAMotifScanX/WorkSpace/1S72_9:70-80_100-108.rmf";
  
  StructuralMotif query(query_file);
  StructuralMotif target(target_file);
  
  MotifGraphMatching align_motifs(query, target, scoring_function);
  align_motifs.align();
  align_motifs.print_alignment(align_motifs.maximal_clique);
  
  return 0;
}
*/



int main(int argc, char** argv)  {
  
  boost_po::options_description desc("Allowed options");
  boost_po::positional_options_description p;
  p.add("query_motif", 1);
  p.add("target_motif", 1);
  
  desc.add_options()
    ("help", "print help message")
    ("fast", "using heuristic alignment approach")
    ("query_motif", boost_po::value<std::string>(&query_file), "query motif whose non-canonical pairs are used as anchors")
    ("target_motif", boost_po::value<std::string>(&target_file), "source RNA structure")
    ("isosteric_matrix", boost_po::value<std::string>(&default_parameters.isosteric_file_name)->default_value("./Data/iso.mat"), "scoring matrix for base-pair substitution with isostericity")
    ("stacking_matrix", boost_po::value<std::string>(&default_parameters.stacking_file_name)->default_value("./Data/stk.mat"), "scoring matrix for base-stacking substitution")
    ("nucleotide_matrix", boost_po::value<std::string>(&default_parameters.nucleotide_file_name)->default_value("./Data/nuc.mat"), "scoring matrix for nucleotide substitution")
  ;
  boost_po::variables_map vm;
  boost_po::store(boost_po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
  boost_po::notify(vm);
  
  if(vm.count("help"))  {
    print_help_message();
    std::cout << desc << std::endl;
    exit(0);
  }
  
  InitParameters();
  Parameters scoring_function(default_parameters);
  //std::cout << default_parameters.gap_extension << std::endl;
  //std::cout << scoring_function.gap_extension << std::endl;
  if(vm.count("query_motif") && vm.count("target_motif"))  {  
    StructuralMotif query(query_file);
    StructuralMotif target(target_file);
    //target.PrintMotifInfo();
    MotifGraphMatching align_motifs(query, target, scoring_function);
    if(vm.count("fast"))  {
      align_motifs.align_heuristic();
    } else {
      align_motifs.align();
    }
    //align_motifs.interpret_maximal_clique();
    align_motifs.print_alignment(align_motifs.maximal_clique);
  } else  {
    std::cout << "Source RNA structure not specified!" << std::endl;
    print_help_message();
    std::cout << desc << std::endl;
    exit(0);
  }
  
  return 0;
}

