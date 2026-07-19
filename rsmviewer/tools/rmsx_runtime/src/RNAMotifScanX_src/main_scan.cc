#include "structural_motif.h"
#include "annotated_structure.h"
#include "motif_graph_matching.h"
#include "simulation.h"
#include "parameters.h"
#include "map_PDB_index.h"
#include "scan_structure.h"

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

static std::string query_file = "";
static std::string structure_file = "";
static std::string map_PDBID_file = "";
static int extend_size = 5;
static int max_num_strands = 2;
static int max_closing_pairs = 1;
static double p_value_cutoff = 0.03;
static int num_threads = 1;

void print_help_message() {
  std::cout << "Usage: scan_structure query_motif RNA_structure" << std::endl;
  std::cout << "  [--map_pdb --extend_size --max_num_strands --max_closing_pairs --p_value_cutoff]" << std::endl; 
  std::cout << "  [--isosteric_matrix --stacking_matrix --nucleotide_matrix --write_alignment]" << std::endl;
  std::cout << std::endl;
  return;
}

int main(int argc, char** argv)  {

  boost_po::options_description desc("Allowed options");
  boost_po::positional_options_description p;
  p.add("query_motif", 1);
  p.add("structure", 1);
  
  desc.add_options()
    ("help", "print help message")
    ("fast", "use heuristic alignment algorithm")
    ("write_alignment", "output complete alignments")
    ("structure", boost_po::value<std::string>(&structure_file), "source RNA structure")
    ("query_motif", boost_po::value<std::string>(&query_file), "query motif whose non-canonical pairs are used as anchors")
    ("map_pdb", boost_po::value<std::string>(&map_PDBID_file), "mapping residue ID to original PDB ID")
    ("extend_size", boost_po::value<int>(&extend_size)->default_value(5), "size of extension to search for motifs")
    ("max_num_strands", boost_po::value<int>(&max_num_strands)->default_value(3), "maximum number of strands allowed")
    ("max_closing_pairs", boost_po::value<int>(&max_closing_pairs)->default_value(1), "maximum number of pairs allowed to close a motif")
    ("pvalue", boost_po::value<double>(&p_value_cutoff), "p-value cutoff for motif identification")
    ("num_threads", boost_po::value<int>(&num_threads)->default_value(1), "number of concurrent threads")
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
  
  if(!vm.count("structure") || !vm.count("query_motif"))  {
    std::cout << "Query motif or target RNA structure is not specified!" << std::endl;
    print_help_message();
    std::cout << desc << std::endl;
    exit(0);
  }
  
  bool use_heuristic = false;
  if(vm.count("fast")) use_heuristic = true;
  
  // define scoring functions
  InitParameters();
  Parameters scoring_function(default_parameters);
  if(vm.count("write_alignment"))  {
    ScanStructure search_motifs(
        query_file, structure_file, map_PDBID_file, scoring_function, 
        max_num_strands, extend_size, max_closing_pairs, p_value_cutoff,
        "./Results", num_threads, true, use_heuristic
    );
    search_motifs.Scan();
  }  else  {
    ScanStructure search_motifs(
        query_file, structure_file, map_PDBID_file, scoring_function, 
        max_num_strands, extend_size, max_closing_pairs, p_value_cutoff,
        "./Results", num_threads, false, use_heuristic
    );
    search_motifs.Scan();
  }  
  return 0;
}

