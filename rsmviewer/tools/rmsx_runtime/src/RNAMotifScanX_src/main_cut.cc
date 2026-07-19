#include "map_PDB_index.h"
#include "structural_motif.h"
#include "annotated_structure.h"

#include <boost/program_options.hpp>
#include <iostream>
#include <list>
#include <string>
#include <vector>
#include <cstdlib>

namespace boost_po = boost::program_options;

static std::string query_file;
static std::string structure_file;
static std::string map_PDBID_file;
static std::string output_directory;
static int extend_size = 5;
static int max_num_strands = 3;
static int max_closing_pairs = 1;

void print_help_message() {
  std::cout << "Usage: cut source_RNA_structure" << std::endl;
  std::cout << "  [--out_dir --query_motif --map_pdb --extend_size --max_num_strands --max_closing_pairs]" << std::endl;
  std::cout << std::endl;
  return;
}

int main(int argc, char** argv)  {
  boost_po::options_description desc("Allowed options");
  boost_po::positional_options_description p;
  p.add("structure", 1);
  p.add("out_dir", 1);
  
  desc.add_options()
    ("help", "print help message")
    ("structure", boost_po::value<std::string>(&structure_file), "source RNA structure")
    ("out_dir", boost_po::value<std::string>(&output_directory)->default_value("./"), "output directory for cut motifs")
    ("query_motif", boost_po::value<std::string>(&query_file), "query motif whose non-canonical pairs are used as anchors")
    ("map_pdb", boost_po::value<std::string>(&map_PDBID_file), "mapping residue ID to original PDB ID")
    ("extend_size", boost_po::value<int>(&extend_size)->default_value(5), "size of extension to search for motifs")
    ("max_num_strands", boost_po::value<int>(&max_num_strands)->default_value(3), "maximum number of strands allowed")
    ("max_closing_pairs", boost_po::value<int>(&max_closing_pairs)->default_value(1), "maximum number of pairs allowed to close a motif")
  ;
  boost_po::variables_map vm;
  boost_po::store(boost_po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
  boost_po::notify(vm);
  
  if(vm.count("help"))  {
    print_help_message();
    std::cout << desc << std::endl;
    exit(0);
  }
  
  if(vm.count("structure"))  {
    if(*(-- output_directory.end()) != '/')  {
      output_directory += "/";
    }
    // get the hash of the PDBIDs
    std::unordered_map<int, PDBIDType> hash_PDBID;
    if(vm.count("map_pdb"))  {
      MapPDBIndex index_ID(map_PDBID_file);
      index_ID.GetPDBMapping(hash_PDBID);
    }
    AnnotatedStructure whole_structure(structure_file, hash_PDBID);
    whole_structure.RemoveStackingPairs(max_closing_pairs);
	  whole_structure.BuildShortestPath(extend_size);
	  if(vm.count("query_motif"))  {
	    StructuralMotif query(query_file);
	    whole_structure.PredictLocalMotifs(query, extend_size, max_num_strands);
	  } else  {
	    whole_structure.PredictLocalMotifs(extend_size, max_num_strands);
	  }
	  std::list<MotifBoundType> target_list;
    whole_structure.CopyMotifInstances(target_list);
    for(auto it = target_list.begin(); it != target_list.end(); ++ it) {
      std::list<MotifBoundType> permuted_targets;
      whole_structure.PermuteStrandOrientation(*it, permuted_targets);
      for(auto it_permuted = permuted_targets.begin(); it_permuted != permuted_targets.end(); ++ it_permuted) {
        StructuralMotif instance;
        whole_structure.GenerateMotif(*it_permuted, instance);
        std::string file_info = instance.GetFileInfo();
        std::string file_out = output_directory + file_info + ".rmf";
        instance.WriteMotifToFile(file_out);
      }
    }
  } else  {
    std::cout << "Source RNA structure not specified!" << std::endl;
    print_help_message();
    std::cout << desc << std::endl;
    exit(0);
  }
  
  return 0;
}
