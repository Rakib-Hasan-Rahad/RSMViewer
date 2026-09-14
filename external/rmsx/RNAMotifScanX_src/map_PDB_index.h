#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

#include <unordered_map>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <cstdlib>

#ifndef _MAP_PDB_INDEX_H_
#define _MAP_PDB_INDEX_H_

#ifndef MAX
#define MAX 30000
#endif

typedef std::string PDBIDType;

class MapPDBIndex {
 public:
  MapPDBIndex();
  MapPDBIndex(std::string index_file);
  ~MapPDBIndex();
  void GetPDBMapping(std::unordered_map<int, PDBIDType>& nuc_map);
 protected:
  std::unordered_map<int, PDBIDType> index_map_;
};

#endif
