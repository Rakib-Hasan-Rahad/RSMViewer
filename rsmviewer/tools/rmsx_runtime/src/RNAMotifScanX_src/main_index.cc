#include "map_PDB_index.h"
#include <iostream>
#include <unordered_map>

using namespace std;

int main()  {
  MapPDBIndex hash_PDB_index("./PDB/1S72_0.rmsx.nch");
  unordered_map<int, PDBIDType> nuc_hash;
  hash_PDB_index.GetPDBMapping(nuc_hash);
  
  for(auto it = nuc_hash.begin(); it != nuc_hash.end(); ++ it) {
    cout << it->first << "  " << it->second << endl; 
  }
  return 0;
}
