#include "map_PDB_index.h"

using namespace std;

MapPDBIndex::MapPDBIndex(void)  {
  return;
}

MapPDBIndex::~MapPDBIndex(void) {
  return;
}

MapPDBIndex::MapPDBIndex(string index_file) {
  ifstream id_fh(index_file.c_str(), ios_base::in);
  if(!id_fh.good())  {
    cout << "Error in reading " << index_file << endl;
    return;
  }
  char* buff = new char[MAX];
  while(!id_fh.eof()) {
    id_fh.getline(buff, MAX - 1);
    if(id_fh.eof())  {
      break;
    }
    string line(buff);
    vector<string> decom;
    boost::split(decom, line, boost::is_any_of( " \t" ), boost::token_compress_on);
    if(decom.size() != 2)  {
      cout << "Error in reading PDB index file:" << endl;
      cout << "content in file: " << line << endl;
      exit(1);
    }
    int seq_id = atoi(decom[1].c_str());
    index_map_[seq_id] = decom[0];
  }
  delete [] buff;
  return;
}

void MapPDBIndex::GetPDBMapping(std::unordered_map<int, PDBIDType>& nuc_map)  {
  nuc_map = index_map_;
  return;
}
