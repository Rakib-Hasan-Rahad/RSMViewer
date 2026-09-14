#include "structural_motif.h"

#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <ctype.h>


using namespace std;

StructuralMotif::StructuralMotif(void)	{
	return;
}

StructuralMotif::StructuralMotif(string file_name)	{
	//	convert string to char
	int i;
	char *file_name_holder = new char [file_name.length() + 1];
	for(i = 0; i < (int) file_name.length(); ++ i)	{
		file_name_holder[i] = file_name[i];
	}
	file_name_holder[file_name.length()] = '\0';
	
	//	open the file
	ifstream input_motif(file_name_holder, ifstream::in);
	if(!input_motif.good())	{
		cout << "Error in reading motif!" << endl;
	}
	
	char *content_holder = new char [MAX];
	
	//	handle file header
	input_motif.getline(content_holder, MAX);
	file_info_ = string(content_holder);
	file_info_ = file_info_.substr(1);
	//	handle sequence
	input_motif.getline(content_holder, MAX);
	sequence_ = string();
	char prev = '.';
	for(i = 0; content_holder[i] != '\0'; ++ i)	{
		if(content_holder[i] != '.')	{
			sequence_ += content_holder[i];
		}	else	{
			if(prev != '.')	{
				break_points_.push_back(sequence_.length());
			}
		}
		prev = content_holder[i];
	}
	
	//	allocate space for basepair and stacking matrices
	chain_stacking_ = vector<int> (sequence_.length(), MAX);
	interaction_adjacency_ = vector< vector<int> > (sequence_.length());
	for(i = 0; i < (int) sequence_.length(); ++ i)	{
		interaction_adjacency_[i] = vector<int> (sequence_.length(), MAX);
	}
	
	string info_block;
	//	handle base pair
	input_motif.getline(content_holder, MAX);
	info_block = string(content_holder);
	if(info_block.compare("#info=basepair"))	{
		cout << "Corrupted motif file:	" << info_block << endl;
	}
	char ch;
	do	{
		input_motif.getline(content_holder, MAX);
		ch = content_holder[0];
		if(isdigit(ch))	{
			parse_basepair(content_holder);
		}	else	{
			info_block = string(content_holder);
		}
	}	while(isdigit(ch) && !input_motif.eof());
	
	if(input_motif.eof())	{
		cout << "No base-stacking information is found. Performing alignment based on base-pairing information." << endl;
		input_motif.close();
		return;
	}
	
	//	handle stacking
	if(info_block.compare("#info=stacking"))	{
		cout << "Corrupted motif file:	" << info_block << endl; 
	}
	while(!input_motif.eof())	{
		input_motif.getline(content_holder, MAX);
		if(isdigit(content_holder[0]))	{
			parse_stacking(content_holder);
		}
	}
	delete [] file_name_holder;
	delete [] content_holder;
	return;
}

StructuralMotif& StructuralMotif::operator=(StructuralMotif &input_motif)	{

	this->sequence_ = input_motif.sequence_;
	this->file_info_ = input_motif.file_info_;
	this->break_points_ = input_motif.break_points_;
	this->chain_stacking_ = input_motif.chain_stacking_;
	
	for(auto it = input_motif.interaction_adjacency_.begin(); 
	    it != input_motif.interaction_adjacency_.end(); ++ it) {
	  this->interaction_adjacency_.push_back(*it);
	}
	
	for(auto it = input_motif.record_interaction_.begin(); it != input_motif.record_interaction_.end(); ++ it) {
	  this->record_interaction_.push_back(*it);
	}
	
	return *this;
}

StructuralMotif::StructuralMotif(
		const string in_file_info, 
		const string in_sequence, 
		const vector<int>& in_break_points, 
		const vector<vector<int> >& in_adjacency, 
		const vector<vector<int> >& in_interaction, 
		const vector<int>& in_chain_stacking
)	{
	file_info_ = in_file_info;
	sequence_ = in_sequence;
	break_points_ = in_break_points;
	chain_stacking_ = in_chain_stacking;
	
	unsigned int i;
	interaction_adjacency_.resize(in_adjacency.size());
	for(i = 0; i < in_adjacency.size(); ++ i)	{
		interaction_adjacency_[i] = in_adjacency[i];
	}	
	for(auto it = in_interaction.begin(); it != in_interaction.end(); ++ it)	{
		record_interaction_.push_back(*it);
	}		
	return;
}

StructuralMotif::StructuralMotif(	// construct the motif by specifying all required information
		const std::string in_file_info, 
		const std::string in_sequence, 
		const std::vector<int>& in_break_points, 
		const std::vector<std::vector<int> >& in_adjacency, 
		const std::vector<std::vector<int> >& in_interaction, 
		const std::vector<int>& in_chain_stacking,
		const std::unordered_map<int, PDBIDType> in_hash_PDBID
) {
  file_info_ = in_file_info;
	sequence_ = in_sequence;
	break_points_ = in_break_points;
	chain_stacking_ = in_chain_stacking;
	
	unsigned int i;
	interaction_adjacency_.resize(in_adjacency.size());
	for(i = 0; i < in_adjacency.size(); ++ i)	{
		interaction_adjacency_[i] = in_adjacency[i];
	}	
	for(auto it = in_interaction.begin(); it != in_interaction.end(); ++ it)	{
		record_interaction_.push_back(*it);
	}		
	hash_PDBID_ = in_hash_PDBID;
	return;
}

StructuralMotif::~StructuralMotif(void)	{
	return;
}

int StructuralMotif::GetUpDistance(int n) {
  int b = 0;
  for(auto it = break_points_.begin(); it != break_points_.end(); ++ it) {
    if(*it >= b && *it <= n)  {
      b = *it;
    }
  }
  return (n - b);
} 

int StructuralMotif::GetDownDistance(int n) {
  int e = sequence_.length() - 1;
  for(auto it = break_points_.begin(); it != break_points_.end(); ++ it) {
    if(*it >= n && *it <= e)  {
      e = *it;
    }
  }
  return (e - n);
} 

void StructuralMotif::WriteMotifToFile(string file_name) {
  ofstream motif_fout(file_name.c_str(), ios_base::out);
  if(!motif_fout.good())  {
    cout << "Error in writing to file:  " << file_name << endl;
    return;
  }
  unsigned int i;
  motif_fout << ">" << file_info_ << endl;
  int b_id = 0;
	for(i = 0; i < sequence_.length(); ++ i) {
	  motif_fout << sequence_[i];
	  if(b_id < (int) break_points_.size() && (int) i + 1 == break_points_[b_id])  {
	    motif_fout << "...";
	    ++ b_id;
	  }
	}
  motif_fout << endl;
  motif_fout << "#info=basepair" << endl;
	for(auto it = record_interaction_.begin(); it != record_interaction_.end(); ++ it) {
	  if((*it)[2] >= 0 && (*it)[2] != MAX)  {
	    motif_fout << (*it)[0] << "-" << (*it)[1] << ",";
	    string info = interpret_hash_value((*it)[2]);
	    motif_fout << info;
	    auto it_begin_ID = hash_PDBID_.find((*it)[0]);
	    auto it_end_ID = hash_PDBID_.find((*it)[1]);
      if(it_begin_ID != hash_PDBID_.end() && it_end_ID != hash_PDBID_.end())  {
        motif_fout << "," << it_begin_ID->second << "-" << it_end_ID->second;
      }
      motif_fout << endl;
	  }
	}
	motif_fout << "#info=stacking" << endl;
	for(auto it = record_interaction_.begin(); it != record_interaction_.end(); ++ it) {
	  if((*it)[2] < 0)  {
	    motif_fout << (*it)[0] << "-" << (*it)[1] << ",";
	    string info = interpret_hash_value((*it)[2]);
	    motif_fout << info;
	    auto it_begin_ID = hash_PDBID_.find((*it)[0]);
	    auto it_end_ID = hash_PDBID_.find((*it)[1]);
      if(it_begin_ID != hash_PDBID_.end() && it_end_ID != hash_PDBID_.end())  {
        motif_fout << "," << it_begin_ID->second << "-" << it_end_ID->second;
      }
      motif_fout << endl;
	  }
	}
	for(i = 0; i < chain_stacking_.size(); ++ i) {
	  if(chain_stacking_[i] != MAX)  {
	    motif_fout << i << "-" << (i + 1) << "," << interpret_hash_value(chain_stacking_[i]);
	    auto it_begin_ID = hash_PDBID_.find(i);
	    auto it_end_ID = hash_PDBID_.find(i + 1);
      if(it_begin_ID != hash_PDBID_.end() && it_end_ID != hash_PDBID_.end())  {
        motif_fout << "," << it_begin_ID->second << "-" << it_end_ID->second;
      }
      motif_fout << endl;
	  }
	}
  return;
}

void StructuralMotif::PrintMotifInfo(void)	{
  unsigned int i;
  cout << ">" << file_info_ << endl;
  int b_id = 0;
	for(i = 0; i < sequence_.length(); ++ i) {
	  cout << sequence_[i];
	  if(b_id < (int) break_points_.size() && (int) i + 1 == break_points_[b_id])  {
	    cout << "...";
	    ++ b_id;
	  }
	}
  cout << endl;
  cout << "#info=basepair" << endl;
	for(auto it = record_interaction_.begin(); it != record_interaction_.end(); ++ it) {
	  if((*it)[2] >= 0 && (*it)[2] != MAX)  {
	    cout << (*it)[0] << "-" << (*it)[1] << ",";
	    string info = interpret_hash_value((*it)[2]);
	    cout << info << endl;
	  }
	}
	cout << "#info=stacking" << endl;
	for(auto it = record_interaction_.begin(); it != record_interaction_.end(); ++ it) {
	  if((*it)[2] < 0)  {
	    cout << (*it)[0] << "-" << (*it)[1] << ",";
	    string info = interpret_hash_value((*it)[2]);
	    cout << info << endl;
	  }
	}
	for(i = 0; i < chain_stacking_.size(); ++ i) {
	  if(chain_stacking_[i] != MAX)  {
	    cout << i << "-" << (i + 1) << "," << interpret_hash_value(chain_stacking_[i]) << endl;
	  }
	}
	return;
}

void StructuralMotif::parse_basepair(char *content)	{

	int loc_a, loc_b;
	char edge_a, edge_b;
	bool is_cis = true;
	
	int i;
	char *holder = new char [MAX];
	int idx = 0;
	for(i = 0; content[i] != '-'; ++ i)	{
		holder[idx ++] = content[i];
	}
	holder[idx] = '\0';
	loc_a = atoi(holder);
	
	idx = 0;
	for(i = i + 1; content[i] != ','; ++ i)	{
		holder[idx ++] = content[i];
	}
	holder[idx] = '\0';
	loc_b = atoi(holder);
	
	edge_a = content[i + 1];
	edge_b = content[i + 3];
	
	if(content[i + 5] == 'c')	{
		is_cis = true;
	}	else if(content[i + 5] == 't')	{
		is_cis = false;
	}
	
	delete [] holder;
			
	if(edge_a != '-' && edge_b != '-')	{
		//	in this case the basepair is a regular base pair
	
		//	record interaction in two-dimensional matrix
		interaction_adjacency_[loc_a][loc_b] = hash_basepair(is_cis, edge_a, edge_b, sequence_[loc_a], sequence_[loc_b]);
		interaction_adjacency_[loc_b][loc_a] = hash_basepair(is_cis, edge_b, edge_a, sequence_[loc_b], sequence_[loc_a]);
		//	record specific edges and orientations for print-out
		vector<int> single_interaction(6);
		single_interaction[0] = loc_a;
		single_interaction[1] = loc_b;
		single_interaction[2] = interaction_adjacency_[loc_a][loc_b];
		single_interaction[3] = (int) edge_a;
		single_interaction[4] = (int) edge_b;
		single_interaction[5] = (int) is_cis;
		record_interaction_.push_back(single_interaction);
	}	else	{
		//	in this case the basepair is a hydrogen bond
		interaction_adjacency_[loc_a][loc_b] = 0;
		interaction_adjacency_[loc_b][loc_a] = 0;
		//	record specific edges and orientations for print-out
		vector<int> single_interaction(6);
		single_interaction[0] = loc_a;
		single_interaction[1] = loc_b;
		single_interaction[2] = interaction_adjacency_[loc_a][loc_b];
		single_interaction[3] = (int) edge_a;
		single_interaction[4] = (int) edge_b;
		single_interaction[5] = -1;
		record_interaction_.push_back(single_interaction);
	}
	return;
}

bool StructuralMotif::is_canonical(const int pair_label) {
  int c = pair_label;
	// 3: AU
	// 12: UA
	// 6: CG
	// 9: GC
	// 11: GU
	// 14: UG
	if(c == 4 || c == 13 || c == 7 || c == 10 || c == 12 || c == 15)	{
		return true;
	}
	return false;
}

bool StructuralMotif::is_canonical(const vector<int>& interaction)	{
	int c = interaction[2];
	// 3: AU
	// 12: UA
	// 6: CG
	// 9: GC
	// 11: GU
	// 14: UG
	if(c == 4 || c == 13 || c == 7 || c == 10 || c == 12 || c == 15)	{
		return true;
	}
	return false;
}


bool StructuralMotif::is_noncanonical(const int pair_label) {
  int c = pair_label;
	// 4: AU
	// 13: UA
	// 7: CG
	// 10: GC
	// 12: GU
	// 15: UG
	if(c > 0 && c != MAX && c != 4 && c != 13 && c != 7 && c != 10 && c != 12 && c != 15)	{
		return true;
	}
	return false;
}

bool StructuralMotif::is_noncanonical(const vector<int>& interaction)	{
	int c = interaction[2];
	// 4: AU
	// 13: UA
	// 7: CG
	// 10: GC
	// 12: GU
	// 15: UG
	if(c > 0 && c != MAX && c != 4 && c != 13 && c != 7 && c != 10 && c != 12 && c != 15)	{
		return true;
	}
	return false;
}

char StructuralMotif::reverse_hash_nucleotide(int value) {
  if(value == 0)  {
    return 'A';
  } else if(value == 1) {
    return 'C';
  } else if(value == 2) {
    return 'G';
  } else if(value == 3) {
    return 'U';
  } else  {
    return 'N';
  }
}


int StructuralMotif::hash_nucleotide(char nuc)	{
	nuc = (char) toupper((int) nuc);
	if(nuc == 'A')	{
		return 0;
	}	else if(nuc == 'C')	{
		return 1;
	}	else if(nuc == 'G')	{
		return 2;
	}	else if(nuc == 'U' || nuc == 'T')	{
		return 3;
	}	else	{
		return MAX;
	}	
}

char StructuralMotif::reverse_hash_edge(int value) {
  if(value == 0)  {
    return 'W';
  } else if(value == 1) {
    return 'H';
  } else if(value == 2) {
    return 'S';
  } else  {
    return ' ';
  }
}

int StructuralMotif::hash_edge(char edge)	{
	edge = (char) toupper((int) edge);
	if(edge == 'W')	{
		return 0;
	}	else if(edge == 'H')	{
		return 1;
	}	else if(edge == 'S')	{
		return 2;
	}	else	{
		return MAX;
	}
}

int StructuralMotif::hash_basepair(bool is_cis, char edge_1, char edge_2, char nuc_1, char nuc_2)	{
	
	int edge_index_1 = hash_edge(edge_1);
	int edge_index_2 = hash_edge(edge_2);
	int nuc_index_1 = hash_nucleotide(nuc_1);
	int nuc_index_2 = hash_nucleotide(nuc_2);
	
	assert(edge_index_1 >= 0);
	assert(edge_index_2 >= 0);
	assert(nuc_index_1 >= 0);
	assert(nuc_index_2 >= 0);
	
	int hash_index = (edge_index_1 * 3 + edge_index_2) * 2;
	if(!is_cis)	{
		hash_index ++;
	}
	hash_index = hash_index * 16 + (nuc_index_1 * 4 + nuc_index_2);
	return hash_index + 1;
}

void StructuralMotif::parse_stacking(char *content)	{
	
	int loc_a, loc_b;
	
	int i;
	char *holder = new char [MAX];
	int idx = 0;
	for(i = 0; content[i] != '-'; ++ i)	{
		holder[idx ++] = content[i];
	}
	holder[idx] = '\0';
	loc_a = atoi(holder);
	
	idx = 0;
	for(i = i + 1; content[i] != ','; ++ i)	{
		holder[idx ++] = content[i];
	}
	holder[idx] = '\0';
	loc_b = atoi(holder);
	char type = content[i + 1];
	char rev_type;
	if(type == 'u')  {
	  rev_type = 'd';
	} else if(type == 'd') {
	  rev_type = 'u';
	} else  {
	  rev_type = type;
	}
	//cout << "hash basepair:	" << type << "	" << loc_a << "	" << loc_b << endl;
	bool is_true_chain = true;
	for(int j = 0; j < (int) break_points_.size(); ++ j) {
	  if(break_points_[j] == (loc_b > loc_a ? loc_b : loc_a))  {
	    is_true_chain = false;
	    break;
	  }
	}
	if(abs(loc_b - loc_a) > 1 || !is_true_chain)	{
		interaction_adjacency_[loc_a][loc_b] = hash_stacking(type);	
		interaction_adjacency_[loc_b][loc_a] = hash_stacking(rev_type);
		vector<int> single_interaction(6);
		single_interaction[0] = loc_a;
		single_interaction[1] = loc_b;
		single_interaction[2] = interaction_adjacency_[loc_a][loc_b];
		single_interaction[3] = (int) type;
		record_interaction_.push_back(single_interaction);
	}	else if(abs(loc_b - loc_a) == 1)	{
	  if(loc_b > loc_a)  {
	    //cout << "adjacent stacking: " << loc_a << " " << loc_b << endl;
		  chain_stacking_[loc_a] = hash_stacking(type);
		} else  {
		  chain_stacking_[loc_b] = hash_stacking(rev_type);
		}
	}
	delete [] holder;
	return;
}

char StructuralMotif::reverse_hash_stacking(int value)  {
  if(value == -1)  {
    return 'u';
  } else if(value == -2) {
    return 'd';
  } else if(value == -3) {
    return 'i';
  } else if(value == -4) {
    return 'o';
  } else  {
    return ' ';
  }
}

int StructuralMotif::hash_stacking(char type)	{
	if(type == 'u')	{
		//	upward
		return -1;
	}	else if(type == 'd')	{
		//	downward
		return -2;
	}	else if(type == 'i')	{
		//	inward
		return -3;
	}	else if(type == 'o')	{
		//	outward
		return -4;
	}	else	{
		return MAX;
	}
}

void StructuralMotif::interpret_base_pair(int value, char& edge_a, char& edge_b, bool& is_cis) {
  assert(value != MAX && value >= 0);
  if(value == 0) {
    // hydrogen boud
    edge_a = '-';
    edge_b = '-';
    is_cis = true;
    return;
  } else  {
    // base-pair
    int edge_and_orientation = floor((value - 1) / 16);
    is_cis = true;
    if(edge_and_orientation % 2 == 1)  {
      is_cis = false;
    }
    int edges = floor(edge_and_orientation / 2);
    int edge_1 = floor(edges / 3);
    int edge_2 = edges % 3;
    edge_a = reverse_hash_edge(edge_1);
    edge_b = reverse_hash_edge(edge_2);
    return;
  }
}

string StructuralMotif::interpret_hash_value(int value)  {
  assert(value != MAX);
  if(value < 0)  {
    // base-stacking
    if(value == -1)  {
      return "upward";
    } else if(value == -2) {
      return "downward";
    } else if(value == -3) {
      return "inward";
    } else if(value == -4) {
      return "outward";
    }
  } else if(value == 0) {
    // hydrogen boud
    return "-/-,hbond";
  } else  {
    // base-pair
    int edge_and_orientation = floor((value - 1) / 16);
    bool is_cis = true;
    if(edge_and_orientation % 2 == 1)  {
      is_cis = false;
    }
    int edges = floor(edge_and_orientation / 2);
    int edge_1 = floor(edges / 3);
    int edge_2 = edges % 3;
    char* edge_holder = new char[5];
    edge_holder[0] = reverse_hash_edge(edge_1);
    edge_holder[1] = '/';
    edge_holder[2] = reverse_hash_edge(edge_2);
    edge_holder[3] = ',';
    edge_holder[4] = '\0';
    string interaction(edge_holder);
    delete []edge_holder;
    if(is_cis)  {
      interaction += "cis";
    } else  {
      interaction += "trans";
    }
    return interaction;
  }
  return "N/A";
}

bool StructuralMotif::IsSameStrand(const int& a, const int& b)  {
  for(auto it = break_points_.begin(); it != break_points_.end(); ++ it) {
    if((a < *it && b >= *it) || (b < *it && a >= *it))  {
      return false;
    }
  }
  return true;
}
