#include "map_PDB_index.h"

#include <cassert>
#include <iostream>
#include <vector>
#include <string>
#include <list>
#include <cmath>
#ifndef _STRUCTURAL_MOTIF_H_
#define _STRUCTURAL_MOTIF_H_

#ifndef MAX
#define MAX 30000
#endif

class AnnotatedStructure;

class StructuralMotif	{

 public:
	StructuralMotif & operator=(StructuralMotif &input_motif);
	// constructor functions
	StructuralMotif(void);
	StructuralMotif(std::string file_name);	// construct the motif from a file
	StructuralMotif(	// construct the motif by specifying all required information
		const std::string in_file_info, 
		const std::string in_sequence, 
		const std::vector<int>& in_break_points, 
		const std::vector<std::vector<int> >& in_adjacency, 
		const std::vector<std::vector<int> >& in_interaction, 
		const std::vector<int>& in_chain_stacking
	);
	StructuralMotif(	// construct the motif by specifying all required information
		const std::string in_file_info, 
		const std::string in_sequence, 
		const std::vector<int>& in_break_points, 
		const std::vector<std::vector<int> >& in_adjacency, 
		const std::vector<std::vector<int> >& in_interaction, 
		const std::vector<int>& in_chain_stacking,
		const std::unordered_map<int, PDBIDType> in_hash_PDBID
	);
	// destructor function
	~StructuralMotif(void);
	inline std::string GetFileInfo(void)  {
	  return file_info_;
	}
	void PrintMotifInfo(void);
	void WriteMotifToFile(std::string file_name);
	int GetUpDistance(int n); // find the distance of the current nucleotide to the begin of its strand
	int GetDownDistance(int n); // find the distance of the current nucleotide to the end of its strand
	bool IsSameStrand(const int& a, const int& b);
	friend class AnnotatedStructure;
	friend class MotifGraphMatching;
	friend class Simulation;
	friend class ScanStructure;
  
 protected:
	std::string file_info_;
	std::string sequence_;
	std::vector<int> break_points_;
	//	used to store basepairs and non-adjacent stackings
	std::vector<std::vector<int> > interaction_adjacency_;	// two dimensional metrix that stores the pairwise interaction
	std::vector<std::vector<int> > record_interaction_;		// contains detailed information for print-out (affiliate to interaction_adjacency)
	//	used to store adjacent stackings
	std::vector<int> chain_stacking_;
	std::unordered_map<int, PDBIDType> hash_PDBID_;
	
	virtual void parse_basepair(char *content);
	virtual void parse_stacking(char *content);
	
	virtual int hash_basepair(bool is_cis, char edge_1, char edge_2, char nuc_1, char nuc_2);
	virtual char reverse_hash_nucleotide(int value);
	virtual int hash_nucleotide(char nuc);
	virtual char reverse_hash_edge(int value);
	virtual int hash_edge(char edge);
	virtual char reverse_hash_stacking(int value);
	virtual int hash_stacking(char type);
	
	virtual bool is_noncanonical(const std::vector<int>& interaction);
	virtual bool is_noncanonical(const int pair_label);
	virtual bool is_canonical(const std::vector<int>& interaction);
	virtual bool is_canonical(const int pair_label);
	
	virtual std::string interpret_hash_value(int value);
	virtual void interpret_base_pair(int value, char& edge_a, char& edge_b, bool& is_cis);
	
};

#endif
