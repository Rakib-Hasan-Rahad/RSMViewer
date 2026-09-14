#include "structural_motif.h"
#include "annotated_structure.h"

#include <iostream>
#include <cstdlib>

using namespace std;

int main()	{

	AnnotatedStructure whole_structure("/home/cczhong/Codes/RNAMotifScanX/PDB/1S72_0.rmsx.in");
	//cout << "END OF CONSTRUCTING WHOLE STRUCTURE" << endl;
	StructuralMotif motif("/home/cczhong/Codes/RNAMotifScanX/Queries/k-turn_consensus.struct");
	
	/*
	list<vector<int> > orders;
	whole_structure.all_permutation(4, orders);
	for(auto it = orders.begin(); it != orders.end(); ++ it) {
	  for(auto it_vector = it->begin(); it_vector != it->end(); ++ it_vector) {
	    cout << *it_vector << "\t";
	  }
	  cout << endl;
	}
	
	return 0;
	*/
	
	whole_structure.RemoveStackingPairs(1);
	whole_structure.BuildShortestPath(10);
	
	list<MotifBoundType> target_list;
	whole_structure.PredictLocalMotifs(motif, 5, 3);
  whole_structure.CopyMotifInstances(target_list);
  for(auto it = target_list.begin(); it != target_list.end(); ++ it) {
    StructuralMotif instance;
    whole_structure.GenerateMotif(*it, instance);
    string file_info = instance.GetFileInfo();
    string file_out = "/home/cczhong/Codes/RNAMotifScanX/WorkSpace/" + file_info + ".rmf";
    instance.WriteMotifToFile(file_out);
  }
	
	return 0;
}
