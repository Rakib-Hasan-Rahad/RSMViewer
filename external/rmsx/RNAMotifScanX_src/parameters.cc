#include "parameters.h"

using namespace std;

Parameters::Parameters(void)	{
	return;
}

Parameters::Parameters(ParameterType input)	{
	
	gap_opening = int (100 * input.gap_opening);
	gap_extension = int (100 * input.gap_extension);
	triple_interaction_bonus = int (100 * input.triple_interaction_bonus);
	
	hbond_match_base = int (100 * input.hbond_match_base);
	hbond_match_bonus = int (100 * input.hbond_match_bonus);
	
	weight_isosteric = int (100 * input.weight_isosteric);
	weight_nonadjacent_stacking = int (100 * input.weight_nonadjacent_stacking);
	weight_adjacent_stacking = int (100 * input.weight_adjacent_stacking);
	weight_sequence = int (100 * input.weight_sequence);
	
	pair_deletion = int (100 * input.pair_deletion); 
	stacking_deletion = int (100 * input.stacking_deletion);
	
	asym_nuc = int (100 * input.asym_nuc); 
	asym_loop = int (100 * input.asym_loop); 
	cons_stacking = int (100 * input.cons_stacking); 
	stack_to_bulge = int (100 * input.stack_to_bulge); 
	stack_to_internal_asym = int (100 * input.stack_to_internal_asym); 
	stack_to_internal_sym = int (100 * input.stack_to_internal_sym);
	
	initialize_basepair_index();
	read_isosteric_scoring_function(input.isosteric_file_name);
	read_stacking_scoring_function(input.stacking_file_name);
	read_nucleotide_scoring_function(input.nucleotide_file_name);
	compute_edge_substitution();
	
	return;
}

Parameters::~Parameters(void)	{
	return;
}	


Parameters & Parameters::operator=(const Parameters &input_parameters)	{
	
	//	copy all scores
	this->gap_opening = input_parameters.gap_opening;
	this->gap_extension = input_parameters.gap_extension;
	this->triple_interaction_bonus = input_parameters.triple_interaction_bonus;
	
	this->hbond_match_base = input_parameters.hbond_match_base;
	this->hbond_match_bonus = input_parameters.hbond_match_bonus;
	
	this->weight_isosteric = input_parameters.weight_isosteric;
	this->weight_nonadjacent_stacking = input_parameters.weight_nonadjacent_stacking;
	this->weight_adjacent_stacking = input_parameters.weight_adjacent_stacking;
	this->weight_sequence = input_parameters.weight_sequence;
	
	this->pair_deletion = input_parameters.pair_deletion; 
	this->stacking_deletion = input_parameters.stacking_deletion;
	
	this->asym_nuc = input_parameters.asym_nuc; 
	this->asym_loop = input_parameters.asym_loop; 
	this->cons_stacking = input_parameters.cons_stacking; 
	this->stack_to_bulge = input_parameters.stack_to_bulge; 
	this->stack_to_internal_asym = input_parameters.stack_to_internal_asym; 
	this->stack_to_internal_sym = input_parameters.stack_to_internal_sym;
	
	this->isosteric_max_score = input_parameters.isosteric_max_score;
	this->stacking_max_score = input_parameters.stacking_max_score;
	this->edge_max_score = input_parameters.edge_max_score;
	
	//	copy all matrices
	int i;
	this->isosteric_substitution = vector< vector<int> > (64);
	for(i = 0; i < 64; ++ i)	{
		this->isosteric_substitution[i] = input_parameters.isosteric_substitution[i];
	}
	this->stacking_substitution = vector< vector<int> > (4);
	for(i = 0; i < 4; ++ i)	{
		this->stacking_substitution[i] = input_parameters.stacking_substitution[i];
	}
	this->nucleotide_substitution = vector< vector<int> > (4);
	for(i = 0; i < 4; ++ i)	{
		this->nucleotide_substitution[i] = input_parameters.nucleotide_substitution[i];
	}
	this->edge_substitution = vector< vector<int> > (12);
	for(i = 0; i < 12; ++ i)	{
		this->edge_substitution[i] = input_parameters.edge_substitution[i];
	}
	
	// copy the basepair_index
	this->basepair_index = input_parameters.basepair_index;
	return *this;
}	

int Parameters::match_isosteric_basepair(int i, int j)	{
	//	return the corresponding score of matching the two basepairs
	//	note that i and j are 1-based while the basepair_index matrix is 0-based
	//	check for whether theoretical impossible basepairs are present
	if(i >= 0 && j >= 0 && basepair_index[i] > 0 && basepair_index[j] > 0)	{
		//	if both base pairs are regular base pairs that can present
		//cout << "within match_isosteric_basepair:	" << basepair_index[i] << "	" << basepair_index[j] << endl;
		int i1 = (i - 1) / 16, j1 = (j - 1) / 16;
		int iori = i1 % 2, jori = j1 % 2;
		int i2 = i1 / 2, j2 = j1 / 2;
		int iedga = i2 % 3, iedgb = i2 / 3, jedga = j2 % 3, jedgb = j2 / 3;
		
		int flip_bonus = 0;
		if(iori == jori && iedga != jedga && iedgb != jedgb && iedga == jedgb && iedgb == jedga)	{
			//cout << "entered!" << endl;
			flip_bonus = (int) (isosteric_substitution[basepair_index[i] - 1][basepair_index[j] - 1] * 0.5);
		}
		return isosteric_substitution[basepair_index[i] - 1][basepair_index[j] - 1] + flip_bonus;
	}	else if(i >= 0 && j >= 0 && (basepair_index[i] == 0 || basepair_index[j] == 0))	{
		//	if either base pair is a hydrogen bond
		return hbond_match_base;
	}	else if(i >= 0 && j >= 0 && (basepair_index[i] == 0 && basepair_index[j] == 0))	{
		//	if both base pair are hydrogen bonds
		return hbond_match_base + hbond_match_bonus;
	}
	return 0;
}

int Parameters::max_pair_score(int i)	{
	//	return the maximal score for matching this basepair
	//	note that i and j are 1-based while the basepair_index matrix is 0-based
	
	if(basepair_index[i] - 1 >= 0)	{
		//	if the base pair is a regular base pair
		return isosteric_max_score[basepair_index[i] - 1];
	}	else if(basepair_index[i] == -1)	{
		//	if the base pair cannot exist theoretically
		int edge_index_i = i / 16;
		return edge_max_score[edge_index_i];
	}	else if(basepair_index[i] == 0)	{
		//	if the base pair is a hydrogen bond
		return hbond_match_base + hbond_match_bonus;
	}
	return 0;
}

int Parameters::match_stacking_basepair(int i, int j)	{
	//	return the corresponding score of matching two basestackings
	if(i < 0 && j < 0)	{
		i = abs(i) - 1;
		j = abs(j) - 1;
		return stacking_substitution[i][j];
	}	else if((i == 0 && j < 0)|| (i < 0 && j == 0))	{
		//	note that we do not consider an hydrogen bond match with a base stacking
		return 0;
	}	else if(i == 0 && j == 0)	{
		return hbond_match_base + hbond_match_bonus;
	}
	return 0;
}

int Parameters::max_stacking_score(int i)	{
	//	return the maximal score of matching this base-stacking
	if(i < 0)	{
		i = abs(i);
		return stacking_max_score[i];
	}	else if(i == 0)	{
		return 0;
	}
	return 0;
}

int Parameters::match_nucleotide(char a, char b)	{
	//	return the corresponding score of matching two nucleotides
	return nucleotide_substitution[hash_nucleotide(a)][hash_nucleotide(b)];
}

void Parameters::read_isosteric_scoring_function(string file_name)	{
	int i, j, k;
	//	allocate space
	isosteric_substitution = vector< vector<int> > (64);
	for(i = 0; i < 64; ++ i)	{
		isosteric_substitution[i] = vector<int> (64);
	}
	//	convert string to char*
	char *file_name_holder = new char [file_name.length() + 1];
	for(i = 0; i < (int) file_name.length(); ++ i)	{
		file_name_holder[i] = file_name[i];
	}
	file_name_holder[file_name.length()] = '\0';
	ifstream input_scoring_function(file_name_holder, ifstream::in);
	if(!input_scoring_function.good())	{
		cout << "Error in reading isosteric scoring function file!" << endl;
		cout << "You can try to specify its location using the -s option." << endl;
		exit(0);
	}
	//	read the file line by line
	int row_index = 0;
	while(!input_scoring_function.eof())	{
		char *a_line = new char [MAX];
		input_scoring_function.getline(a_line, MAX);
		if(a_line[0] == '#' || a_line[0] == '\0')	{
			//	if the line is a comment or contains nothing
			continue;
		}	else	{
			//	if the line contains the scoring function
			string parse_line(a_line);
			if(parse_line[parse_line.length() - 1] != '\t')	{
				parse_line.append(1, '\t');
			}
			string number_holder("");
			int column_index = 0;
			for(k = 0; k < (int) parse_line.length(); ++ k)	{
				if(parse_line[k] != '\t')	{
					number_holder.append(1, parse_line[k]);
				}	else	{
					//	transform string to char*
					char *from_string = new char [MAX];
					for(j = 0; j < (int) number_holder.length(); ++ j)	{
						from_string[j] = number_holder[j];
					}
					from_string[number_holder.length()] = '\0';
					//	record the number
					isosteric_substitution[row_index][column_index] = (int) (atof(from_string) * 100);
					++ column_index;
					//	clear existing information
					delete [] from_string;
					number_holder = string("");
				}
			}
			//	double check the numbers in each row is 41
			if(column_index != 64)	{
				cout << "Error in column-wise isosteric basepair substitution file!" << endl;
			}
			++ row_index;
		}
		delete [] a_line;
	}
	//	double check the number of rows is 41
	if(row_index != 64)	{
		cout << "Error in row-wise isosteric basepair substitution file!" << endl;
	}
	input_scoring_function.close();
	delete [] file_name_holder;
	
	//	compute the maximum score for each column
	isosteric_max_score = vector<int> (64);
	int total_score = 0;
	for(i = 0; i < 64; ++ i)	{
		int max_score = -MAX;
		for(j = 0; j < 64; ++ j)	{
			if(isosteric_substitution[i][j] > max_score)	{
				max_score = isosteric_substitution[i][j];
			}
			total_score += isosteric_substitution[i][j];
		}
		isosteric_max_score[i] = max_score;
		//cout << isosteric_max_score[i] << endl;
	}
	
	return;
}

void Parameters::read_stacking_scoring_function(string file_name)	{
	int i, j, k;
	//	allocate space
	stacking_substitution = vector< vector<int> > (4);
	for(i = 0; i < 4; ++ i)	{
		stacking_substitution[i] = vector<int> (4);
	}
	//	convert string to char*
	char *file_name_holder = new char [file_name.length() + 1];
	for(i = 0; i < (int) file_name.length(); ++ i)	{
		file_name_holder[i] = file_name[i];
	}
	file_name_holder[file_name.length()] = '\0';
	ifstream input_scoring_function(file_name_holder, ifstream::in);
	if(!input_scoring_function.good())	{
		cout << "Error in reading isosteric scoring function file!" << endl;
		cout << "You can try to specify its location using the -t option." << endl;
		exit(0);
	}
	//	read the file line by line
	int row_index = 0;
	while(!input_scoring_function.eof())	{
		char *a_line = new char [MAX];
		input_scoring_function.getline(a_line, MAX);
		if(a_line[0] == '#' || a_line[0] == '\0')	{
			//	if the line is a comment or contains nothing
			continue;
		}	else	{
			//	if the line contains the scoring function
			string parse_line(a_line);
			if(parse_line[parse_line.length() - 1] != '\t')	{
				parse_line.append(1, '\t');
			}
			string number_holder("");
			int column_index = 0;
			for(k = 0; k < (int) parse_line.length(); ++ k)	{
				if(parse_line[k] != '\t')	{
					number_holder.append(1, parse_line[k]);
				}	else	{
					//	transform string to char*
					char *from_string = new char [MAX];
					for(j = 0; j < (int) number_holder.length(); ++ j)	{
						from_string[j] = number_holder[j];
					}
					from_string[number_holder.length()] = '\0';
					//	record the number
					stacking_substitution[row_index][column_index] = (int) (atof(from_string) * 100);
					++ column_index;
					//	clear existing information
					delete [] from_string;
					number_holder = string("");
				}
			}
			//	double check the numbers in each row is 5
			if(column_index != 4)	{
				cout << "Error in column-wise base stacking substitution file!" << endl;
			}
			++ row_index;
		}
		delete [] a_line;
	}
	//	double check the number of rows is 4
	if(row_index != 4)	{
		cout << "Error in row-wise base stacking substitution file!" << endl;
	}
	input_scoring_function.close();
	delete [] file_name_holder;
	
	//	compute maximum score for each column
	int total_score = 0;
	stacking_max_score = vector<int> (4);
	for(i = 0; i < 4; ++ i)	{
		int max_score = -MAX;
		for(j = 0; j < 4; ++ j)	{
			if(stacking_substitution[i][j] > max_score)	{
				max_score = stacking_substitution[i][j];
			}
			total_score += stacking_substitution[i][j];
		}
		stacking_max_score[i] = max_score;
		//cout << stacking_max_score[i] << endl;
	}
	
	return;
}


void Parameters::read_nucleotide_scoring_function(string file_name)	{
	int i, j, k;
	//	allocate space
	nucleotide_substitution = vector< vector<int> > (4);
	for(i = 0; i < 4; ++ i)	{
		nucleotide_substitution[i] = vector<int> (4);
	}
	//	convert string to char*
	char *file_name_holder = new char [file_name.length() + 1];
	for(i = 0; i < (int) file_name.length(); ++ i)	{
		file_name_holder[i] = file_name[i];
	}
	file_name_holder[file_name.length()] = '\0';
	ifstream input_scoring_function(file_name_holder, ifstream::in);
	if(!input_scoring_function.good())	{
		cout << "Error in reading isosteric scoring function file!" << endl;
		cout << "You can try to specify its location using the -n option." << endl;
		exit(0);
	}
	//	read the file line by line
	int row_index = 0;
	while(!input_scoring_function.eof())	{
		char *a_line = new char [MAX];
		input_scoring_function.getline(a_line, MAX);
		if(a_line[0] == '#' || a_line[0] == '\0')	{
			//	if the line is a comment or contains nothing
			continue;
		}	else	{
			//	if the line contains the scoring function
			string parse_line(a_line);
			if(parse_line[parse_line.length() - 1] != '\t')	{
				parse_line.append(1, '\t');
			}
			string number_holder("");
			int column_index = 0;
			for(k = 0; k < (int) parse_line.length(); ++ k)	{
				if(parse_line[k] != '\t')	{
					number_holder.append(1, parse_line[k]);
				}	else	{
					//	transform string to char*
					char *from_string = new char [MAX];
					for(j = 0; j < (int) number_holder.length(); ++ j)	{
						from_string[j] = number_holder[j];
					}
					from_string[number_holder.length()] = '\0';
					//	record the number
					nucleotide_substitution[row_index][column_index] = (int) (atof(from_string) * 100);
					++ column_index;
					//	clear existing information
					delete [] from_string;
					number_holder = string("");
				}
			}
			//	double check the numbers in each row is 41
			if(column_index != 4)	{
				cout << "Error in column-wise isosteric basepair substitution file!" << endl;
			}
			++ row_index;
		}
		delete [] a_line;
	}
	//	double check the number of rows is 41
	if(row_index != 4)	{
		cout << "Error in row-wise isosteric basepair substitution file!" << endl;
	}
	input_scoring_function.close();
	delete [] file_name_holder;
	return;
}

void Parameters::compute_edge_substitution(void)	{
	//	set up the number of isosteric groups in each edge-combination
	vector<int> edge_size(12);
	edge_size[0] = 6; edge_size[1] = 6; edge_size[2] = 4; edge_size[3] = 5;
	edge_size[4] = 5; edge_size[5] = 4; edge_size[6] = 2; edge_size[7] = 3;
	edge_size[8] = 1; edge_size[9] = 2; edge_size[10] = 1; edge_size[11] = 2;
	//	computes the edge substitution score by averaging the basepair substitution score
	int i, j, k, l;
	edge_substitution = vector< vector<int> > (12);
	for(i = 0; i < 12; ++ i)	{
		edge_substitution[i] = vector<int> (12, 0);
	}
	
	int i_start = 0;
	for(i = 0; i < 12; ++ i)	{
		int j_start = 0;
		for(j = 0; j < 12; ++ j)	{
			int sum_score = 0;
			int total_size = 0;
			for(k = 0; k < edge_size[i]; ++ k)	{
				for(l = 0; l < edge_size[j]; ++ l)	{
					sum_score += isosteric_substitution[i_start + k][j_start + l];
					++ total_size;
				}
			}
			edge_substitution[i][j] = sum_score / total_size;	
			j_start += edge_size[j];
		}
		i_start += edge_size[i];
	}
	
	//	compute maximum matching score for each edge combination
	edge_max_score = vector<int> (12, 0);
	for(i = 0; i < 12; ++ i)	{
		int max_score = -MAX;
		for(j = 0; j < 12; ++ j)	{
			if(edge_substitution[i][j] > max_score)	{
				max_score = edge_substitution[i][j];
			}
		}
		edge_max_score[i] = max_score;
	}
	return;
}

int Parameters::hash_nucleotide(char nuc)	{
	nuc = (char) toupper((int) nuc);
	if(nuc == 'A')	{
		return 0;
	}	else if(nuc == 'C')	{
		return 1;
	}	else if(nuc == 'G')	{
		return 2;
	}	else if(nuc == 'U')	{
		return 3;
	}	else	{
		return -1;
	}	
}

int Parameters::hash_edge(char edge)	{
	edge = (char) toupper((int) edge);
	if(edge == 'W')	{
		return 0;
	}	else if(edge == 'H')	{
		return 1;
	}	else if(edge == 'S')	{
		return 2;
	}	else	{
		return -1;
	}
}



void Parameters::initialize_basepair_index(void)	{

	basepair_index = vector<int> (289);
	
	//	<*****matrix_0*****>
	basepair_index[0] = 0;
	
	//	the theoretically impossible base pairs are treated as hydrogen bonds
	//	Note that if basepair_index[i] > 0, then it means it is a possible base pair
	//	the score should be located at (basepair_index[i] - 1) of the isosteric base pair 
	//	substitution matrix
	
	
	//	<*****matrix_1*****>
	basepair_index[1] = 4;		//	cis	W/W	AA	(0)
	basepair_index[2] = 2;		//	cis	W/W	AC	(1)
	basepair_index[3] = 3;		//	cis	W/W	AG	(2)
	basepair_index[4] = 1;		//	cis	W/W	AU	(3)
	basepair_index[5] = 2;		//	cis	W/W	CA	(4)
	basepair_index[6] = 6;		//	cis	W/W	CC	(5)
	basepair_index[7] = 1;		//	cis	W/W	CG	(6)
	basepair_index[8] = 5;		//	cis	W/W	CU	(7)
	basepair_index[9] = 3;		//	cis	W/W	GA	(8)
	basepair_index[10] = 1;		//	cis	W/W	GC	(9)
	basepair_index[11] = 0;		//	cis	W/W	GG	(10)
	basepair_index[12] = 2;		//	cis	W/W	GU	(11)
	basepair_index[13] = 1;		//	cis	W/W	UA	(12)
	basepair_index[14] = 5;		//	cis	W/W	UC	(13)
	basepair_index[15] = 2;		//	cis	W/W	UG	(14)
	basepair_index[16] = 6;		//	cis	W/W	UU	(15)
	
	//	<*****matrix_2*****>
	basepair_index[17] = 10;	//	trans	W/W	AA	(16)
	basepair_index[18] = 9;		//	trans	W/W	AC	(17)
	basepair_index[19] = 0;		//	trans	W/W	AG	(18)
	basepair_index[20] = 7;		//	trans	W/W	AU	(19)
	basepair_index[21] = 9;		//	trans	W/W	CA	(20)
	basepair_index[22] = 12;	//	trans	W/W	CC	(21)
	basepair_index[23] = 8;		//	trans	W/W	CG	(22)
	basepair_index[24] = 11;	//	trans	W/W	CU	(23)
	basepair_index[25] = 0;		//	trans	W/W	GA	(24)
	basepair_index[26] = 8;		//	trans	W/W	GC	(25)
	basepair_index[27] = 10;	//	trans	W/W	GG	(26)
	basepair_index[28] = 9;		//	trans	W/W	GU	(27)
	basepair_index[29] = 7;		//	trans	W/W	UA	(28)
	basepair_index[30] = 11;	//	trans	W/W	UC	(29)
	basepair_index[31] = 9;		//	trans	W/W	UG	(30)
	basepair_index[32] = 12;	//	trans	W/W	UU	(31)

	//	<*****matrix_3*****>
	basepair_index[33] = 0;		//	cis	W/H	AA	(32)
	basepair_index[34] = 0;		//	cis	W/H	AC	(33)
	basepair_index[35] = 15;	//	cis	W/H	AG	(34)
	basepair_index[36] = 15;	//	cis	W/H	AU	(35)
	basepair_index[37] = 0;		//	cis	W/H	CA	(36)
	basepair_index[38] = 14;	//	cis	W/H	CC	(37)
	basepair_index[39] = 13;	//	cis	W/H	CG	(38)
	basepair_index[40] = 13;	//	cis	W/H	CU	(39)
	basepair_index[41] = 15;	//	cis	W/H	GA	(40)
	basepair_index[42] = 0;		//	cis	W/H	GC	(41)
	basepair_index[43] = 16;	//	cis	W/H	GG	(42)
	basepair_index[44] = 0;		//	cis	W/H	GU	(43)
	basepair_index[45] = 13;	//	cis	W/H	UA	(44)
	basepair_index[46] = 0;		//	cis	W/H	UC	(45)
	basepair_index[47] = 13;	//	cis	W/H	UG	(46)
	basepair_index[48] = 14;	//	cis	W/H	UU	(47)

	//	<*****matrix_4*****>
	basepair_index[49] = 20;	//	trans	W/H	AA	(48)
	basepair_index[50] = 0;		//	trans	W/H	AC	(49)
	basepair_index[51] = 20;	//	trans	W/H	AG	(50)
	basepair_index[52] = 0;		//	trans	W/H	AU	(51)
	basepair_index[53] = 18;	//	trans	W/H	CA	(52)
	basepair_index[54] = 17;	//	trans	W/H	CC	(53)
	basepair_index[55] = 18;	//	trans	W/H	CG	(54)
	basepair_index[56] = 0;		//	trans	W/H	CU	(55)
	basepair_index[57] = 0;		//	trans	W/H	GA	(56)
	basepair_index[58] = 0;		//	trans	W/H	GC	(57)
	basepair_index[59] = 21;	//	trans	W/H	GG	(58)
	basepair_index[60] = 20;	//	trans	W/H	GU	(59)
	basepair_index[61] = 17;	//	trans	W/H	UA	(60)
	basepair_index[62] = 0;		//	trans	W/H	UC	(61)
	basepair_index[63] = 19;	//	trans	W/H	UG	(62)
	basepair_index[64] = 18;	//	trans	W/H	UU	(63)

	//	<*****matrix_5*****>
	basepair_index[65] = 22;	//	cis	W/S	AA	(64)
	basepair_index[66] = 22;	//	cis	W/S	AC	(65)
	basepair_index[67] = 22;	//	cis	W/S	AG	(66)
	basepair_index[68] = 22;	//	cis	W/S	AU	(67)
	basepair_index[69] = 23;	//	cis	W/S	CA	(68)
	basepair_index[70] = 23;	//	cis	W/S	CC	(69)
	basepair_index[71] = 23;	//	cis	W/S	CG	(70)
	basepair_index[72] = 23;	//	cis	W/S	CU	(71)
	basepair_index[73] = 24;	//	cis	W/S	GA	(72)
	basepair_index[74] = 24;	//	cis	W/S	GC	(73)
	basepair_index[75] = 26;	//	cis	W/S	GG	(74)
	basepair_index[76] = 24;	//	cis	W/S	GU	(75)
	basepair_index[77] = 25;	//	cis	W/S	UA	(76)
	basepair_index[78] = 25;	//	cis	W/S	UC	(77)
	basepair_index[79] = 25;	//	cis	W/S	UG	(78)
	basepair_index[80] = 25;	//	cis	W/S	UU	(79)

	//	<*****matrix_6*****>
	basepair_index[81] = 27;	//	trans	W/S	AA	(80)
	basepair_index[82] = 27;	//	trans	W/S	AC	(81)
	basepair_index[83] = 27;	//	trans	W/S	AG	(82)
	basepair_index[84] = 27;	//	trans	W/S	AU	(83)
	basepair_index[85] = 27;	//	trans	W/S	CA	(84)
	basepair_index[86] = 27;	//	trans	W/S	CC	(85)
	basepair_index[87] = 27;	//	trans	W/S	CG	(86)
	basepair_index[88] = 27;	//	trans	W/S	CU	(87)
	basepair_index[89] = 0;		//	trans	W/S	GA	(88)
	basepair_index[90] = 28;	//	trans	W/S	GC	(89)
	basepair_index[91] = 0;		//	trans	W/S	GG	(90)
	basepair_index[92] = 28;	//	trans	W/S	GU	(91)
	basepair_index[93] = 29;	//	trans	W/S	UA	(92)
	basepair_index[94] = 29;	//	trans	W/S	UC	(93)
	basepair_index[95] = 30;	//	trans	W/S	UG	(94)
	basepair_index[96] = 29;	//	trans	W/S	UU	(95)

	//	<*****matrix_7*****>
	basepair_index[97] = 0;		//	cis	H/W	AA	(96)
	basepair_index[98] = 0;		//	cis	H/W	AC	(97)
	basepair_index[99] = 33;	//	cis	H/W	AG	(98)
	basepair_index[100] = 31;	//	cis	H/W	AU	(99)
	basepair_index[101] = 0;	//	cis	H/W	CA	(100)
	basepair_index[102] = 32;	//	cis	H/W	CC	(101)
	basepair_index[103] = 0;	//	cis	H/W	CG	(102)
	basepair_index[104] = 0;	//	cis	H/W	CU	(103)
	basepair_index[105] = 33;	//	cis	H/W	GA	(104)
	basepair_index[106] = 31;	//	cis	H/W	GC	(105)
	basepair_index[107] = 34;	//	cis	H/W	GG	(106)
	basepair_index[108] = 31;	//	cis	H/W	GU	(107)
	basepair_index[109] = 33;	//	cis	H/W	UA	(108)
	basepair_index[110] = 31;	//	cis	H/W	UC	(109)
	basepair_index[111] = 0;	//	cis	H/W	UG	(110)
	basepair_index[112] = 32;	//	cis	H/W	UU	(111)

	//	<*****matrix_8*****>
	basepair_index[113] = 38;	//	trans	H/W	AA	(112)
	basepair_index[114] = 36;	//	trans	H/W	AC	(113)
	basepair_index[115] = 0;	//	trans	H/W	AG	(114)
	basepair_index[116] = 35;	//	trans	H/W	AU	(115)
	basepair_index[117] = 0;	//	trans	H/W	CA	(116)
	basepair_index[118] = 35;	//	trans	H/W	CC	(117)
	basepair_index[119] = 0;	//	trans	H/W	CG	(118)
	basepair_index[120] = 0;	//	trans	H/W	CU	(119)
	basepair_index[121] = 38;	//	trans	H/W	GA	(120)
	basepair_index[122] = 36;	//	trans	H/W	GC	(121)
	basepair_index[123] = 39;	//	trans	H/W	GG	(122)
	basepair_index[124] = 37;	//	trans	H/W	GU	(123)
	basepair_index[125] = 0;	//	trans	H/W	UA	(124)
	basepair_index[126] = 0;	//	trans	H/W	UC	(125)
	basepair_index[127] = 38;	//	trans	H/W	UG	(126)
	basepair_index[128] = 36;	//	trans	H/W	UU	(127)

	//	<*****matrix_9*****>
	basepair_index[129] = 0;	//	cis	H/H	AA	(128)
	basepair_index[130] = 0;	//	cis	H/H	AC	(129)
	basepair_index[131] = 41;	//	cis	H/H	AG	(130)
	basepair_index[132] = 0;	//	cis	H/H	AU	(131)
	basepair_index[133] = 0;	//	cis	H/H	CA	(132)
	basepair_index[134] = 0;	//	cis	H/H	CC	(133)
	basepair_index[135] = 40;	//	cis	H/H	CG	(134)
	basepair_index[136] = 0;	//	cis	H/H	CU	(135)
	basepair_index[137] = 41;	//	cis	H/H	GA	(136)
	basepair_index[138] = 40;	//	cis	H/H	GC	(137)
	basepair_index[139] = 40;	//	cis	H/H	GG	(138)
	basepair_index[140] = 0;	//	cis	H/H	GU	(139)
	basepair_index[141] = 0;	//	cis	H/H	UA	(140)
	basepair_index[142] = 0;	//	cis	H/H	UC	(141)
	basepair_index[143] = 0;	//	cis	H/H	UG	(142)
	basepair_index[144] = 0;	//	cis	H/H	UU	(143)
	
	//	<*****matrix_10*****>
	basepair_index[145] = 42;	//	trans	H/H	AA	(144)
	basepair_index[146] = 42;	//	trans	H/H	AC	(145)
	basepair_index[147] = 43;	//	trans	H/H	AG	(146)
	basepair_index[148] = 43;	//	trans	H/H	AU	(147)
	basepair_index[149] = 42;	//	trans	H/H	CA	(148)
	basepair_index[150] = 0;	//	trans	H/H	CC	(149)
	basepair_index[151] = 42;	//	trans	H/H	CG	(150)
	basepair_index[152] = 43;	//	trans	H/H	CU	(151)
	basepair_index[153] = 43;	//	trans	H/H	GA	(152)
	basepair_index[154] = 42;	//	trans	H/H	GC	(153)
	basepair_index[155] = 44;	//	trans	H/H	GG	(154)
	basepair_index[156] = 0;	//	trans	H/H	GU	(155)
	basepair_index[157] = 43;	//	trans	H/H	UA	(156)
	basepair_index[158] = 43;	//	trans	H/H	UC	(157)
	basepair_index[159] = 0;	//	trans	H/H	UG	(158)
	basepair_index[160] = 0;	//	trans	H/H	UU	(159)

	//	<*****matrix_11*****>
	basepair_index[161] = 45;	//	cis	H/S	AA	(160)
	basepair_index[162] = 45;	//	cis	H/S	AC	(161)
	basepair_index[163] = 45;	//	cis	H/S	AG	(162)
	basepair_index[164] = 45;	//	cis	H/S	AU	(163)
	basepair_index[165] = 45;	//	cis	H/S	CA	(164)
	basepair_index[166] = 45;	//	cis	H/S	CC	(165)
	basepair_index[167] = 45;	//	cis	H/S	CG	(166)
	basepair_index[168] = 45;	//	cis	H/S	CU	(167)
	basepair_index[169] = 45;	//	cis	H/S	GA	(168)
	basepair_index[170] = 0;	//	cis	H/S	GC	(169)
	basepair_index[171] = 45;	//	cis	H/S	GG	(170)
	basepair_index[172] = 0;	//	cis	H/S	GU	(171)
	basepair_index[173] = 46;	//	cis	H/S	UA	(172)
	basepair_index[174] = 45;	//	cis	H/S	UC	(173)
	basepair_index[175] = 45;	//	cis	H/S	UG	(174)
	basepair_index[176] = 45;	//	cis	H/S	UU	(175)

	//	<*****matrix_12*****>
	basepair_index[177] = 47;	//	trans	H/S	AA	(176)
	basepair_index[178] = 47;	//	trans	H/S	AC	(177)
	basepair_index[179] = 47;	//	trans	H/S	AG	(178)
	basepair_index[180] = 47;	//	trans	H/S	AU	(179)
	basepair_index[181] = 47;	//	trans	H/S	CA	(180)
	basepair_index[182] = 47;	//	trans	H/S	CC	(181)
	basepair_index[183] = 0;	//	trans	H/S	CG	(182)
	basepair_index[184] = 47;	//	trans	H/S	CU	(183)
	basepair_index[185] = 0;	//	trans	H/S	GA	(184)
	basepair_index[186] = 0;	//	trans	H/S	GC	(185)
	basepair_index[187] = 48;	//	trans	H/S	GG	(186)
	basepair_index[188] = 0;	//	trans	H/S	GU	(187)
	basepair_index[189] = 48;	//	trans	H/S	UA	(188)
	basepair_index[190] = 0;	//	trans	H/S	UC	(189)
	basepair_index[191] = 48;	//	trans	H/S	UG	(190)
	basepair_index[192] = 0;	//	trans	H/S	UU	(191)

	//	<*****matrix_13*****>
	basepair_index[193] = 49;	//	cis	S/W	AA	(192)
	basepair_index[194] = 50;	//	cis	S/W	AC	(193)
	basepair_index[195] = 51;	//	cis	S/W	AG	(194)
	basepair_index[196] = 52;	//	cis	S/W	AU	(195)
	basepair_index[197] = 49;	//	cis	S/W	CA	(196)
	basepair_index[198] = 50;	//	cis	S/W	CC	(197)
	basepair_index[199] = 51;	//	cis	S/W	CG	(198)
	basepair_index[200] = 52;	//	cis	S/W	CU	(199)
	basepair_index[201] = 49;	//	cis	S/W	GA	(200)
	basepair_index[202] = 50;	//	cis	S/W	GC	(201)
	basepair_index[203] = 53;	//	cis	S/W	GG	(202)
	basepair_index[204] = 52;	//	cis	S/W	GU	(203)
	basepair_index[205] = 49;	//	cis	S/W	UA	(204)
	basepair_index[206] = 50;	//	cis	S/W	UC	(205)
	basepair_index[207] = 51;	//	cis	S/W	UG	(206)
	basepair_index[208] = 52;	//	cis	S/W	UU	(207)

	//	<*****matrix_14*****>
	basepair_index[209] = 54;	//	trans	S/W	AA	(208)
	basepair_index[210] = 54;	//	trans	S/W	AC	(209)
	basepair_index[211] = 0;	//	trans	S/W	AG	(210)
	basepair_index[212] = 56;	//	trans	S/W	AU	(211)
	basepair_index[213] = 54;	//	trans	S/W	CA	(212)
	basepair_index[214] = 54;	//	trans	S/W	CC	(213)
	basepair_index[215] = 55;	//	trans	S/W	CG	(214)
	basepair_index[216] = 56;	//	trans	S/W	CU	(215)
	basepair_index[217] = 54;	//	trans	S/W	GA	(216)
	basepair_index[218] = 54;	//	trans	S/W	GC	(217)
	basepair_index[219] = 0;	//	trans	S/W	GG	(218)
	basepair_index[220] = 57;	//	trans	S/W	GU	(219)
	basepair_index[221] = 54;	//	trans	S/W	UA	(220)
	basepair_index[222] = 54;	//	trans	S/W	UC	(221)
	basepair_index[223] = 55;	//	trans	S/W	UG	(222)
	basepair_index[224] = 56;	//	trans	S/W	UU	(223)

	//	<*****matrix_15*****>
	basepair_index[225] = 58;	//	cis	S/H	AA	(224)
	basepair_index[226] = 58;	//	cis	S/H	AC	(225)
	basepair_index[227] = 58;	//	cis	S/H	AG	(226)
	basepair_index[228] = 59;	//	cis	S/H	AU	(227)
	basepair_index[229] = 58;	//	cis	S/H	CA	(228)
	basepair_index[230] = 58;	//	cis	S/H	CC	(229)
	basepair_index[231] = 0;	//	cis	S/H	CG	(230)
	basepair_index[232] = 58;	//	cis	S/H	CU	(231)
	basepair_index[233] = 58;	//	cis	S/H	GA	(232)
	basepair_index[234] = 58;	//	cis	S/H	GC	(233)
	basepair_index[235] = 58;	//	cis	S/H	GG	(234)
	basepair_index[236] = 58;	//	cis	S/H	GU	(235)
	basepair_index[237] = 58;	//	cis	S/H	UA	(236)
	basepair_index[238] = 58;	//	cis	S/H	UC	(237)
	basepair_index[239] = 0;	//	cis	S/H	UG	(238)
	basepair_index[240] = 58;	//	cis	S/H	UU	(239)

	//	<*****matrix_16*****>
	basepair_index[241] = 60;	//	trans	S/H	AA	(240)
	basepair_index[242] = 60;	//	trans	S/H	AC	(241)
	basepair_index[243] = 0;	//	trans	S/H	AG	(242)
	basepair_index[244] = 61;	//	trans	S/H	AU	(243)
	basepair_index[245] = 60;	//	trans	S/H	CA	(244)
	basepair_index[246] = 60;	//	trans	S/H	CC	(245)
	basepair_index[247] = 0;	//	trans	S/H	CG	(246)
	basepair_index[248] = 0;	//	trans	S/H	CU	(247)
	basepair_index[249] = 60;	//	trans	S/H	GA	(248)
	basepair_index[250] = 0;	//	trans	S/H	GC	(249)
	basepair_index[251] = 61;	//	trans	S/H	GG	(250)
	basepair_index[252] = 61;	//	trans	S/H	GU	(251)
	basepair_index[253] = 60;	//	trans	S/H	UA	(252)
	basepair_index[254] = 60;	//	trans	S/H	UC	(253)
	basepair_index[255] = 0;	//	trans	S/H	UG	(254)
	basepair_index[256] = 0;	//	trans	S/H	UU	(255)

	//	<*****matrix_17*****>
	basepair_index[257] = 62;	//	cis	S/S	AA	(256)
	basepair_index[258] = 62;	//	cis	S/S	AC	(257)
	basepair_index[259] = 62;	//	cis	S/S	AG	(258)
	basepair_index[260] = 62;	//	cis	S/S	AU	(259)
	basepair_index[261] = 62;	//	cis	S/S	CA	(260)
	basepair_index[262] = 62;	//	cis	S/S	CC	(261)
	basepair_index[263] = 62;	//	cis	S/S	CG	(262)
	basepair_index[264] = 62;	//	cis	S/S	CU	(263)
	basepair_index[265] = 62;	//	cis	S/S	GA	(264)
	basepair_index[266] = 62;	//	cis	S/S	GC	(265)
	basepair_index[267] = 62;	//	cis	S/S	GG	(266)
	basepair_index[268] = 62;	//	cis	S/S	GU	(267)
	basepair_index[269] = 62;	//	cis	S/S	UA	(268)
	basepair_index[270] = 62;	//	cis	S/S	UC	(269)
	basepair_index[271] = 62;	//	cis	S/S	UG	(270)
	basepair_index[272] = 62;	//	cis	S/S	UU	(271)

	//	<*****matrix_18*****>
	basepair_index[273] = 63;	//	trans	S/S	AA	(272)
	basepair_index[274] = 63;	//	trans	S/S	AC	(273)
	basepair_index[275] = 63;	//	trans	S/S	AG	(274)
	basepair_index[276] = 63;	//	trans	S/S	AU	(275)
	basepair_index[277] = 0;	//	trans	S/S	CA	(276)
	basepair_index[278] = 0;	//	trans	S/S	CC	(277)
	basepair_index[279] = 0;	//	trans	S/S	CG	(278)
	basepair_index[280] = 0;	//	trans	S/S	CU	(279)
	basepair_index[281] = 64;	//	trans	S/S	GA	(280)
	basepair_index[282] = 64;	//	trans	S/S	GC	(281)
	basepair_index[283] = 64;	//	trans	S/S	GG	(282)
	basepair_index[284] = 64;	//	trans	S/S	GU	(283)
	basepair_index[285] = 0;	//	trans	S/S	UA	(284)
	basepair_index[286] = 0;	//	trans	S/S	UC	(285)
	basepair_index[287] = 0;	//	trans	S/S	UG	(286)
	basepair_index[288] = 0;	//	trans	S/S	UU	(287)
	
	
	return;
}
