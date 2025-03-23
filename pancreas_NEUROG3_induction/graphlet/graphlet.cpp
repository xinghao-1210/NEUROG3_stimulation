#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <map>
#include <vector>
#include <set>
#include <algorithm>

using namespace std;
typedef long long int64;


/**
 *  Original code by Anida Sarajlic
 *
 *  Updated by Noel Malod-Dognin
 *  - removed searches in predecessors and successor: O(d) -> cst; requieres additional O(n^2) space
 *  - to be added: parallel version
 *
 */


//main program

int main ( int argc, char *argv[] ){

	vector<pair<int64, int64> > edgelista;
	vector<string> dictionary;
	string x1, x2, line;
	vector<vector<int64> > pred, succ, signatures;
	vector<vector<bool> >  adj_matrix;
	vector<int64> graphlets;
	int nb_nodes(0);
	int ind1, ind2;
	map<string,int> node_sorter;

	//reading in the edgelist (clean self loops) and creating edgelist of nodes (nodes will be named using numbering from 0 to (numofnodes-1) so we create a dictionary for nodenames too)
	if (argc != 2){
		cout<<"Usage: "<< argv[0] <<" <input_file_name> "<<endl;
	}
    else{
		ifstream data_file(argv[1]);
		if (!data_file.is_open()){
            cout<<"Could not open file\n"<<endl;
		}
		else{
			while(getline(data_file,line)){
				stringstream string_linije(line);
				string_linije >> x1 >> x2;
				if (x1 != x2){
					map<string,int>::iterator fx1(node_sorter.find(x1));
					if(fx1!=node_sorter.end()){
						ind1=fx1->second;
					}
					else{
						node_sorter[x1]=nb_nodes;
						ind1=nb_nodes;
						++nb_nodes;
						dictionary.push_back(x1);
					}
					map<string,int>::iterator fx2(node_sorter.find(x2));
					if(fx2!=node_sorter.end()){
						ind2=fx2->second;
					}
					else{
						node_sorter[x2]=nb_nodes;
						ind2=nb_nodes;
						++nb_nodes;
						dictionary.push_back(x2);
					}
					edgelista.push_back(make_pair(ind1,ind2));
				}
			}
		}
		data_file.close();
	}

	//creating a container (vector of vectors) for nodes' predecessors, and similar container for nodes' succesors : these two are length of the  dictionary, and 'inner' vectors are still unknown length. Elements in succ[i] (same goes for pred[i]) have to be different, and this ensures that the duplicated edges are removed. This is important for updateing orbits later - not to duplicate them. Also node will not be appearing among its own sucessors or predecessors as selfloops are removed when creating edgelist.
	// Also, a similar container of vectors signatures is created, length is the same as node dictionary but 'inner vectors are length 129(there are 129 orbits from 0 to 128). Vector 'graphlets' - size 40, will store number of graphlets in the network (graphlet0 to graphlet 39).

	pred.resize(dictionary.size());
	succ.resize(dictionary.size());
	signatures.resize(dictionary.size());
	adj_matrix.resize(dictionary.size());
	for(size_t i(0); i<dictionary.size(); ++i){
		adj_matrix[i].resize(dictionary.size());
		fill(adj_matrix[i].begin(), adj_matrix[i].end(), false);
	}
	graphlets.resize(40);
	for (vector<vector<int64> >::iterator it = signatures.begin(); it !=signatures.end(); ++it) {
		(*it).resize(129);
		fill((*it).begin(), (*it).end(), 0);
	}
	for (vector<pair<int64,int64> >::iterator it = edgelista.begin(); it != edgelista.end(); ++it){
		if(adj_matrix[(*it).first][(*it).second]==false){
			adj_matrix[(*it).first][(*it).second] = true;
			succ[(*it).first].push_back((*it).second);
			pred[(*it).second].push_back((*it).first);
		}
	}



//PROBA PRint64ANJA SUCCESORA I PREDECESSORA
/*
for (int64 i(0); i<succ.size();++i){
	cout<<"node "<<i<<" has succesors: "<<endl;
		for (int64 j(0); j<succ[i].size();++j) { cout<<succ[i][j]<<endl;}
}
for (int64 i(0); i<pred.size();++i){
	cout<<"node "<<i<<" has predeccesors: "<<endl;
		for (int64 j(0); j<pred[i].size();++j) { cout<<pred[i][j]<<endl;}
}
*/
cout<<"pred and succ matrixes over"<<endl;

//UPDATE-ING ORBITS AND GRAPHLETS -  some graphlets/orbits are grouped together because for them the for loops and neighbourhood checks can be reused. Ofcourse they can be grouped differently, maybe even more efficiently, but this is one of the ways to do it...Also, there are many options how to 'walk' through each graphlet, which node is the 'first' to be checked. It is IMPORTANT to take into account how this was done when correcting for overcounts later in the end!

for (int64 i(0); i<dictionary.size(); ++i){

	//orb0, orb1, faster counting
	signatures[i][0]=succ[i].size();
	signatures[i][1]=pred[i].size();
	graphlets[0]+=succ[i].size();

	//orb2,3,4
	for (vector<int64> ::iterator it1 = pred[i].begin(); it1 !=pred[i].end(); ++it1) {
		for (vector<int64> ::iterator it2 = succ[i].begin(); it2 !=succ[i].end(); ++it2) {
			//if ((*it1)!=(*it2) and (is_it_in_vector(*it1,succ[*it2]) == 0) and (is_it_in_vector(*it1,pred[*it2]) == 0)) {
			if ((*it1)!=(*it2) and (!adj_matrix[*it2][*it1]) and (!adj_matrix[*it1][*it2])) {
				++signatures[*it1][2];
				++signatures[*it2][4];
				++signatures[i][3];
				++graphlets[1];

				//orb 13,14,15,16 and orb 53,54,55,56 and orb 39 and orb 81,82,83,84 adn orb 40,41,42,43 and orb 105,106,107,108 and orb 109,110,111,112
				for (vector<int64> ::iterator it3 = succ[*it2].begin(); it3 !=succ[*it2].end(); ++it3) {
					//if ((*it3)!=(*it1) and (*it3)!=(i) and (is_it_in_vector(*it3,succ[*it1]) == 0) and (is_it_in_vector(*it3,pred[*it1]) == 0) and (is_it_in_vector(*it3,succ[i]) == 0) and 	(is_it_in_vector(*it3,pred[i]) == 0)) {
					if ((*it3)!=(*it1) and (*it3)!=(i) and (!adj_matrix[*it1][*it3]) and (!adj_matrix[*it3][*it1]) and (!adj_matrix[i][*it3]) and (!adj_matrix[*it3][i])) {
						++signatures[*it1][13];
						++signatures[*it2][15];
						++signatures[i][14];
						++signatures[*it3][16];
						++graphlets[6];
					}
					//if ((*it3)!=(*it1) and (*it3)!=(i) and (is_it_in_vector(*it3,succ[*it1]) == 0) and (is_it_in_vector(*it3,pred[*it1]) == 0)  and (is_it_in_vector(*it3,pred[i]) == 1)) {
					if ((*it3)!=(*it1) and (*it3)!=(i) and (!adj_matrix[*it1][*it3]) and (!adj_matrix[*it3][*it1]) and (adj_matrix[*it3][i]) ) {
						++signatures[*it1][56];
						++signatures[*it2][53];
						++signatures[i][55];
						++signatures[*it3][54];
						++graphlets[19];
					}
					//if ((*it3)!=(*it1) and (*it3)!=(i) and (is_it_in_vector(*it3,pred[*it1]) == 1) and (is_it_in_vector(*it3,succ[i]) == 0) and (is_it_in_vector(*it3,pred[i]) == 0)) {
					if ((*it3)!=(*it1) and (*it3)!=(i) and (adj_matrix[*it3][*it1]) and (!adj_matrix[i][*it3]) and (!adj_matrix[*it3][i])) {
						++signatures[*it1][39];
						++signatures[*it2][39];
						++signatures[i][39];
						++signatures[*it3][39];
						++graphlets[14];
					}
					//if ((*it3)!=(*it1) and (*it3)!=(i) and (is_it_in_vector(*it3,pred[*it1]) == 1) and (is_it_in_vector(*it3,pred[i]) == 1)) {
					if ((*it3)!=(*it1) and (*it3)!=(i) and (adj_matrix[*it3][*it1]) and (adj_matrix[*it3][i])) {
						++signatures[*it1][82];
						++signatures[*it2][83];
						++signatures[i][81];
						++signatures[*it3][84];
						++graphlets[26];
					}
					//if ((*it3)!=(*it1) and (*it3)!=(i) and (is_it_in_vector(*it3,succ[*it1]) == 1) and (is_it_in_vector(*it3,succ[i]) == 0) and (is_it_in_vector(*it3,pred[i]) == 0)) {
					if ((*it3)!=(*it1) and (*it3)!=(i) and (adj_matrix[*it1][*it3]) and (!adj_matrix[i][*it3]) and (!adj_matrix[*it3][i])) {
						++signatures[*it1][42];
						++signatures[*it2][41];
						++signatures[i][40];
						++signatures[*it3][43];
						++graphlets[15];
					}
					//if ((*it3)!=(*it1) and (*it3)!=(i) and (is_it_in_vector(*it3,succ[*it1]) == 1) and (is_it_in_vector(*it3,succ[i]) == 1)) {
					if ((*it3)!=(*it1) and (*it3)!=(i) and (adj_matrix[*it1][*it3]) and (adj_matrix[i][*it3])) {
						++signatures[*it1][107];
						++signatures[*it2][106];
						++signatures[i][108];
						++signatures[*it3][105];
						++graphlets[33];
					}
					//if ((*it3)!=(*it1) and (*it3)!=(i) and (is_it_in_vector(*it3,succ[*it1]) == 1)  and (is_it_in_vector(*it3,pred[i]) == 1)) {
					if ((*it3)!=(*it1) and (*it3)!=(i) and (adj_matrix[*it1][*it3])  and (adj_matrix[*it3][i])) {
						++signatures[*it1][111];
						++signatures[*it2][110];
						++signatures[i][112];
						++signatures[*it3][109];
						++graphlets[34];
					}
				}

			//orb 17,18,19,20 and orb 77,78,79,80 and orb 69,70,71,72 and 113,114,115,116
			for (vector<int64> ::iterator it4 = pred[*it2].begin(); it4 !=pred[*it2].end(); ++it4) {
				//if ((*it4)!=(*it1) and (*it4)!=(i) and (is_it_in_vector(*it4,succ[*it1]) == 0) and (is_it_in_vector(*it4,pred[*it1]) == 0) and (is_it_in_vector(*it4,succ[i]) == 0) and (is_it_in_vector(*it4,pred[i]) == 0)) {
				if ((*it4)!=(*it1) and (*it4)!=(i) and (!adj_matrix[*it1][*it4]) and (!adj_matrix[*it4][*it1]) and (!adj_matrix[i][*it4]) and (!adj_matrix[*it4][i])) {
					++signatures[*it1][17];
					++signatures[*it2][19];
					++signatures[i][18];
					++signatures[*it4][20];
					++graphlets[7];
				}
				//if ((*it4)!=(*it1) and (*it4)!=(i) and (is_it_in_vector(*it4,succ[*it1]) == 0) and (is_it_in_vector(*it4,pred[*it1]) == 0) and  (is_it_in_vector(*it4,pred[i]) == 1)) {
				if ((*it4)!=(*it1) and (*it4)!=(i) and (!adj_matrix[*it1][*it4]) and (!adj_matrix[*it4][*it1]) and (adj_matrix[*it4][i])) {
					++signatures[*it1][80];
					++signatures[*it2][78];
					++signatures[i][79];
					++signatures[*it4][77];
					++graphlets[25];
				}
				//if ((*it4)!=(*it1) and (*it4)!=(i) and (is_it_in_vector(*it4,succ[*it1]) == 0) and (is_it_in_vector(*it4,pred[*it1]) == 0) and  (is_it_in_vector(*it4,succ[i]) == 1)) {
				if ((*it4)!=(*it1) and (*it4)!=(i) and (!adj_matrix[*it1][*it4]) and (!adj_matrix[*it4][*it1]) and (adj_matrix[i][*it4])) {
					++signatures[*it1][72];
					++signatures[*it2][70];
					++signatures[i][71];
					++signatures[*it4][69];
					++graphlets[23];
				}
				//if ((*it4)!=(*it1) and (*it4)!=(i) and (is_it_in_vector(*it4,succ[*it1]) == 1)  and  (is_it_in_vector(*it4,succ[i]) == 1)) {
				if ((*it4)!=(*it1) and (*it4)!=(i) and (adj_matrix[*it1][*it4])  and  (adj_matrix[i][*it4])) {
					++signatures[*it1][113];
					++signatures[*it2][115];
					++signatures[i][114];
					++signatures[*it4][116];
					++graphlets[35];
				}
			}
		}
	}
}

//orb5,6
for (vector<int64> ::iterator it1 = succ[i].begin(); it1 !=succ[i].end(); ++it1) {
	for (vector<int64> ::iterator it2 = succ[i].begin(); it2 !=succ[i].end(); ++it2) {
		//if ((*it1)!=(*it2) and (is_it_in_vector(*it1,succ[*it2]) == 0) and (is_it_in_vector(*it1,pred[*it2]) == 0)) {
		if ((*it1)!=(*it2) and (!adj_matrix[*it2][*it1]) and (!adj_matrix[*it1][*it2])) {
			++signatures[*it1][5];
			++signatures[*it2][5];
			++signatures[i][6];
			++graphlets[2];

			//orb 21,22,23,24 and orb 73,74,75,76 and orb 97,98,99,100 and orb 101,102,103,104 and orb 44,45 and orb 94,95,96
			for (vector<int64> ::iterator it3 = pred[*it1].begin(); it3 !=pred[*it1].end(); ++it3) {
				//if ((*it3)!=(*it2) and (*it3)!=(i) and (is_it_in_vector(*it3,succ[*it2]) == 0) and (is_it_in_vector(*it3,pred[*it2]) == 0) and (is_it_in_vector(*it3,succ[i]) == 0) and (is_it_in_vector(*it3,pred[i]) == 0)) {
				if ((*it3)!=(*it2) and (*it3)!=(i) and (!adj_matrix[*it2][*it3]) and (!adj_matrix[*it3][*it2]) and (!adj_matrix[i][*it3]) and (!adj_matrix[*it3][i])) {
					++signatures[*it1][22];
					++signatures[*it2][24];
					++signatures[i][23];
					++signatures[*it3][21];
                                        ++graphlets[8];
				}
				//if ((*it3)!=(*it2) and (*it3)!=(i) and (is_it_in_vector(*it3,succ[*it2]) == 0) and (is_it_in_vector(*it3,pred[*it2]) == 0) and (is_it_in_vector(*it3,pred[i]) == 1)) {
				if ((*it3)!=(*it2) and (*it3)!=(i) and (!adj_matrix[*it2][*it3]) and (!adj_matrix[*it3][*it2]) and (adj_matrix[*it3][i])) {
					++signatures[*it1][74];
					++signatures[*it2][76];
					++signatures[i][75];
					++signatures[*it3][73];
                                        ++graphlets[24];
				}
				//if ((*it3)!=(*it2) and (*it3)!=(i) and (is_it_in_vector(*it3,succ[*it2]) == 1)  and (is_it_in_vector(*it3,succ[i]) == 1)) {
				if ((*it3)!=(*it2) and (*it3)!=(i) and (adj_matrix[*it2][*it3]) and (adj_matrix[i][*it3])) {
					++signatures[*it1][99];
					++signatures[*it2][98];
					++signatures[i][100];
					++signatures[*it3][97];
                                        ++graphlets[31];
				}
				//if ((*it3)!=(*it2) and (*it3)!=(i) and (is_it_in_vector(*it3,succ[*it2]) == 1)  and (is_it_in_vector(*it3,pred[i]) == 1)) {
				if ((*it3)!=(*it2) and (*it3)!=(i) and (adj_matrix[*it2][*it3]) and (adj_matrix[*it3][i])) {
					++signatures[*it1][103];
					++signatures[*it2][102];
					++signatures[i][104];
					++signatures[*it3][101];
                                        ++graphlets[32];
				}
				//if ((*it3)!=(*it2) and (*it3)!=(i) and (is_it_in_vector(*it3,pred[*it2]) == 1)  and (is_it_in_vector(*it3,succ[i]) == 0) and (is_it_in_vector(*it3,pred[i]) == 0)) {
				if ((*it3)!=(*it2) and (*it3)!=(i) and (adj_matrix[*it3][*it2])  and (!adj_matrix[i][*it3]) and (!adj_matrix[*it3][i])) {
					++signatures[*it1][44];
					++signatures[*it2][44];
					++signatures[i][45];
					++signatures[*it3][45];
                                        ++graphlets[16];
				}
				//if ((*it3)!=(*it2) and (*it3)!=(i) and (is_it_in_vector(*it3,pred[*it2]) == 1)  and (is_it_in_vector(*it3,succ[i]) == 1)) {
				if ((*it3)!=(*it2) and (*it3)!=(i) and (adj_matrix[*it3][*it2])  and (adj_matrix[i][*it3])) {
					++signatures[*it1][95];
					++signatures[*it2][95];
					++signatures[i][96];
					++signatures[*it3][94];
                                        ++graphlets[30];
				}
			}

			//orb 25,26,27,28 and orb 49,50,51,52 and orb 65,66,67,68 and orb 46,47,48 and orb 85,86,87 and 88,89,90
			for (vector<int64> ::iterator it4 = succ[*it1].begin(); it4 !=succ[*it1].end(); ++it4) {
				//if ((*it4)!=(*it2) and (*it4)!=(i) and (is_it_in_vector(*it4,succ[*it2]) == 0) and (is_it_in_vector(*it4,pred[*it2]) == 0) and (is_it_in_vector(*it4,succ[i]) == 0) and (is_it_in_vector(*it4,pred[i]) == 0)) {
				if ((*it4)!=(*it2) and (*it4)!=(i) and (!adj_matrix[*it2][*it4]) and (!adj_matrix[*it4][*it2]) and (!adj_matrix[i][*it4]) and (!adj_matrix[*it4][i])) {
					++signatures[*it1][26];
					++signatures[*it2][28];
					++signatures[i][27];
					++signatures[*it4][25];
					++graphlets[9];
				}
				//if ((*it4)!=(*it2) and (*it4)!=(i) and (is_it_in_vector(*it4,succ[*it2]) == 0) and (is_it_in_vector(*it4,pred[*it2]) == 0)  and (is_it_in_vector(*it4,pred[i]) == 1)) {
				if ((*it4)!=(*it2) and (*it4)!=(i) and (!adj_matrix[*it2][*it4]) and (!adj_matrix[*it4][*it2])  and (adj_matrix[*it4][i])) {
					++signatures[*it1][49];
					++signatures[*it2][52];
					++signatures[i][51];
					++signatures[*it4][50];
					++graphlets[18];
				}
				//if ((*it4)!=(*it2) and (*it4)!=(i) and (is_it_in_vector(*it4,succ[*it2]) == 0) and (is_it_in_vector(*it4,pred[*it2]) == 0)  and (is_it_in_vector(*it4,succ[i]) == 1)) {
				if ((*it4)!=(*it2) and (*it4)!=(i) and (!adj_matrix[*it2][*it4]) and (!adj_matrix[*it4][*it2]) and (adj_matrix[i][*it4])) {
					++signatures[*it1][65];
					++signatures[*it2][68];
					++signatures[i][67];
					++signatures[*it4][66];
					++graphlets[22];
				}
				//if ((*it4)!=(*it2) and (*it4)!=(i) and (is_it_in_vector(*it4,succ[*it2]) == 1)  and (is_it_in_vector(*it4,succ[i]) == 0) and (is_it_in_vector(*it4,pred[i]) == 0)) {
				if ((*it4)!=(*it2) and (*it4)!=(i) and (adj_matrix[*it2][*it4])  and (!adj_matrix[i][*it4]) and (!adj_matrix[*it4][i])) {
					++signatures[*it1][46];
					++signatures[*it2][46];
					++signatures[i][48];
					++signatures[*it4][47];
					++graphlets[17];
				}
				//if ((*it4)!=(*it2) and (*it4)!=(i) and (is_it_in_vector(*it4,succ[*it2]) == 1)  and (is_it_in_vector(*it4,succ[i]) == 1)) {
				if ((*it4)!=(*it2) and (*it4)!=(i) and (adj_matrix[*it2][*it4])  and (adj_matrix[i][*it4])) {
					++signatures[*it1][86];
					++signatures[*it2][86];
					++signatures[i][87];
					++signatures[*it4][85];
					++graphlets[27];
				}
				//if ((*it4)!=(*it2) and (*it4)!=(i) and (is_it_in_vector(*it4,succ[*it2]) == 1)  and (is_it_in_vector(*it4,pred[i]) == 1)) {
				if ((*it4)!=(*it2) and (*it4)!=(i) and (adj_matrix[*it2][*it4]) and (adj_matrix[*it4][i])) {
					++signatures[*it1][89];
					++signatures[*it2][89];
					++signatures[i][90];
					++signatures[*it4][88];
					++graphlets[28];
				}
			}

			//orb 31,32
			for (vector<int64> ::iterator it5 = succ[i].begin(); it5 !=succ[i].end(); ++it5) {
				//if ((*it5)!=(*it2) and (*it5)!=(*it1) and (is_it_in_vector(*it5,succ[*it2]) == 0) and (is_it_in_vector(*it5,pred[*it2]) == 0) and (is_it_in_vector(*it5,succ[*it1]) == 0) and (is_it_in_vector(*it5,pred[*it1]) == 0)) {
				if ((*it5)!=(*it2) and (*it5)!=(*it1) and (!adj_matrix[*it2][*it5]) and (!adj_matrix[*it5][*it2]) and (!adj_matrix[*it1][*it5]) and (!adj_matrix[*it5][*it1])) {
					++signatures[*it1][31];
					++signatures[*it2][31];
					++signatures[i][32];
					++signatures[*it5][31];
					++graphlets[11];

				}
			}

			//orb 33,34,35
			for (vector<int64> ::iterator it6 = pred[i].begin(); it6 !=pred[i].end(); ++it6) {
				//if ((*it6)!=(*it2) and (*it6)!=(*it1) and (is_it_in_vector(*it6,succ[*it2]) == 0) and (is_it_in_vector(*it6,pred[*it2]) == 0) and (is_it_in_vector(*it6,succ[*it1]) == 0) and (is_it_in_vector(*it6,pred[*it1]) == 0)) {
				if ((*it6)!=(*it2) and (*it6)!=(*it1) and (!adj_matrix[*it2][*it6]) and (!adj_matrix[*it6][*it2]) and (!adj_matrix[*it1][*it6]) and (!adj_matrix[*it6][*it1])) {
					++signatures[*it1][35];
					++signatures[*it2][35];
					++signatures[i][34];
					++signatures[*it6][33];
					++graphlets[12];
				}
			}
		}
	}
}

//orb7,8
for (vector<int64> ::iterator it1 = pred[i].begin(); it1 !=pred[i].end(); ++it1) {
	for (vector<int64> ::iterator it2 = pred[i].begin(); it2 !=pred[i].end(); ++it2) {
		//if ((*it1)!=(*it2) and (is_it_in_vector(*it1,succ[*it2]) == 0) and (is_it_in_vector(*it1,pred[*it2]) == 0)) {
		if ((*it1)!=(*it2) and (!adj_matrix[*it2][*it1]) and (!adj_matrix[*it1][*it2])) {
			++signatures[*it1][7];
			++signatures[*it2][7];
			++signatures[i][8];
                        ++graphlets[3];
			//orb 29,30 and orb 61,62,63,64

			for (vector<int64> ::iterator it3 = pred[i].begin(); it3 !=pred[i].end(); ++it3) {
				//if ((*it3)!=(*it2) and (*it3)!=(*it1) and (is_it_in_vector(*it3,succ[*it2]) == 0) and (is_it_in_vector(*it3,pred[*it2]) == 0) and (is_it_in_vector(*it3,succ[*it1]) == 0) and (is_it_in_vector(*it3,pred[*it1]) == 0)) {
				if ((*it3)!=(*it2) and (*it3)!=(*it1) and (!adj_matrix[*it2][*it3]) and (!adj_matrix[*it3][*it2]) and (!adj_matrix[*it1][*it3]) and (!adj_matrix[*it3][*it1])) {
					++signatures[*it1][29];
					++signatures[*it2][29];
					++signatures[i][30];
					++signatures[*it3][29];
					++graphlets[10];

				}
				//if ((*it3)!=(*it2) and (*it3)!=(*it1) and  (is_it_in_vector(*it3,pred[*it2]) == 1) and (is_it_in_vector(*it3,succ[*it1]) == 0) and (is_it_in_vector(*it3,pred[*it1]) == 0)) {
				if ((*it3)!=(*it2) and (*it3)!=(*it1) and  (adj_matrix[*it3][*it2]) and (!adj_matrix[*it1][*it3]) and (!adj_matrix[*it3][*it1])) {
					++signatures[*it1][64];
					++signatures[*it2][62];
					++signatures[i][63];
					++signatures[*it3][61];
					++graphlets[21];

				}
			}

			//orb 36,37,38 and orb 91,92,93
			for (vector<int64> ::iterator it4 = succ[i].begin(); it4 !=succ[i].end(); ++it4) {
				//if ((*it4)!=(*it2) and (*it4)!=(*it1) and (is_it_in_vector(*it4,succ[*it2]) == 0) and (is_it_in_vector(*it4,pred[*it2]) == 0) and (is_it_in_vector(*it4,succ[*it1]) == 0) and (is_it_in_vector(*it4,pred[*it1]) == 0)) {
				if ((*it4)!=(*it2) and (*it4)!=(*it1) and (!adj_matrix[*it2][*it4]) and (!adj_matrix[*it4][*it2]) and (!adj_matrix[*it1][*it4]) and (!adj_matrix[*it4][*it1])) {
					++signatures[*it1][38];
					++signatures[*it2][38];
					++signatures[i][37];
					++signatures[*it4][36];
					++graphlets[13];
				}
				//if ((*it4)!=(*it2) and (*it4)!=(*it1) and (is_it_in_vector(*it4,succ[*it2]) == 1)  and (is_it_in_vector(*it4,succ[*it1]) == 1)) {
				if ((*it4)!=(*it2) and (*it4)!=(*it1) and (adj_matrix[*it2][*it4])  and (adj_matrix[*it1][*it4])) {
					++signatures[*it1][92];
					++signatures[*it2][92];
					++signatures[i][93];
					++signatures[*it4][91];
					++graphlets[29];
				}
			}

		}
	}

}

//orb9
for (vector<int64> ::iterator it1 = pred[i].begin(); it1 !=pred[i].end(); ++it1) {
	for (vector<int64> ::iterator it2 = succ[i].begin(); it2 !=succ[i].end(); ++it2) {
		//if ((*it1)!=(*it2) and (is_it_in_vector(*it1,succ[*it2]) == 1)) {
		if ((*it1)!=(*it2) and (adj_matrix[*it2][*it1])) {
			++signatures[*it1][9];
			++signatures[*it2][9];
			++signatures[i][9];
			++graphlets[4];
		}
	}
}

//orb10,11,12
for (vector<int64> ::iterator it1 = pred[i].begin(); it1 !=pred[i].end(); ++it1) {
	for (vector<int64> ::iterator it2 = pred[i].begin(); it2 !=pred[i].end(); ++it2) {
		//if ((*it1)!=(*it2) and (is_it_in_vector(*it1,pred[*it2]) == 1)) {
		if ((*it1)!=(*it2) and (adj_matrix[*it1][*it2])) {
			++signatures[*it1][11];
			++signatures[*it2][12];
			++signatures[i][10];
			++graphlets[5];

			//orbits 57,58,59,60 and orb 121,122,123,124 and orb 119,120 and orb 125,126,127,128
			for (vector<int64> ::iterator it3 = succ[i].begin(); it3 !=succ[i].end(); ++it3) {
				//if ((*it3)!=(*it2) and (*it3)!=(*it1) and (is_it_in_vector(*it3,succ[*it2]) == 0) and (is_it_in_vector(*it3,pred[*it2]) == 0) and (is_it_in_vector(*it3,succ[*it1]) == 0) and (is_it_in_vector(*it3,pred[*it1]) == 0)) {
				if ((*it3)!=(*it2) and (*it3)!=(*it1) and (!adj_matrix[*it2][*it3]) and (!adj_matrix[*it3][*it2]) and (!adj_matrix[*it1][*it3]) and (!adj_matrix[*it3][*it1])) {
					++signatures[*it1][57];
					++signatures[*it2][58];
					++signatures[i][59];
					++signatures[*it3][60];
					++graphlets[20];
				}
				//if ((*it3)!=(*it2) and (*it3)!=(*it1) and (is_it_in_vector(*it3,pred[*it2]) == 1) and (is_it_in_vector(*it3,pred[*it1]) == 1)) {
				if ((*it3)!=(*it2) and (*it3)!=(*it1) and (adj_matrix[*it3][*it2]) and (adj_matrix[*it3][*it1])) {
					++signatures[*it1][124];
					++signatures[*it2][123];
					++signatures[i][122];
					++signatures[*it3][121];
					++graphlets[38];
				}
				//if ((*it3)!=(*it2) and (*it3)!=(*it1) and (is_it_in_vector(*it3,pred[*it2]) == 1) and (is_it_in_vector(*it3,succ[*it1]) == 1)) {
				if ((*it3)!=(*it2) and (*it3)!=(*it1) and (adj_matrix[*it3][*it2]) and (adj_matrix[*it1][*it3])) {
					++signatures[*it1][120];
					++signatures[*it2][119];
					++signatures[i][119];
					++signatures[*it3][119];
					++graphlets[37];
				}
				//if ((*it3)!=(*it2) and (*it3)!=(*it1) and (is_it_in_vector(*it3,succ[*it2]) == 1) and (is_it_in_vector(*it3,succ[*it1]) == 1)) {
				if ((*it3)!=(*it2) and (*it3)!=(*it1) and (adj_matrix[*it2][*it3]) and (adj_matrix[*it1][*it3])) {
					++signatures[*it1][128];
					++signatures[*it2][126];
					++signatures[i][127];
					++signatures[*it3][125];
					++graphlets[39];
				}
			}
			//orb 117,118
			for (vector<int64> ::iterator it4 = pred[i].begin(); it4 !=pred[i].end(); ++it4) {
				//if ((*it4)!=(*it2) and (*it4)!=(*it1) and  (is_it_in_vector(*it4,pred[*it1]) == 1) and (is_it_in_vector(*it4,succ[*it2]) == 1)) {
				if ((*it4)!=(*it2) and (*it4)!=(*it1) and  (adj_matrix[*it4][*it1]) and (adj_matrix[*it2][*it4])) {
					++signatures[*it1][117];
					++signatures[*it2][117];
					++signatures[i][118];
					++signatures[*it4][117];
					++graphlets[36];
				}
			}

		}
	}
}
//cout<<"node "<<i<<" over"<<endl;
}

//CORRECTING OVERCOUNTS - I also doublechecked how I walked through each graphlet in the loops above, just in case it induced aditional overcounts (besides the intuitive value that can be derived observing automorphism orbits). If there were any uneccesary overcounts, this is corrected here too.

for (int64 i(0); i<dictionary.size(); ++i){
	signatures[i][5]=signatures[i][5]/2;
	signatures[i][6]=signatures[i][6]/2;
	signatures[i][7]=signatures[i][7]/2;
	signatures[i][8]=signatures[i][8]/2;
	signatures[i][9]=signatures[i][9]/3;
	signatures[i][29]=signatures[i][29]/6;
	signatures[i][30]=signatures[i][30]/6;
	signatures[i][31]=signatures[i][31]/6;
	signatures[i][32]=signatures[i][32]/6;
	signatures[i][33]=signatures[i][33]/2;
	signatures[i][34]=signatures[i][34]/2;
	signatures[i][35]=signatures[i][35]/2;
	signatures[i][36]=signatures[i][36]/2;
	signatures[i][37]=signatures[i][37]/2;
	signatures[i][38]=signatures[i][38]/2;
	signatures[i][39]=signatures[i][39]/4;
	signatures[i][44]=signatures[i][44]/4;
	signatures[i][45]=signatures[i][45]/4;
	signatures[i][46]=signatures[i][46]/2;
	signatures[i][47]=signatures[i][47]/2;
	signatures[i][48]=signatures[i][48]/2;
	signatures[i][94]=signatures[i][94]/2;
	signatures[i][95]=signatures[i][95]/2;
	signatures[i][96]=signatures[i][96]/2;
	signatures[i][85]=signatures[i][85]/2;
	signatures[i][86]=signatures[i][86]/2;
	signatures[i][87]=signatures[i][87]/2;
	signatures[i][88]=signatures[i][88]/2;
	signatures[i][89]=signatures[i][89]/2;
	signatures[i][90]=signatures[i][90]/2;
	signatures[i][91]=signatures[i][91]/2;
	signatures[i][92]=signatures[i][92]/2;
	signatures[i][93]=signatures[i][93]/2;
	signatures[i][119]=signatures[i][119]/3;
	signatures[i][120]=signatures[i][120]/3;
	signatures[i][117]=signatures[i][117]/3;
	signatures[i][118]=signatures[i][118]/3;

}
graphlets[2]=graphlets[2]/2;
graphlets[3]=graphlets[3]/2;
graphlets[4]=graphlets[4]/3;
graphlets[10]=graphlets[10]/6;
graphlets[11]=graphlets[11]/6;
graphlets[12]=graphlets[12]/2;
graphlets[13]=graphlets[13]/2;
graphlets[14]=graphlets[14]/4;
graphlets[16]=graphlets[16]/4;
graphlets[17]=graphlets[17]/2;
graphlets[30]=graphlets[30]/2;
graphlets[27]=graphlets[27]/2;
graphlets[28]=graphlets[28]/2;
graphlets[29]=graphlets[29]/2;
graphlets[37]=graphlets[37]/3;
graphlets[36]=graphlets[36]/3;

//outputs


    //output signatures

    ofstream out_file1;
    out_file1.open(string(string(argv[1])+".signatures.txt").c_str());
    for (vector<vector<int64> >::iterator it = signatures.begin(); it != signatures.end(); ++it){
    	for (int64 i(0); i<129; ++i){
        	out_file1<<(*it)[i]<<" ";
        }
        out_file1<<endl;
    }


    //output dictionary (could only print the second column since the first one just corresponds to ordered numbers starting from 0... but I left it like this, just to print out the dictionary as it is)

    ofstream out_file2;
    out_file2.open(string(string(argv[1])+".dictionary.txt").c_str());
    for (int64 i(0); i<dictionary.size(); ++i)
       out_file2<<i<<" "<<dictionary[i]<<endl;


    //output graphlet counts
    ofstream out_file3;
    out_file3.open(string(string(argv[1])+".graphletcounts.txt").c_str()) ;
    for (int64 i(0); i<graphlets.size();++i)
       out_file3<<i<<" "<<graphlets[i]<<endl;
    out_file1.close();
    out_file2.close();
    out_file3.close();

//At one stage I also had all the redundancy formulas checked in the end of the script (it helps with testing too - in case we only use networks without antiparralel pairs of arcs), but this was creating lots of excess output so I removed that at one point...

//cout<<"INPUT: a text file representing network edgelist (entry A	B in a row represents a directed edge from node A to node B)"<<endl;
//cout<<"NOTE: Selfloops and multiple edges will be removed/ignored because the counter counts directed graphlets and orbits in networks without selfloops and multiple edges. Antiparallel pairs of arcs will be taken into account when counting. "<<endl;
//cout<<"OUTPUT FILES: "<<endl;
//cout<<"1. input_file_name.signatures.txt - each line represents a DGDV of a node (for the node names refer to the second column of input_file_name.dictionary.txt). Number of lines corresponds to number of nodes in the network minus the nodes that were only involved in selfloops."<<endl;
//cout<<"2. input_file_name.dictionary.txt - the first column corresponds to the line number from the input_file_name.signatures.txt (numbered from 0) and the second column is the name of the node for that line. "<<endl;
//cout<<"3. input_file_name.graphletcounts.txt - the first column denotes the graphlet (from 0 to 39 for the 40 directed graphlets) and the second column represents the count of the corresponding graphlet in the network."<<endl;
}
