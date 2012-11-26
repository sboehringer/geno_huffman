// main.cpp : main project file.

#include "stdafx.h"


#include <functional>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <fstream>
#include <vector>
#include <string>
#include <queue>
#include <map>
#include <stdio.h>
#include "math.h"
#include <time.h>
#include <ctime>
#include <stdlib.h>
#include <string.h>

using namespace std;

class Node;

template <typename DataType, typename Frequency> 
class Huff_tree{
    
    
  public:
    
  typedef typename std::map< DataType,  Frequency>	Huff_map;
  typedef typename std::map< DataType,  Frequency>::iterator Huff_iterator;
  Huff_tree(){
    root;
    body;
  }
  ~Huff_tree(){
   // delete root;
   // delete body;
  }
  Huff_tree(const Huff_tree&);
  
  
//======================================================
//======================================================
 class Node {
    public:
    Frequency frequency;
    DataType data; 
  
   std::vector<vector<int>> quaternary_data; //small matrix with genotypes, which are coded by 0,1,2,3
    Node* left_child;
    Node* right_child;
    bool exist_child;
    std::vector<bool> encoded;
    bool is_root;
   
	Node (Node* left, Node* right){ //for internal nodes
      left_child = left;
      right_child = right;
      exist_child = true;
      this->frequency = left_child->frequency + right_child->frequency; //frequency of node is sum of left child frequency and right child frequency
      data = 0;         
      quaternary_data;
      encoded;
      is_root = false;
    }
    
    Node (Frequency f, DataType d, std::vector<vector<int>> quater){ //for leafs
      frequency = f;
      data  = d;
      left_child = NULL;
      right_child = NULL;
      exist_child = false;
      quaternary_data = quater;
      encoded;
      is_root = false;
    }
    
    Node(){
      frequency = 0;
	  exist_child = false;
	  is_root = false;
    }
    
    ~Node(){
    
	 if(exist_child){
	 delete right_child;
	 delete left_child;
	 }
           
     }

   bool operator()(Node* a, Node* b){ //compare 2 nodes using their frequencies
      if (a->frequency < b->frequency){
	return false;
      }
	  
      else {
		  if (a->frequency == b->frequency){
			  if (!a->exist_child && !b->exist_child){
			  if (a->data < b->data){
				  return false;
			  }
			  else {
				  return true;
			  }
			  }
			  else {
				  if (a->left_child > b->left_child){
					  return true;}
				  else{
				  return false;
				  }
			  }
		  }
		  else{
			  return true;
	     }
      }
   }
    
  
  
  //==========================================
  void fill(std::vector<bool> prefix, std::map<DataType, std::vector<bool> >& code){ //fill the tree; for every left child edge will be 0, for right child - 1. The way from root to leaf is the code.
           
    if (exist_child){
      prefix.push_back(0);
      left_child->fill(prefix, code);
      prefix.pop_back();
      prefix.push_back(1);
      right_child->fill(prefix, code);
	
    }
    else{
      if (is_root){
	
	encoded.push_back(0);
      }
      else{
    encoded = prefix;
      }
    }
        
  }
  
  //============================================
  void print(std::map< std::vector<vector<int>>, std::vector<bool> > &code_1){ //print the leafs and their codes on the screen
    
       if (exist_child){
         left_child->print(code_1);
         right_child->print(code_1);
	   }

       else{
         
	 code_1[quaternary_data] = encoded;
  
	 
	 std::vector<bool>::iterator it;
	 std::vector<vector<int>>::iterator it2;
	 std::vector<int>::iterator it3;
         for ( it2=quaternary_data.begin() ; it2 != quaternary_data.end(); it2++ ){
	    
		for( int l = 0; l < (int)(*it2).size(); l++){
			std::cout << (*it2)[l]; 
		}
		cout << "__";
         } 
	std::cout << "\t";
	 
         for ( it=code_1[quaternary_data].begin() ; it !=code_1[quaternary_data].end(); it++ ){
	   std::cout << *it; 
         } 
     std::cout << "\t"<< frequency << "\n";
	  
	 std::cout <<"\n"; 


      }
    }
    
   
}; // end of class Node
//======================================================
//=======================================================
    
  
 
 Node* root;
 std::vector<Node*> body;
 typename std::vector<Node*>::iterator body_iterator;    
 
  //=====================================================================


void make_leafs_for_2_snp(map < vector<vector <int>> , double> &prob_matrix){ //for every element from prob_matrix create leaf
  body.clear();
  map < vector<vector <int>> , double>::iterator it;
  int d = 0;
  for (it=prob_matrix.begin(); it != prob_matrix.end(); it++ ) {

	 if ((*it).second != 0.0){
	  Node* dataNode = new Node((*it).second, d, (*it).first); 
     body.push_back(dataNode);
	 d++;
	 }
  }
      
}

//=======================================================================
Node* construct_tree(){ //join all leafs into tree
	priority_queue<Node*, vector<Node*>, Node> pqueue;
        
  int size1 = (int)body.size();
  for( int i = 0; i < size1; i++){
  pqueue.push(body[i]);
  }

  while (!pqueue.empty()){
    
       Node* top = pqueue.top();
       pqueue.pop();
       if (pqueue.empty()){
          root = top;
	  root->is_root = true;
       }
       else {
        Node* top2 = pqueue.top();
        pqueue.pop();
        pqueue.push(new Node(top, top2));
       }
   }
   std::cout << "Leafs:  " << body.size() << std::endl;
   int size = (int)pqueue.size();
   
       return root; 
}
//=========================================================================== 
void bite_cost (vector<vector <double>> &costs_cor, vector <vector <int>> &set_cor_snp, vector <long double> &individ_costs_add, const int Number_of_ind, const int blocks) {//calculate bit costs for correlated SNPs.
	//Bit cost characterizes contribution of given genotype into weight of whole block. Weight of block is the number of bits required for encoding this block. 
    vector <double> weight_old_row (4, 1.0); //
	int Number_of_corr_snp = (int)set_cor_snp.size();
	vector< vector <double>> weight_old (Number_of_corr_snp, weight_old_row);
   	
    const int Niteration = 30; 
    
    double size = (double) body.size();
    for (int i = 0; i<Niteration; i++) { //iterative algorithm; number of iteration equals 30; usually it is enough for accurate calculation.
		vector <double> weight_row (4, 0.0);
		vector< vector <double>> weight (Number_of_corr_snp, weight_row);
  
		for (body_iterator=body.begin() ; body_iterator != body.end(); body_iterator++ ){

			double block_weight = 0;
			for ( int k = 0; k < Number_of_corr_snp ; k++){ //calculate the weight of block, based on previous estimations
				for (int l = 0; l < blocks; l++){ 
					block_weight +=  weight_old[k][(*body_iterator)->quaternary_data[k][l]];
				}
			}

			  	  
			for (int k   = 0; k < Number_of_corr_snp ; k++){// calculate new weights, based on previous estimations 
				for (int l = 0; l < blocks; l++){
					weight[k][(*body_iterator)->quaternary_data[k][l]] +=  (1/size)*(weight_old[k][(*body_iterator)->quaternary_data[k][l]]/block_weight)*(*body_iterator)->encoded.size();
				} 
			}
      
		}
  
		weight_old = weight;
	  
    }
	costs_cor = weight_old;
   
	for (int k = 0; k  < Number_of_corr_snp; k++){// increase individual costs
		for (int l = 0; l < Number_of_ind; l++) {
			individ_costs_add[l] += (long double)costs_cor[k][set_cor_snp[k][l]];
		}
	}
 
}
    

void string_encode_for_corr (vector <vector <int> > &set_cor_snp, vector <int> &cor_snp_number,  std::map< std::vector <vector<int>>, std::vector<bool> > &code_1, int &cur_length, const int Number_of_ind, const int blocks, vector <bool> &cur_string) {
    // encode snp string using calculated code 
	std::vector<bool>::iterator it2;
    std::vector <vector<int>> part;
   vector< vector <int>> :: iterator it_snps;
    
 	cout << "\n";
    cout << "  length = " << Number_of_ind*2*(int)cor_snp_number.size()<<  "\n";
    
      for ( int i = 0; i <Number_of_ind; i+=blocks){
          for (it_snps = set_cor_snp.begin(); it_snps != set_cor_snp.end(); it_snps++){
	   vector <int> r((*it_snps).begin()+i,(*it_snps).begin()+i+blocks);
		   part.push_back(r);
	   }
		
       
      for ( it2=code_1[part].begin() ; it2 != code_1[part].end(); it2++ ){
 	   cur_string.push_back(*it2);
	   cur_length++;
	}
	 
     part.clear();
              	          
    }
    
    cout << "  length = " << cur_length<<  "\n";
    
  }

  
};// end of class Huff_tree;


class Genotypes{
public:
  char* buffer;
  unsigned char start[3];
  int Nid, Nsnps;
  
 Genotypes(char *filename, int Nid1, int Nsnps1) : Nid(Nid1), Nsnps(Nsnps1) {//read input file
    char *file = filename;
   FILE *in = fopen(file, "r");
  if (!in)
    cout << "Couln't open input file: " << "\t" << file <<"\n";
  
  if (fread(start, 1, 3, in)!=3)
    cout << "Failed to read first 3 bytes" << "\n";
  if (start[0]!='\x6C' || start[1]!='\x1B')
    cout << "Input file does not appear to be a .bed file (%X, %X)" << "\n";
    fseek (in, 0, SEEK_END);
    long filesize = ftell (in);
     buffer  = (char*)malloc(filesize);
     fseek (in, 0, SEEK_SET);
    fread (buffer, filesize, 1, in); 
  }
 //========================================================== 
  int genotype(const int snp, const int id) {//return genotype (0, 1, 2 or 3) from given position (number of snp string, number of individual)
    const unsigned char recode[4] = {'\x01', '\x00', '\x02', '\x03'};
    int bites_by_snp = ((Nid+3)/4);
    int start = bites_by_snp*snp +3;
    int ind = id/4;
    
    int gt = (buffer[start+ind] >> (2*(id%4))) & 3;
    return recode[gt];
  }
  //==================================================
  
  void vector_for_snp (int snp1, std::vector<int> &row1) {//fill snp string
    
    for (int i = 0; i < Nid; i++) {
     row1.push_back(genotype(snp1 ,i));
    // cout<< genotype(snp1 ,i);
      
    }
    
  }
  
//===========================================================
  
   void string_frequency (vector <vector <double>> &symbol_probability, vector <vector<int>> &set_cor_snp, const int Number_of_ind) {//estimate frequency of each genotype in given snp
   
	for (int a = 0; a < (int)set_cor_snp.size(); a++){
		  for  (int b = 0; b < Number_of_ind; b++){
    symbol_probability[a][set_cor_snp[a][b]]+=1.0/Number_of_ind;
		  }
	}
	
  }
  
  
  //===================================================================================
  
  bool correlation (const int snp1, const int snp2,  std::vector<int> &row1, std::vector<int> &row2, const int Number_of_ind) {// calculate correlation betweeen two snps
  //if R^2 > 0,81 return true
	  const double Threshold = 0.81;
    double mean1 = 0.0;
    double mean2 = 0.0;
    int length = 0;
    for (int i=0; i<Number_of_ind; i++) {
      if (row1[i]*row2[i]) {
       mean1+=row1[i];
       mean2+=row2[i];
       length++;
    }
    }
    mean1/=length;
    mean2/=length;

    double covariance = 0.0;
    double correlation = 0.0;
    double standdiv1 = 0.0;
    double standdiv2 = 0.0;
    
    for (int i=0; i<Number_of_ind; i++) {//calculate standart deviations
      if (row1[i]*row2[i]) {
	covariance += (row1[i] - mean1)*(row2[i] - mean2);
	standdiv1+=(row1[i] - mean1)*(row1[i] - mean1);
        standdiv2+=(row2[i] - mean2)*(row2[i] - mean2);        
	
      }

     }
    
   correlation = covariance/sqrt(standdiv1*standdiv2);
   covariance/=(Number_of_ind -1);
    if (correlation*correlation > Threshold) {
    cout << snp1 << "\t" << snp2 << "\t" << correlation << "\n";
    return true;
    }
    
    else{ return false;
    }
   
  cout << endl;
  cout << endl;
  
  
    
}


  void calculate_prob(map < vector < vector <int> >, double> &prob_matrix, vector< vector <double>> &symbol_probability, const int max_a, const int max_b, vector< vector <int>> &key, vector <int> &row,  double matrix_freq, int a, int b, int s ){
	  //in prob_matrix set of genotypes corresponds with probability, that these ganotypes appear together. 
	 //method calculate_prob help to list all variants of sets and calculate their probabilities.
	   if (symbol_probability[a][s] != 0.0){
			 row.push_back(s);
			 matrix_freq*=symbol_probability[a][s];
			b++;
			 if (b == max_b){
				 key.push_back(row);
				 a++;
				if (a == max_a){
					prob_matrix[key]=matrix_freq;
					matrix_freq /= symbol_probability[a-1][row.back()];
					key.pop_back();
					a--;
					row.pop_back();
					b = max_b - 1;
								
				}
				else {
					
					row.clear();
					b = 0;
					
				calculate_prob(prob_matrix, symbol_probability, max_a, max_b, key, row, matrix_freq, a, b, 0);
				vector <int> row_b = key.back();
				for (int bbb = 0; bbb < max_b; bbb++){
					matrix_freq /= symbol_probability[a][row_b[bbb]];
				}
				a--;
				key.pop_back();
				b = max_b - 1;
				row = row_b;
				row.pop_back();

				}
			
			 }
			 else {
				 
				 calculate_prob(prob_matrix, symbol_probability, max_a, max_b, key, row, matrix_freq, a, b, 0);
				 b--;
				 int h = row.back();
				 matrix_freq /= symbol_probability[a][h];
				 row.pop_back();
			 }
			 
		 }
		 
	if (s < 3){	
 s++;
 calculate_prob(prob_matrix, symbol_probability, max_a, max_b, key, row, matrix_freq, a, b, s);
	}
	
 }

  bool check_allele (vector <int> &row1, double maf_max, double maf_min, const int Number_of_ind) {// calculate MAF
  vector <int>::iterator count_1;
  int minor_allele = 0;
  for (count_1 = row1.begin(); count_1 != row1.end(); count_1++){
	  if ((*count_1) == 1){
		  minor_allele += 2 ;
	  }
	  else{
		  if((*count_1) == 3){
			  minor_allele += 1;
		  }
	  }

  }
  double MAF = (double)minor_allele/(Number_of_ind*2);
  if ((MAF > maf_min) && MAF < maf_max){
	  return true;
  }
  else{
	  return false;
  }

  }

  void blocking_across_snps (const int window_size, const int Number_of_ind, const int Number_of_snps, Huff_tree <int, double> &tree, Huff_tree <int, double> &min_tree, 
	  double maf_max, double maf_min, ofstream &result, ofstream &ind_result ){ //
	
	  int position_snp = 0; //number of current snp string in .bim file 
	  vector <int> row1; //current snp vector, each element corresponds to individual genotype: 0 - missing genotype, 1 - minor homozygote, 2 - major homozygote, 3 - heterozygote 
	  vector <vector <int>> window; //set of neighboured snps.
	  vector <int> snp_number; //numbers of snps, which are included in window
	  
	  vector <vector <int>> set_cor_snp; // set of correlated snps
	  vector <int> cor_snp_number; // numbers of snps, which are included in set_cor_snp
	  vector < vector <int>> :: iterator it_seq; //
	  vector <int>:: iterator it_number; //
	  std::vector<long double> individ_costs(Number_of_ind, 0.0); //
	  std::vector<long double> zero(Number_of_ind, 0.0); //
	  
	  while ((position_snp < Number_of_ind ) || ((int) window.size() > 0)){ //
		
		  while (((int) window.size() < window_size) && (position_snp < Number_of_ind))  { //fill window
			row1.clear ();
			
			vector_for_snp (position_snp, row1);
			if (check_allele(row1, maf_max, maf_min, Number_of_ind)){
			window.push_back(row1);
			snp_number.push_back(position_snp);
			}
			position_snp++;
		}
		 
		std::vector<long double> individ_costs_add(Number_of_ind, 0.0); //
		std::vector<long double> min_individ_costs_add(Number_of_ind, 0.0); //

		int current_snp = 0; //number of snp in window
		
		set_cor_snp.push_back(window[0]); //
		cor_snp_number.push_back(snp_number[0]);
		for (int i=current_snp+1; i<(int)window.size(); i++) { //fill set_cor_snp by correlated snps
		  

	 	    if( correlation (current_snp, i, window[current_snp], window[i], Number_of_ind)){
				set_cor_snp.push_back(window[i]);
				cor_snp_number.push_back(snp_number[i]);
				current_snp = i;
				it_seq = window.begin() + i;
				it_number = snp_number.begin() + i;
				window.erase(it_seq);
				snp_number.erase(it_number);
		
			}
		}
		window.erase(window.begin());
		snp_number.erase(snp_number.begin());

	
		map < vector < vector <int>> , double> prob_matrix;

		int optimal_blocks;
		int min_length = Number_of_ind*4*(int)set_cor_snp.size();
		int cur_length = 0;
		int max_block_size = 5 - (int)set_cor_snp.size()*(int)set_cor_snp.size();
		if (max_block_size < 1) {
			max_block_size = 1;
		
		}
		const int min_block_size = 1;
		vector <bool> cur_string;
		vector <bool> min_string;

	//	int freq_symb;
		vector <double> r(4,0.0);
		vector <vector <double>> symbol_probability (set_cor_snp.size(),r);
		string_frequency (symbol_probability, set_cor_snp, Number_of_ind);
        std::vector<bool> prefix;
		std::map<int, std::vector<bool> > code;
		std::map< std::vector <vector<int>>, std::vector<bool> > code_1;

		cout << "Correlated snps:  ";
    std::vector<int>::iterator it;
		for ( it=cor_snp_number.begin() ; it != cor_snp_number.end(); it++ ){
      cout << *it << "\t";
    }
	cout << "\n";

		for (int blocks = min_block_size; blocks <= max_block_size; blocks++){
			
			cout <<"blocks: " << "\t" << blocks << "\n";
			tree.body.clear();
			vector <int> row;
			vector <vector <int>> key;
			calculate_prob(prob_matrix, symbol_probability, (int)set_cor_snp.size(), blocks, key, row, 1.0, 0, 0, 0);
			
			tree.make_leafs_for_2_snp(prob_matrix); 
			tree.construct_tree()->fill(prefix,code);
	
			tree.root->print(code_1);
			tree.string_encode_for_corr(set_cor_snp, cor_snp_number, code_1, cur_length, Number_of_ind, blocks, cur_string);
			vector <vector <double>> costs_cor;
			tree.bite_cost(costs_cor, set_cor_snp, individ_costs_add, Number_of_ind, blocks);

      
			if (cur_length < min_length)  {
				min_length = cur_length;
				optimal_blocks = blocks;
				 cur_length = 0;
				min_string = cur_string;
				min_individ_costs_add = individ_costs_add;
			
     
			}
	
   cur_length = 0;
   cur_string.clear();
   individ_costs_add=zero;

delete tree.root;
   tree.body.clear();
    prefix.clear();
   code.clear();
   code_1.clear();
   prob_matrix.clear();
	 

		}
	
	 std::vector<vector<int>>::iterator it2;
	 std::vector<int>::iterator it3;
	 cout << "Source string:\n";
         for ( it2=set_cor_snp.begin() ; it2 != set_cor_snp.end(); it2++ ){
	    
		for( it3 = (*it2).begin(); it3 != (*it2).end(); it3++){
			std::cout << *it3; 
		}
		cout << "\n";
         } 
		
   vector <bool>::iterator it1;
   cout << "Min_string:    ";
   for (it1=min_string.begin(); it1 != min_string.end(); it1++){
	   cout << (*it1);
   }
   cout << "\n" << "Optimal_blocks:  "<< optimal_blocks << "\n" << "Min_length:    " << min_length << "\n";
   for (int l = 0; l < Number_of_ind; l++){
	   individ_costs[l]+=min_individ_costs_add[l];
   }
  
   
   
   result << cor_snp_number[0];
   for (it3 = cor_snp_number.begin()+1; it3 != cor_snp_number.end();  it3++){
	   result <<", " <<  *it3;
   }

   result << "\t";

   std::vector<vector<double>>::iterator it4;
	 std::vector<double>::iterator it5;
	 
         for ( it4=symbol_probability.begin() ; it4 !=symbol_probability.end(); it4++ ){
	    
		for( it5 = (*it4).begin(); it5 != (*it4).end(); it5++){
			result << *it5 << " "; 
		}
		result << ";";
         } 
	result <<"\t"<< optimal_blocks <<"\t" << min_length << "\n";
   symbol_probability.clear();
   set_cor_snp.clear();
   cor_snp_number.clear();
   min_individ_costs_add.clear();
   
   
	  // min_tree.root.clear();
   cout << "=========================\n";
	  }

	  for (int h = 0; h < Number_of_ind; h++){
		  ind_result << h << "\t" << individ_costs[h] << "\n";
	  }
  }

  

};// end of class Genotypes


//==================================================================================================


 
 
int main(int argc, char **argv) {
  
	 char* myFile;
		char* maf_max; 
		char* maf_min;
		char* snp;
		char* ind;
   if (argc < 11) { // Check the value of argc. If not enough parameters have been passed, inform user and exit.
        std::cout << "Usage is -in <input file name> -maf_min <> -maf_max <> -snp <number of snps> -ind <number of individuals>\n"; 
        std::cin.get();
        exit(0);
    } else { 
       
        std::cout << argv[0];
        for (int i = 1; i < argc; i+=2) { 
            if (i + 1 != argc){ 
               
				if (string(argv[i]) =="-in") {
                    
               	myFile = argv[i + 1];
                } else if (string(argv[i]) == "-maf_min") {
                    maf_min = argv[i + 1];
                } else if (string(argv[i]) == "-maf_max") {
                    maf_max = argv[i + 1];
				} else if (string(argv[i]) == "-snp") {
                    snp = argv[i + 1];
				} else if (string(argv[i]) == "-ind") {
                    ind = argv[i + 1];
                } else {
                    std::cout << " Not enough or invalid arguments, please try again.\n";
                    
                    exit(0);
            }
            std::cout <<" "<< argv[i] << " ";
        }
		}
   }
//const int Number_of_snps = 2239393;
 // const int Number_of_ind = 60;
 int Number_of_snps = (int)atof(snp);
 int Number_of_ind = (int)atof(ind);

  int window_size = 50;
  if (Number_of_snps >50 ){
	  window_size =50;
  }
  else {
	  window_size = Number_of_snps;
  }
  
  
  

  long double whole_length = 0.0; 
  time_t start,end;
  start = time(NULL);
  
   ofstream result("result_test.txt");  
   //make a head of output file
   result << "Correlated snps\tfrequencies\tblocks\tlength\n";
   ofstream ind_result("ind_result.txt");
   std::cout << "Hello, world!" << std::endl;
   
   
   std::vector<int> row1; //will be string from file
   std::vector<int> row2; //
    
  
   
  Genotypes* genotype = new Genotypes (myFile, Number_of_ind, Number_of_snps);
  Huff_tree <int, double>* tree = new Huff_tree <int, double>;

   Huff_tree <int, double>* min_tree = new  Huff_tree <int, double>; 

   genotype->blocking_across_snps (window_size, Number_of_ind, Number_of_snps, *tree, *min_tree, atof(maf_max), atof(maf_min), result, ind_result);

  end = time(NULL);
 
  cout <<  difftime(end,start) << "\n";
  
  return 0;
    
}