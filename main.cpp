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
  
  
//======================================================
//======================================================
 class Node {
    public:
    Frequency frequency;
    DataType data; //symbol of alphabet
  
   std::vector<int> quaternary_data;
  // int code_length;
    Node* left_child;
    Node* right_child;
    bool exist_child;
    std::vector<bool> encoded;
    bool is_root;
   
	Node (Node* left, Node* right){
      left_child = left;
      right_child = right;
      exist_child = true;
      this->frequency = left_child->frequency + right_child->frequency;
      data = 0;         
      quaternary_data;
      encoded;
      is_root = false;
    }
    
    Node (Frequency f, DataType d, std::vector<int> quater){
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
    // cout << data << "\n";
      if(exist_child){
	 delete right_child;
	 delete left_child;
	
       }
       
       
     }

   bool operator()(Node* a, Node* b){
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
  void fill(std::vector<bool> prefix, std::map<DataType, std::vector<bool> >& code){
        
      
    if (exist_child){
    //  std::cout << "Hello, world!" << std::endl;  
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
  void print(std::map< std::vector<int>, std::vector<bool> > &code_1){
    
       if (exist_child){
         left_child->print(code_1);
         right_child->print(code_1);
	 
       }

       else{
         
	 code_1[quaternary_data] = encoded;
  


	 std::vector<bool>::iterator it;
	 std::vector<int>::iterator it2;
         for ( it2=quaternary_data.begin() ; it2 != quaternary_data.end(); it2++ ){
	    std::cout << *it2; 
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
  void make_leafs(typename Huff_tree::Huff_map letter_probability, int blocks, int &freq_symb) {
      
	  int length = 0;
    //Huff_tree::Huff_map &leafs = * new Huff_tree::Huff_map();
   const int max = (int)pow(4.0, (double)blocks);
   double max_freq = 0;
   for (int k =0; k < 4; k++){
     if (letter_probability[k] > max_freq){
       max_freq = letter_probability[k];
       freq_symb = k;
     }
   }
    
    Frequency freq;
    for (int i = 0 ; i < max; i++){
      bool zero_prob = false;
      length = 0;
      int copy_i = i;
      freq = 1;
      while (length < blocks){
	if (copy_i > 0){
	  if (letter_probability[(copy_i % 4)]== 0){
	    zero_prob = true;
	    break;
	  }
	    
	    freq *= letter_probability[(copy_i % 4)];
	    
	    copy_i /= 4;
         }	
	else{
	  if (letter_probability[0] == 0){
	    zero_prob = true;
	    break;
	  }
	 freq *= letter_probability[0];
   
	}
	length++;
     }
     if (!zero_prob){
    Node* dataNode = new Node(freq, i, quaternary_convertion(i , blocks)); 
     body.push_back(dataNode);
     
     }
    }   
   
}  


void make_leafs_for_2_snp(map < vector <int> , double> &prob_matrix){
  body.clear();
  map < vector <int> , double>::iterator it;

  for (it=prob_matrix.begin(); it != prob_matrix.end(); it++ ) {
	  vector <int>::iterator it_int;
	  vector <int> vector1 = (*it).first;
	  int sum = 0;
	  int power = (int) (*it).first.size();
	  for (it_int=vector1.begin(); it_int != vector1.end(); it_int++ ){
		  sum+=(int)pow (4.0, power)*(*it_int);
		  power--;
	  }
	  Node* dataNode = new Node((*it).second, sum, (*it).first); 
     body.push_back(dataNode);
  }
      
}

//=======================================================================
Node* construct_tree(){
    priority_queue<Node*, vector<Node*>, Node> pqueue;
        
   //  cout << "in construct tree()" << "\n";
    //cout << body.size() << "\n";
  int size1 = (int)body.size();
  for( int i = 0; i < size1; i++){
  pqueue.push(body[i]);
  }

//  for (body_iterator=body.begin(); body_iterator != body.end(); body_iterator++ ){
  //    pqueue.push(*body_iterator);
//	  std::cout << "Hello, world!111" << std::endl; 
     // cout<< (*body_iterator)->data << "\n";
 // }
    
 //std::cout << "Hello, world!111" << std::endl; 
  while (!pqueue.empty()){
    
       Node* top = pqueue.top();
       pqueue.pop();
       if (pqueue.empty()){
          root = top;
	  root->is_root = true;
	//  cout << "find_root" << "\n";
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
  
  
  
  void calculate_expectation (const int blocks){
    
    double expectation = 0;   
  for ( body_iterator=body.begin() ; body_iterator != body.end(); body_iterator++ ){
   // if (!(*body_iterator)->exist_child){
    expectation += (*body_iterator)->frequency*(*body_iterator)->encoded.size();
  //   }
  }
   std::cout << "=========================================" << "\n";
    cout << "Expectation is" << "\t" << expectation << "\n";
    cout << "before Huffman" << "\t" << blocks*2 << "\n";
    cout << (blocks*2)/expectation << "\n";
    std::cout << "=========================================" << "\n";
  }
  
//========================================================================== 
 
 std::vector<int> quaternary_convertion(const int i, const int blocks) {
  
    std::vector<int> quaternary;
    
      int copy_i = i;
      while (blocks > (int)quaternary.size()){
	if (copy_i > 0){
	  int r = copy_i % 4;
	  quaternary.push_back(r);
	  copy_i /= 4;
         }	
	 else{
	 quaternary.push_back(0);
   
         }
       }
       return quaternary;
  }

//=========================================================================

  void bite_cost_2 (const int blocks, vector <double> &costs) {
    double weight_old [4] = {1.0, 1.0, 1.0, 1.0};
    double weight [4];
    const int Niteration = 20;
    
    double size  =(double) body.size();
    for (int i = 0; i<Niteration; i++) {
      weight[0] = 0.0;
      weight[1] = 0.0;
      weight[2] = 0.0;
      weight[3] = 0.0;
      int diff_blocks = 0;
      std::map< std::vector<int>, bool >::iterator iter1;
      
      
      for (body_iterator=body.begin() ; body_iterator != body.end(); body_iterator++ ){
	//if (!(*body_iterator)->exist_child){
	//  if (is_in_string[(*body_iterator)->quaternary_data] == 1){
	  double block_weight = 0;
	  std::vector<int>::iterator it2;
	  for ( it2=(*body_iterator)->quaternary_data.begin() ; it2 != (*body_iterator)->quaternary_data.end(); it2++ ){
	      block_weight +=  weight_old[(*it2)];
	  }
	  
	   
	   
	   
	 //  weight[(*body_iterator)->quaternary_data[0]]+= (1/size)*(weight_old[(*body_iterator)->quaternary_data[0]]/block_weight)*(*body_iterator)->encoded.size()*blocks;
	     
	   
	   
	  for ( it2=(*body_iterator)->quaternary_data.begin() ; it2 != (*body_iterator)->quaternary_data.end(); it2++ ){
	      weight[(*it2)] +=  (1/size)*(weight_old[(*it2)]/block_weight)*(*body_iterator)->encoded.size();
	  }
      //  }
      
    //  }
      }
   //  std::cout << weight[0] << "\t" << "\t" << weight[1] << "\t" <<"\t" <<  weight[2] << "\t" << "\t" << weight[3] << "\n";
      weight_old[0]= weight [0];
      weight_old[1]= weight [1];
      weight_old[2]= weight [2];
      weight_old[3]= weight [3];
    }
    costs.push_back(weight [0]);
    costs.push_back(weight [1]);
    costs.push_back(weight [2]);
    costs.push_back(weight [3]);
  
  }
  
  
 //======================================================================== 
   void string_encode (std::vector<int> string_snp, const int blocks, std::map< std::vector<int>, std::vector<bool> > &code_1, int &cur_length,  int &freq_symb) {
    std::vector<int>::iterator it;
    std::vector<bool>::iterator it2;
    std::vector<int> part;
    int i = 0;
    
    cout << "Source string:  ";
    for ( it=string_snp.begin() ; it != string_snp.end(); it++ ){
     cout << *it;
    }
    cout << "  length = " << string_snp.size()<<  "\n";
  //  cout << "Encoded string:  " ;
     for ( it=string_snp.begin() ; it != string_snp.end(); it++ ){
       //cout << *it;
       part.push_back(*it);
       i++;
     
       if (i == blocks) {
	// cout << code_1[part].size()<< "\t";
	 for ( it2=code_1[part].begin() ; it2 != code_1[part].end(); it2++ ){
      	// cout << *it2;
	 cur_length++;
	 }
	 i = 0;
	 
	 part.clear();
      }
      
     	          
    }
    if (part.size() > 0){
      while (blocks > (int)part.size()){
	part.push_back(freq_symb);
      }
      for ( it2=code_1[part].begin() ; it2 != code_1[part].end(); it2++ ){
      	// cout << *it2;
	 cur_length++;
	 }
      
    }
    
    cout << "  length = " << cur_length<<  "\n";
    
  }
  
     void string_encode_for_corr (vector <vector <int> > &set_cor_snp, vector <int> &cor_snp_number,  std::map< std::vector<int>, std::vector<bool> > &code_1, int &cur_length, const int Number_of_ind) {
    std::vector<int>::iterator it;
    std::vector<bool>::iterator it2;
    std::vector<int> part;
   vector< vector <int>> :: iterator it_snps;
    
    cout << "Correlated snps:  ";
    for ( it=cor_snp_number.begin() ; it != cor_snp_number.end(); it++ ){
      cout << *it << "\t";
    }
	cout << "\n";
    cout << "  length = " << Number_of_ind*2*(int)cor_snp_number.size()<<  "\n";
    
   
    
    cout << "Encoded string:  " ;
     for ( int i = 0; i <Number_of_ind; i++){
       //cout << *it;
       for (it_snps = set_cor_snp.begin(); it_snps != set_cor_snp.end(); it_snps++){
	   part.push_back((*it_snps)[i]);
	   }
		
       
      for ( it2=code_1[part].begin() ; it2 != code_1[part].end(); it2++ ){
       cout << *it2;
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
  
 Genotypes(char *filename, int Nid1, int Nsnps1) : Nid(Nid1), Nsnps(Nsnps1) {
    char *file = filename;
   FILE *in = fopen(file, "r");
  if (!in)
    cout << "Couln't open input file: %s" << "\t" << file <<"\n";
  
  if (fread(start, 1, 3, in)!=3)
    cout << "Failed to read first 3 bytes" << "\n";
  if (start[0]!='\x6C' || start[1]!='\x1B')
    cout << "Input file does not appear to be a .bed file (%X, %X)" << "\n";
    fseek (in, 0, SEEK_END);
    long filesize = ftell (in);
     buffer  = (char*)malloc(filesize);
     fseek (in, 0, SEEK_SET);
    fread (buffer, filesize, 1, in); //fread ( void * ptr, size_t size, size_t count, FILE * stream );
  }
 //========================================================== 
  int genotype(const int snp, const int id) {
    const unsigned char recode[4] = {'\x01', '\x00', '\x02', '\x03'};
    int bites_by_snp = ((Nid+3)/4);
    int start = bites_by_snp*snp +3;
    int ind = id/4;
    
    int gt = (buffer[start+ind] >> (2*(id%4))) & 3;
    return recode[gt];
  }
  //==================================================
  
  void vector_for_snp (int snp1, std::vector<int> &row1) {
    
    for (int i = 0; i < Nid; i++) {
     row1.push_back(genotype(snp1 ,i));
    // cout<< genotype(snp1 ,i);
      
    }
    
  }
  
//===========================================================
  
   std::map<int, double> string_frequency (vector<int> &row, ofstream &result, const int blocks) {
   std::vector<int>::iterator it;
   
   std::map<int, double> symbol_probability;
   
    symbol_probability[0] = 0.0;
    symbol_probability[1] = 0.0;
    symbol_probability[2] = 0.0;
    symbol_probability[3] = 0.0;
    
    double length = (double)row.size();
   // cout << "!!!!" << length << "\n";
    for ( it=row.begin() ; it != row.end(); it++ ){
      symbol_probability[*it] += 1/length;
	          
    }
  
    cout << "============================================"<< "\n" ;
    if (blocks ==2){
   result << symbol_probability[0] << "\t" <<  symbol_probability[1] << "\t" <<  symbol_probability[2] << "\t" <<  symbol_probability[3] << "\t" ;}
    cout << "============================================"<< "\n" ;
   return symbol_probability;		 
  
     
  }
  
  
  
  //===================================================================================
  
  bool correlation (const int snp1, const int snp2,  std::vector<int> &row1, std::vector<int> &row2, const int Number_of_ind) {
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
   // cout << "length: " << length << "\n";
   //cout << mean1 << "\n";
    mean1/=length;
    mean2/=length;
    
   // cout << mean1 << "\n";
    double covariance = 0.0;
    double correlation = 0.0;
    double standdiv1 = 0.0;
    double standdiv2 = 0.0;
    
    for (int i=0; i<Number_of_ind; i++) {
      if (row1[i]*row2[i]) {
	covariance += (row1[i] - mean1)*(row2[i] - mean2);
	standdiv1+=(row1[i] - mean1)*(row1[i] - mean1);
        standdiv2+=(row2[i] - mean2)*(row2[i] - mean2);        
	
      }
//      double k = (prob_matrix.at(row1[i])).at(row2[i]);
  //    vector < double> kk = prob_matrix.at(row1[i]);
      
//      prob_matrix[row1[i]][row2[i]] += (1.0/Number_of_ind);
     //cout << (prob_matrix.at(row1[i])).at(row2[i]);
     // cout << correlation << "\t" << row1[i] << "\t" << row2[i]<< "\n";
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

  void width_search ( vector < vector<int> > &neighbours, vector < vector<int> > &result, const int snp, const int number_of_group, const int window_size, const int start) {
    vector <int> :: iterator iter;
 
    result[number_of_group][(neighbours[snp][0]-start)%window_size] = neighbours[snp][0];
    if (neighbours[snp].size() < 2){
       }
    else{
       
      for (iter=neighbours[snp].begin()+1; iter!=neighbours[snp].end(); iter++){// +++++it = 89098+1
	  width_search ( neighbours,result, (*iter-start)%window_size, number_of_group, window_size, start);
  
       }
   }
    neighbours[snp].clear();
   
}

  void calculate_prob(map < vector <int> , double> &prob_matrix, vector <vector <int>> &set_cor_snp, const int Number_of_ind){
	  
	  for (int j = 0; j < Number_of_ind; j++){
	  vector <int> key;
	  for (int k = 0; k < (int)set_cor_snp.size(); k++){
		   key.push_back(set_cor_snp[k][j]);
	   }

	   if ( prob_matrix[key]){
	   prob_matrix[key]+=1.0/Number_of_ind;
	   }
	   else {
       prob_matrix[key] = 1.0/Number_of_ind;
	   }
	  
	  }

	  }



  void blocking_across_snps (const int window_size, const int Number_of_ind, Huff_tree <int, double> &tree){
	  int position_snp = 0;
	  vector <int> row1;
	  vector <int> snp_number;
	  vector <vector <int>> window;
	  vector <vector <int>> set_cor_snp;
	  vector <int> cor_snp_number;
	  vector < vector <int>> :: iterator it_seq;
	  vector <int>:: iterator it_number;
	  while (position_snp < 300){
	  while ((int) window.size() < window_size)  {
        row1.clear ();
        vector_for_snp (position_snp, row1);
        window.push_back(row1);
		snp_number.push_back(position_snp);
		position_snp++;
	  }
	  int current_snp = 0;
	  set_cor_snp.push_back(window[0]);
	  cor_snp_number.push_back(snp_number[0]);
	  for (int i=current_snp+1; i<(int)window.size(); i++) {
		  

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

	
	map < vector <int> , double> prob_matrix;
	calculate_prob(prob_matrix, set_cor_snp, Number_of_ind);
	tree.make_leafs_for_2_snp(prob_matrix);
	 std::vector<bool> prefix;
   std::map<int, std::vector<bool> > code;
   std::map< std::vector<int>, std::vector<bool> > code_1;
   int cur_length = 0;
  
   

     //leafs already exist
    tree.construct_tree()->fill(prefix,code);
    
    tree.root->print(code_1);
    tree.string_encode_for_corr(set_cor_snp, cor_snp_number, code_1, cur_length, Number_of_ind);
	cout << "==============================\n";
   delete tree.root;
   tree.body.clear();
    prefix.clear();
   code.clear();
   code_1.clear();
   set_cor_snp.clear();
   cor_snp_number.clear();
	  }
  }


  
  void find_all_neighbours(const int Number_of_snps, const int Number_of_ind) {
      
    vector < vector<int> > neighbours;
    vector < vector<int> > result;
    vector<int> row1;
    vector<int> row2;
    const int start = 51;
    const int window_size = 50;
    vector<int> roww (window_size);
    vector<int> rowww;
    vector < vector <double> > prob_matrix;
    vector<double> row;
	 for (int k = 0; k <4; k++){
	   row.push_back(0.0);
	 }
	 for (int l = 0; l<4; l++){
	   prob_matrix.push_back(row);
	 }
   

   
   for (int i = start; i < start + window_size && i < Number_of_snps; i++) {
      row1.clear ();
      vector_for_snp (i, row1);
      cout << i << "\n";
      neighbours.push_back(rowww);
      for (int j = i; j < start + window_size && j < Number_of_snps; j++) { //++++ j=i
	row2.clear ();
        
	vector_for_snp (j, row2);
        if (correlation(i,j ,row1,row2, Number_of_ind)) {
	   
	  neighbours[i-start].push_back (j);
	  
      }
      }
    }
       
    
    vector< vector <int> >::iterator it2;
    int number_of_group= -1;
    int i = -1;
    cout << "===========================" <<"\n";
   
    
    for (it2=neighbours.begin(); it2!=neighbours.end(); it2++){
      
      i++;
      
        if ((*it2).size() != 0){
	  
	number_of_group++;
	 
         result.push_back(roww);
	 cout << "i: " << i << "  ";
	width_search (neighbours,result,i, number_of_group, window_size, start);
	cout << "\n";
      }
    }
    
    for (it2=result.begin(); it2!=result.end(); it2++){
      vector<int>::iterator iter;
      for (iter=(*it2).begin(); iter!=(*it2).end(); iter++){
	if ((*iter)!=0){
	  cout << (*iter) << "\t";
	}
	
    }
    cout << "\n";
    
  }
  }
};// end of class Genotypes


//==================================================================================================


  void construct_tree_for_blocks(Genotypes genotype, Huff_tree <int, double> tree, ofstream &result, int snp_number, vector<double> &costs, std::vector<int> &row1, const int Number_of_ind) {
    std::vector<bool> prefix;
   std::map<int, std::vector<bool> > code;
   std::map< std::vector<int>, std::vector<bool> > code_1;
   
   
  int optimal_blocks;
  int min_length = Number_of_ind*4;
  int cur_length = 0;
  const int max_block_size = 4;
  const int min_block_size = 2;
  int freq_symb;
     for (int blocks = min_block_size; blocks < max_block_size; blocks++){
    cout <<"blocks: " << "\t" << blocks << "\n";
    tree.body.clear();
    //tree->root = NULL;

    tree.make_leafs((genotype.string_frequency(row1, result,blocks)), blocks,  freq_symb); 
    tree.construct_tree()->fill(prefix,code);
	
    tree.root->print(code_1);
    tree.string_encode(row1, blocks, code_1, cur_length,  freq_symb);
 
    
    if (cur_length < min_length)  {
    min_length = cur_length;
    optimal_blocks = blocks;
    costs.clear();
    tree.bite_cost_2(blocks, costs);
     
    }
  else {
    result <<  optimal_blocks  << "\t" << min_length  ;
    vector<double>::iterator it;
    for (it = costs.begin()  ; it != costs.end(); it++ ){
      
      result << "\t" << (*it);
    
    }
    costs.clear();
    result << "\n";
      
    
      //<< costs[0] <<"\t" <<costs[1] <<"\t"<< costs[2] <<"\t"<< costs[3] <<"\n";
 // tree.calculate_expectation(blocks);
    cur_length =0;

   prefix.clear();
   code.clear();
   code_1.clear();
   delete tree.root;
    
    break;
    
  }

  if (blocks == max_block_size -1){
    result << optimal_blocks  << "\t" << min_length  ;
    vector<double>::iterator it;
    for (it = costs.begin()  ; it != costs.end(); it++ ){
         result << "\t" << (*it);
  }
    result << "\n";
      //<< costs[0] <<"\t" <<costs[1] <<"\t"<< costs[2] <<"\t"<< costs[3] <<"\n";
    //66666 tree->calculate_expectation(blocks);
  }
  

cur_length =0;
 
   prefix.clear();
   code.clear();
   code_1.clear();
   delete tree.root;
 // tree->root->clear_tree();
  
  }
  }
  
  
  //=======================================================================================
  
  void construct_tree_for_correlated_snps (Genotypes genotype, Huff_tree <int, double> tree, ofstream &result, int snp_number, vector<double> &costs, std::vector<int> &row1, std::vector<int> &row2, const int Number_of_ind) {
      
   std::vector<bool> prefix;
   std::map<int, std::vector<bool> > code;
   std::map< std::vector<int>, std::vector<bool> > code_1;
   int cur_length = 0;
  
   

     //leafs already exist
    tree.construct_tree()->fill(prefix,code);
    
    tree.root->print(code_1);
  //  tree.string_encode_for_corr(row1,row2, code_1, cur_length, Number_of_ind);
   // tree.bite_cost_2(1, costs);
    vector<double>:: iterator it;
     for (it = costs.begin()  ; it != costs.end(); it++ ){
    //  cout << "!!!!\t" << (*it);
      result << "\t" << (*it);
    
    }
    costs.clear();
   
   delete tree.root;
   tree.body.clear();
  }
  
 
 
 //=======================================================================================
 
 
int main(int argc, char **argv) {
  
  const int Number_of_snps = 2239393;
  const int Number_of_ind = 60;
  const int window_size = 50;
  
  

  long double whole_length = 0.0; 
  time_t start,end;
  start = time(NULL);
  
   ofstream result("result_test.txt");  
   //make a head of output file
   result << "Snp\tfreq[0]\tfreq[1]\tfreq[2]\tfreq[3]\tblocks\tlength\tcost[0]\tcost[1]\tcost[2]\tcost[3]\n";
   
   std::cout << "Hello, world!" << std::endl;
   
   
   std::vector<int> row1; //will be string from file
   std::vector<int> row2; //
   std::vector<long double> individ_costs(Number_of_ind, 0.0); //contained average number of bites for each individual
  
   
  // for(int i =0; i < Number_of_ind; i++){
   //   individ_costs.push_back(0.0);
   //}
  
  
   
  Genotypes* genotype = new Genotypes ("hapmap-ceu.bed", Number_of_ind, Number_of_snps);
  Huff_tree <int, double>* tree = new Huff_tree <int, double>;
  genotype->blocking_across_snps (window_size, Number_of_ind, *tree);
#if 0
  vector<bool> snp_list (Number_of_snps, true); // 0 if snp was already processed, 1 if not
//  genotype->find_all_neighbours(Number_of_snps, Number_of_ind);
  
  
  //for (int i = 0; i < Number_of_snps; i++) {
  //  snp_list.push_back(1) ;
	 
  //}
    
 
  std::cout << "Hello, world!" << std::endl; 
   for (int snp1 = 0; snp1 <5 ; snp1++){
     result << snp1 << "\t";
     
   if (snp_list[snp1]) {
     std::cout << "SNP:   " << snp1 << "\n";
     genotype->vector_for_snp (snp1, row1); //recieve data for SNP and put it in the vector<int> row1
     bool correlated = false; // = true, if we found correlating SNP
     
     vector<double> costs; //estimated bite costs of genotypes
     
       for (int snp2 = snp1+1; ((snp2 < snp1+window_size )&&( snp2 < Number_of_snps)); snp2++) {
	
	 if (snp_list[snp2]){
	
	 
	 //fill prob_matrix
	 vector<double> row (4, 0.0);
	// for (int k = 0; k <4; k++){
	//   row.push_back(0.0);
	 //}
	  vector< vector <double> > prob_matrix (4, row);
	 //for (int l = 0; l<4; l++){
	 //  prob_matrix.push_back(row);
	 //}
	 
	 
	 row2.clear();
	 genotype->vector_for_snp (snp2, row2); //take data for second snp
         
	 if(genotype-> correlation (snp1, snp2,  row1, row2, Number_of_ind)){
	//    tree->make_leafs_for_2_snp(prob_matrix);
	 
	 correlated = true;
	 
	 snp_list[snp2] = 0;
	
	 break;
          
	 }
      
       }
	 
      }
       
       if (!correlated){
		   
       construct_tree_for_blocks(*genotype, *tree, result, snp1, costs, row1, Number_of_ind);
	   
       }
       else{
	 
	 construct_tree_for_correlated_snps (*genotype, *tree, result, snp1, costs, row1, row2, Number_of_ind);
	 
       }
       snp_list[snp1] = 0;
        for(int i =0; i < Number_of_ind; i++){
        individ_costs[i]+=costs[row1[i]];
    
        }
       row1.clear();
      
    } 
    
  }
  
   
     
  ofstream result_ind("result_ind_plink.txt"); 
  for(int i =0; i < Number_of_ind; i++){
  result_ind << i<< "\t" << individ_costs[i] <<"\n";
    
  }
  cout << "Whole_length:  " << whole_length << "\n";

#endif 
  end = time(NULL);
 
  cout <<  difftime(end,start) << "\n";
  
  return 0;
    
}
