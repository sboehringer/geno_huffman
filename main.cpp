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
      this->left_child = 0;
      this->right_child = 0;
      exist_child = false;
      quaternary_data = quater;
      encoded;
      is_root = false;
    }
    
    Node(){
      frequency = 0;            
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
      else {return true;}
    }
    
  
  
  //==========================================
  void fill(std::vector<bool> prefix, std::map<DataType, std::vector<bool> >& code){
       
      
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
  void print(std::map< std::vector<int>, std::vector<bool> > &code_1){
    
       if (exist_child){
         left_child->print(code_1);
         right_child->print(code_1);
	 
       }

       else{
         
	 code_1[quaternary_data] = encoded;
  

#if 0
	 std::vector<bool>::iterator it;
	 std::vector<int>::iterator it2;
         for ( it2=quaternary_data.begin() ; it2 != quaternary_data.end(); it2++ ){
	    std::cout << *it2; 
         } 
	std::cout << "\t";
	 
         for ( it=code_1[quaternary_data].begin() ; it !=code_1[quaternary_data].end(); it++ ){
	   std::cout << *it; 
         } 
     // std::cout << "\t"<< frequency << "\n";
	  

         
std::cout <<"\n"; 
	 
#endif
      }
    }
    
 void clear_tree(){
 // cout << "aaa\n";
  if(right_child->exist_child){
    right_child->clear_tree();
    left_child->clear_tree();
    
  }
  else{
    //cout << "aaa\n";
    delete right_child;
    delete left_child;
  }
  
}
   
}; // end of class Node
//======================================================
//=======================================================
    
 
 
 
 Node* root;
 std::vector<Node*> body;
 typename std::vector<Node*>::iterator body_iterator;    
 
  //=====================================================================
  void make_leafs(Huff_tree::Huff_map letter_probability, int blocks) {
       int length = 0;
    //Huff_tree::Huff_map &leafs = * new Huff_tree::Huff_map();
    int max = (int)pow(4.0, (double)blocks);
    
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

//=======================================================================
Node* construct_tree(){
    priority_queue<Node*, vector<Node*>, Node> pqueue;
        
     
    
     
  for (body_iterator=body.begin(); body_iterator != body.end(); body_iterator++ ){
      pqueue.push(*body_iterator);
  }
    
 
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
   std::cout << body.size() << std::endl;
   int size = pqueue.size();
       return root; 
}
//=========================================================================== 
  
  
  
  void calculate_expectation (int blocks){
    
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
 
 std::vector<int> quaternary_convertion(int i, int blocks) {
  
    std::vector<int> quaternary;
    
      int copy_i = i;
      while (quaternary.size() < blocks){
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
#if 0
  void bite_cost (vector <double> &costs ) {
    double weight_old [4] = {1.0, 1.0, 1.0, 1.0};
    double weight[4];
    for (int i = 0; i<21; i++) {
      weight[0] = 0.0;
      weight[1] = 0.0;
      weight[2] = 0.0;
      weight[3] = 0.0;
 //     = {0.0, 0.0, 0.0, 0.0};
      typename std::vector<Node*>::iterator iter;
      for ( iter=body.begin() ; iter != body.end(); iter++ ){
	//if (!(*iter)->exist_child){
	  double block_weight = 0;
	  std::vector<int>::iterator it2;
	  for ( it2=(*iter)->quaternary_data.begin() ; it2 != (*iter)->quaternary_data.end(); it2++ ){
	      block_weight +=  weight_old[(*it2)];
	  }
	  for ( it2=(*iter)->quaternary_data.begin() ; it2 != (*iter)->quaternary_data.end(); it2++ ){
	    weight[(*it2)] += (!block_weight)? 1e-6: (*iter)->frequency*(weight_old[(*it2)]/block_weight)*(*iter)->encoded.size();
	  }
      //  }
      }
  // std::cout << weight[0] << "\t" << "\t" << weight[1] << "\t" <<"\t" <<  weight[2] << "\t" << "\t" << weight[3] << "\n";
      weight_old[0]= weight [0];
      weight_old[1]= weight [1];
      weight_old[2]= weight [2];
      weight_old[3]= weight [3];
           
    }
   
    costs.push_back(weight [0]);
    costs.push_back(weight [1]);
    costs.push_back(weight [2]);
    costs.push_back(weight [3]);
  
           
  // result << weight[0] << "\t" << weight[1] << "\t" << weight[2] << "\t" << weight[3] << "\n";
   //std::cout << "=========================================" << "\n";
  }
#endif 
  //=============================================================================
  void bite_cost_2 (int blocks, vector <double> &costs) {
    double weight_old [4] = {1.0, 1.0, 1.0, 1.0};
    double weight [4];
    
    double size  =(double) body.size();
    for (int i = 0; i<20; i++) {
      weight[0] = 0.0;
      weight[1] = 0.0;
      weight[2] = 0.0;
      weight[3] = 0.0;
      
      
      for (body_iterator=body.begin() ; body_iterator != body.end(); body_iterator++ ){
	//if (!(*body_iterator)->exist_child){
	  double block_weight = 0;
	  std::vector<int>::iterator it2;
	  for ( it2=(*body_iterator)->quaternary_data.begin() ; it2 != (*body_iterator)->quaternary_data.end(); it2++ ){
	      block_weight +=  weight_old[(*it2)];
	  }
	  
	//  weight[(*body_iterator)->quaternary_data[0]]+= (!block_weight)? 1e-12: (*body_iterator)->frequency*(weight_old[(*body_iterator)->quaternary_data[0]]*(1/block_weight))*(*body_iterator)->encoded.size()*blocks;
	   weight[(*body_iterator)->quaternary_data[0]]+= (1/size)*(weight_old[(*body_iterator)->quaternary_data[0]]/block_weight)*(*body_iterator)->encoded.size()*blocks;
	     //cout << weight_old[(*body_iterator)->quaternary_data[0]]/block_weight << "\t"<< (max/size) <<"\n";
	  //for ( it2=(*body_iterator)->quaternary_data.begin() ; it2 != (*body_iterator)->quaternary_data.end(); it2++ ){
	   //   weight[(*it2)] +=  (*body_iterator)->frequency*(weight_old[(*it2)]/block_weight)*(*body_iterator)->encoded.size();
	 // }
      //  }
      }
     // std::cout << weight[0] << "\t" << "\t" << weight[1] << "\t" <<"\t" <<  weight[2] << "\t" << "\t" << weight[3] << "\n";
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
   void string_encode (std::vector<int> string_snp, int blocks, std::map< std::vector<int>, std::vector<bool> > &code_1, int &cur_length) {
    std::vector<int>::iterator it;
    std::vector<bool>::iterator it2;
    std::vector<int> part;
    int i = 0;
    
    cout << "Source string:  ";
    for ( it=string_snp.begin() ; it != string_snp.end(); it++ ){
      cout << *it;
    }
    cout << "  length = " << string_snp.size()<<  "\n";
    cout << "Encoded string:  " ;
     for ( it=string_snp.begin() ; it != string_snp.end(); it++ ){
       //cout << *it;
       part.push_back(*it);
       i++;
     
       if (i == blocks) {
	// cout << code_1[part].size()<< "\t";
	 for ( it2=code_1[part].begin() ; it2 != code_1[part].end(); it2++ ){
      	 cout << *it2;
	 cur_length++;
	 }
	 i = 0;
	 part.clear();
      }
      
     	          
    }
    
    cout << "  length = " << cur_length<<  "\n";
    
  }
#if 0 
  void clear_leafs(){
  typename std::vector<Node*>::iterator iter;
      for ( iter=body.begin() ; iter != body.end(); iter++ ){
	delete (*iter);
      }
    
  }
#endif 
  
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
  int genotype(int snp, int id) {
    const unsigned char recode[4] = {'\x01', '\x00', '\x02', '\x03'};
    int bites_by_snp = 3+((Nid+3)/4);
    int start = bites_by_snp*snp +3;
    int ind = id/4;
    int gt = (buffer[start+ind] >> (2*(id%4))) & 3;
    return recode[gt];
  }
  //==================================================
  
  void vector_for_snp (int snp1, std::vector<int> &row) {
    
    for (int i = 0; i < Nid; i++) {
     row.push_back(genotype(snp1 ,i));
    // cout<< genotype(snp1 ,i);
      
    }
    
  }
  
//===========================================================
  
   std::map<int, double> string_frequency (vector<int> &row) {
   std::vector<int>::iterator it;
   
   std::map<int, double> symbol_probability;
   // double	symbol_probability[4];
    symbol_probability[0] = 0.0;
    symbol_probability[1] = 0.0;
    symbol_probability[2] = 0.0;
    symbol_probability[3] = 0.0;
    
    double length = (double)row.size();
    cout << "!!!!" << length << "\n";
    for ( it=row.begin() ; it != row.end(); it++ ){
      symbol_probability[*it] += 1/length;
	          
    } 
    cout << "============================================"<< "\n" ;
   cout << symbol_probability[0] << "\t" <<  symbol_probability[1] << "\t" <<  symbol_probability[2] << "\t" <<  symbol_probability[3] << "\n" ;
    cout << "============================================"<< "\n" ;
   return symbol_probability;		 
  //return *new std::map<int, double>;
     
  }
  
};// end of class Genotypes

  
int main(int argc, char **argv) {
  long double whole_length = 0.0; 
  time_t start,end;
   start = time(NULL);
   ofstream result("result1.txt");  
  result << "Snp" << "\t"<< "blocks" << "\t"<< "length" <<"\t"<< "cost[0]" << "\t"<< "cost[1]" << "\t"<< "cost[2]" << "\t"<< "cost[3]" << "\n";
  std::cout << "Hello, world!" << std::endl;
    
   // std::map<int, double> letter_probability;
   
  //  letter_probability[0] = 0.4;
   // letter_probability[1] = 0.3;
    //letter_probability[2] = 0.2;
    //letter_probability[3] = 0.1;
   std::vector<bool> prefix;
   std::map<int, std::vector<bool> > code;
   std::map< std::vector<int>, std::vector<bool> > code_1;
   std::vector<int> row;
   std::vector<long double> individ_costs;
  for(int i =0; i < 60; i++){
  individ_costs.push_back(0.0);
    
  }
  int blocks;
  
  Genotypes* genotype = new Genotypes ("hapmap-ceu.bed", 60, 2239392);
  Huff_tree <int, double>* tree = new Huff_tree <int, double>;
  
  for (int snp_number = 0; snp_number < 2239392; snp_number++){
  cout << "snp  " << snp_number << "\n";
  int optimal_blocks;
  int min_length = 400;
  int cur_length = 0;
  vector<double> costs;
  
   genotype->vector_for_snp(snp_number, row);
    
 for (blocks = 2; blocks < 7; blocks++){
    cout <<"blocks: " << "\t" << blocks << "\n";
    tree->body.clear();
    //tree->root = NULL;

    tree->make_leafs((genotype->string_frequency(row)), blocks); 
    tree->construct_tree()->fill(prefix,code);
    tree->root->print(code_1);
    tree->string_encode(row, blocks, code_1, cur_length);
  
    
    if (cur_length < min_length)  {
    min_length = cur_length;
    optimal_blocks = blocks;
    costs.clear();
    tree->bite_cost_2(blocks, costs);
     
    }
  else {
    result << snp_number << "\t" << optimal_blocks  << "\t" << min_length  ;
    vector<double>::iterator it;
    for (it = costs.begin()  ; it != costs.end(); it++ ){
      
      result << "\t" << (*it);
    
    }
    costs.clear();
    result << "\n";
      
    
      //<< costs[0] <<"\t" <<costs[1] <<"\t"<< costs[2] <<"\t"<< costs[3] <<"\n";
    tree->calculate_expectation(blocks);
    cur_length =0;

   prefix.clear();
   code.clear();
   code_1.clear();
   delete tree->root;
    
    break;
    
  }

  if (blocks == 6){
    result << snp_number << "\t" << optimal_blocks  << "\t" << min_length  ;
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
   delete tree->root;
 // tree->root->clear_tree();
  
  }
  
  whole_length+=min_length;
 // for(int i =0; i < 60; i++){
 // individ_costs[i]+=costs[row[i]];
    
 // }
 
 
  
   row.clear();
  }
  
  ofstream result_ind("result_ind.txt"); 
  for(int i =0; i < 60; i++){
  result_ind << i<< "\t" << individ_costs[i] <<"\n";
    
  }
  end = time(NULL);
 cout << "Whole_length:  " << whole_length << "\n";
  cout <<  difftime(end,start) << "\n";
  
  return 0;
    
}
