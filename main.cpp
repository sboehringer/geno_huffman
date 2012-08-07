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
	  

         
std::cout <<"\t"; 


      }
    }
    

   
}; // end of class Node
//======================================================
//=======================================================
    
 
 
 
 Node* root;
 std::vector<Node*> body;
 typename std::vector<Node*>::iterator body_iterator;    
 
  //=====================================================================
  void make_leafs(Huff_tree::Huff_map letter_probability, int blocks, int &freq_symb) {
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


void make_leafs_for_2_snp(vector< vector< double> > &prob_matrix){
  body.clear();
 // cout << "in make_leafs()" << "\n";
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      vector <int> quater;
      quater.push_back(i);
      quater.push_back(j);
     // cout << prob_matrix[i][j] << "\n";
   //   cout << (i*4)+j << "\n";
    if(prob_matrix[i][j]){
      Node* dataNode = new Node(prob_matrix[i][j], (i*4)+j, quater); 
     body.push_back(dataNode);
    }
    }
    }
  
}

//=======================================================================
Node* construct_tree(){
    priority_queue<Node*, vector<Node*>, Node> pqueue;
        
   //  cout << "in construct tree()" << "\n";
    
     
  for (body_iterator=body.begin(); body_iterator != body.end(); body_iterator++ ){
      pqueue.push(*body_iterator);
     // cout<< (*body_iterator)->data << "\n";
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
   std::cout << "Leafs:  " << body.size() << std::endl;
   int size = pqueue.size();
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
      while (part.size() < blocks){
	part.push_back(freq_symb);
      }
      for ( it2=code_1[part].begin() ; it2 != code_1[part].end(); it2++ ){
      	// cout << *it2;
	 cur_length++;
	 }
      
    }
    
    cout << "  length = " << cur_length<<  "\n";
    
  }
  
     void string_encode_for_corr (std::vector<int> string_snp, std::vector<int> corr_snp,  std::map< std::vector<int>, std::vector<bool> > &code_1, int &cur_length, const int Number_of_ind) {
    std::vector<int>::iterator it;
    std::vector<bool>::iterator it2;
    std::vector<int> part;
   
    
    cout << "\nSource string1:  ";
    for ( it=string_snp.begin() ; it != string_snp.end(); it++ ){
      cout << *it;
    }
    cout << "  length = " << string_snp.size()<<  "\n";
    cout << "Source string2:  ";
    for ( it=corr_snp.begin() ; it != corr_snp.end(); it++ ){
      cout << *it;
    }
    cout << "  length = " << corr_snp.size()<<  "\n";
   
    
    cout << "Encoded string:  " ;
     for ( int i = 0; i <Number_of_ind; i++){
       //cout << *it;
       part.push_back(string_snp[i]);
       part.push_back(corr_snp[i]);
       
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
  
  bool correlation (const int snp1, const int snp2, vector < vector <double> > &prob_matrix, std::vector<int> &row1, std::vector<int> &row2, const int Number_of_ind) {
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
      double k = (prob_matrix.at(row1[i])).at(row2[i]);
      vector < double> kk = prob_matrix.at(row1[i]);
      
      prob_matrix[row1[i]][row2[i]] += (1.0/Number_of_ind);
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
   

   
   for (int i = start; i < start + window_size and i < Number_of_snps; i++) {
      row1.clear ();
      vector_for_snp (i, row1);
      cout << i << "\n";
      neighbours.push_back(rowww);
      for (int j = i; j < start + window_size and j < Number_of_snps; j++) { //++++ j=i
	row2.clear ();
        
	vector_for_snp (j, row2);
        if (correlation(i,j ,prob_matrix,row1,row2, Number_of_ind)) {
	   
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
    tree.string_encode_for_corr(row1,row2, code_1, cur_length, Number_of_ind);
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
  
  const int Number_of_snps = 545080;
  const int Number_of_ind = 2062;
  const int window_size = 0;
  


  long double whole_length = 0.0; 
  time_t start,end;
  start = time(NULL);
  
   ofstream result("result_plink.txt");  
   //make a head of output file
   result << "Snp\tfreq[0]\tfreq[1]\tfreq[2]\tfreq[3]\tblocks\tlength\tcost[0]\tcost[1]\tcost[2]\tcost[3]\n";
   
   std::cout << "Hello, world!" << std::endl;
   
   
   std::vector<int> row1; //will be string from file
   std::vector<int> row2; //
   std::vector<long double> individ_costs; //contained average number of bites for each individual
  
   
   for(int i =0; i < Number_of_ind; i++){
      individ_costs.push_back(0.0);
   }
  
  
   
  Genotypes* genotype = new Genotypes ("narac-plink.bed", Number_of_ind, Number_of_snps);
  Huff_tree <int, double>* tree = new Huff_tree <int, double>;
 
  vector<bool> snp_list; // 0 if snp was already processed, 1 if not
//  genotype->find_all_neighbours(Number_of_snps, Number_of_ind);
  
  for (int i = 0; i < Number_of_snps; i++) {
    snp_list.push_back(1) ;
  }
    
    

   for (int snp1 = 0; snp1 <5 ; snp1++){
     result << snp1 << "\t";
     
   if (snp_list[snp1]) {
     cout << "SNP:   " << snp1 << "\n";
     genotype->vector_for_snp (snp1, row1); //recieve data for SNP and put it in the vector<int> row1
     bool correlated = false; // = true, if we found correlating SNP
     
     vector<double> costs; //estimated bite costs of genotypes
     
       for (int snp2 = snp1+1; ((snp2 < snp1+window_size )and( snp2 < Number_of_snps)); snp2++) {
	
	 if (snp_list[snp2]){
	 vector< vector <double> > prob_matrix;
	 
	 //fill prob_matrix
	 vector<double> row;
	 for (int k = 0; k <4; k++){
	   row.push_back(0.0);
	 }
	 for (int l = 0; l<4; l++){
	   prob_matrix.push_back(row);
	 }
	 
	 
	 row2.clear();
	 genotype->vector_for_snp (snp2, row2); //take data for second snp
         
	 if(genotype-> correlation (snp1, snp2, prob_matrix, row1, row2, Number_of_ind)){
	    tree->make_leafs_for_2_snp(prob_matrix);
	 
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

    
  end = time(NULL);
 
  cout <<  difftime(end,start) << "\n";
  
  return 0;
    
}
