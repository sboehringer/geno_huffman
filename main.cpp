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
using namespace std;

template <typename DataType, typename Frequency> 
class Huff_tree{
  
    
  public:
    
  typedef typename std::map< DataType,  Frequency>	Huff_map;
  typedef typename std::map< DataType,  Frequency>::iterator Huff_iterator;
  
  
  
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
    
    Node (Node* left, Node* right){
      left_child = left;
      right_child = right;
      exist_child = true;
      this->frequency = left_child->frequency + right_child->frequency;
      data = 0;         
      quaternary_data;
      encoded;
      
    }
    
    Node (Frequency f, DataType d, std::vector<int> quater){
      frequency = f;
      data  = d;
      this->left_child = 0;
      this->right_child = 0;
      exist_child = false;
      quaternary_data = quater;
      encoded;
    }
    
    Node(){
      frequency = 0;            
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
      
    encoded = prefix;
    }
         
  }
  
  //============================================
  void print(){
    
       if (exist_child){
         left_child->print();
         right_child->print();
	 
       }
       else{
         std::vector<bool>::iterator it;
	 std::vector<int>::iterator it2;
         for ( it2=quaternary_data.begin() ; it2 != quaternary_data.end(); it2++ ){
	    std::cout << *it2; 
         } 
	std::cout << "\t";
	 
         for ( it=encoded.begin() ; it != encoded.end(); it++ ){
	    std::cout << *it; 
         } 
         std::cout << "\t"<< frequency << "\n";
	// std::cout <<"\n";  
    
       }   
    }
    
 
   
}; // end of class Node
//======================================================
//=======================================================
    
 
 
 
 Node* root;
  std::vector<Node*> body;
     
 
  //=====================================================================
  void make_leafs(Huff_tree::Huff_map &letter_probability, int blocks) {
       int length = 0;
    //Huff_tree::Huff_map &leafs = * new Huff_tree::Huff_map();
    int max = (int)pow(4.0, (double)blocks);
    
    Frequency freq;
    for (int i = 0 ; i < max; i++){
      length = 0;
      int copy_i = i;
      freq = 1;
      while (length < blocks){
	if (copy_i > 0){
	    freq *= letter_probability[(copy_i % 4)];
	    
	    copy_i /= 4;
         }	
	else{
	 freq *= letter_probability[0];
   
	}
	length++;
     }
    Node* dataNode = new Node(freq, i, quaternary_convertion(i , blocks)); 
     body.push_back(dataNode);
      
    }   
   
}  

//=======================================================================
Node* construct_tree(){
    priority_queue<Node*, vector<Node*>, Node> pqueue;
        
     Huff_tree::Huff_iterator it;
    typename std::vector<Node*>::iterator iter;
     
  for (iter=body.begin(); iter != body.end(); iter++ ){
      pqueue.push(*iter);
  }
    
 
  while (!pqueue.empty()){
    
       Node* top = pqueue.top();
       pqueue.pop();
       if (pqueue.empty()){
          root = top;
       }
       else {
        Node* top2 = pqueue.top();
        pqueue.pop();
        pqueue.push(new Node(top, top2));
       }
   }
  // std::cout << "find_root" << std::endl;
   int size = pqueue.size();
       return root; 
}
//=========================================================================== 
  
  
  
  void calculate_expectation (){
    typename std::vector<Node*>::iterator iter;
    double expectation = 0;   
  for ( iter=body.begin() ; iter != body.end(); iter++ ){
    if (!(*iter)->exist_child){
    expectation += (*iter)->frequency*(*iter)->encoded.size();
     }
  }
   std::cout << "=========================================" << "\n";
    cout << "Expectation is" << "\t" << expectation << "\n";
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

  void bite_cost () {
    double weight_old [4] = {1.0, 1.0, 1.0, 1.0};
    double weight [4];
    for (int i = 0; i<20; i++) {
      double weight [4] = {0.0, 0.0, 0.0, 0.0};
      typename std::vector<Node*>::iterator iter;
      for ( iter=body.begin() ; iter != body.end(); iter++ ){
	if (!(*iter)->exist_child){
	  double block_weight = 0;
	  std::vector<int>::iterator it2;
	  for ( it2=(*iter)->quaternary_data.begin() ; it2 != (*iter)->quaternary_data.end(); it2++ ){
	      block_weight +=  weight_old[(*it2)];
	  }
	  for ( it2=(*iter)->quaternary_data.begin() ; it2 != (*iter)->quaternary_data.end(); it2++ ){
	      weight[(*it2)] +=  (*iter)->frequency*(weight_old[(*it2)]/block_weight)*(*iter)->encoded.size();
	  }
        }
      }
      std::cout << weight[0] << "\t" << "\t" << weight[1] << "\t" <<"\t" <<  weight[2] << "\t" << "\t" << weight[3] << "\n";
      weight_old[0]= weight [0];
      weight_old[1]= weight [1];
      weight_old[2]= weight [2];
      weight_old[3]= weight [3];
      weight_old= weight; 
      
    }
     std::cout << "=========================================" << "\n";
  }
  
  void bite_cost_2 (int blocks) {
    double weight_old [4] = {1.0, 1.0, 1.0, 1.0};
    double weight [4];
    for (int i = 0; i<20; i++) {
      double weight [4] = {0.0, 0.0, 0.0, 0.0};
      typename std::vector<Node*>::iterator iter;
      for ( iter=body.begin() ; iter != body.end(); iter++ ){
	if (!(*iter)->exist_child){
	  double block_weight = 0;
	  std::vector<int>::iterator it2;
	  for ( it2=(*iter)->quaternary_data.begin() ; it2 != (*iter)->quaternary_data.end(); it2++ ){
	      block_weight +=  weight_old[(*it2)];
	  }
	  weight[(*iter)->quaternary_data[0]]+= (*iter)->frequency*(weight_old[(*iter)->quaternary_data[0]]/block_weight)*(*iter)->encoded.size()*blocks;
	  //for ( it2=(*iter)->quaternary_data.begin() ; it2 != (*iter)->quaternary_data.end(); it2++ ){
	   //   weight[(*it2)] +=  (*iter)->frequency*(weight_old[(*it2)]/block_weight)*(*iter)->encoded.size();
	 // }
        }
      }
      std::cout << weight[0] << "\t" << "\t" << weight[1] << "\t" <<"\t" <<  weight[2] << "\t" << "\t" << weight[3] << "\n";
      weight_old[0]= weight [0];
      weight_old[1]= weight [1];
      weight_old[2]= weight [2];
      weight_old[3]= weight [3];
    }
  }
  
};// end of class Huff_tree;




  
int main(int argc, char **argv) {
    std::cout << "Hello, world!" << std::endl;
    Huff_tree <int, double>* tree = new Huff_tree <int, double>;
    std::map<int, double> letter_probability;
   
    letter_probability[0] = 0.4;
    letter_probability[1] = 0.3;
    letter_probability[2] = 0.2;
    letter_probability[3] = 0.1;
   std::vector<bool> prefix;
   std::map<int, std::vector<bool> > code ;
   
   
  int blocks = 5;
  tree->make_leafs(letter_probability, blocks); 
  tree->construct_tree()->fill(prefix,code);
  tree->root->print();
  tree->calculate_expectation();
  tree->bite_cost(); 
  tree->bite_cost_2(blocks);
  return 0;
    
}
