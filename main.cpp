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
  
  
  void print(){
    
       if (exist_child){
         left_child->print();
         right_child->print();
	 
       }
       else{
         std::vector<bool>::iterator it;
         std::cout << data << "\t";
	 
         for ( it=encoded.begin() ; it != encoded.end(); it++ ){
	    std::cout << *it; 
         } 
	 std::cout <<"\n";  
    
       }   
    }
    
 }; // end of class Node

    
  Node* root;
  std::vector<Node*> body;
     
 
  
  void make_leafs(Huff_tree::Huff_map &letter_probability, int blocks) {
       int length = 0;
    //Huff_tree::Huff_map &leafs = * new Huff_tree::Huff_map();
    int max = (int)pow(4.0, (double)blocks);
    
    Frequency freq;
    for (int i = 0 ; i < max; i++){
      length = 0;
      int copy_i = i;
      freq = 1;
      while (length < blocks+1){
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
 
 #if 0
  void map_convertion (std::map<DataType, std::vector<bool> >& code, std::map<std::string, std::string>& code_in_string, int blocks) {
       std::vector<bool>::iterator it;
       int max = (int)pow(4.0, (double)blocks);
       std::string pref;
       
       for (int i = 0; i < max; i++) {
	 pref="";
	 for ( it=code[i].begin() ; it != code[i].end(); it++ ){
	   if (*it) {
	   pref.append("1");
	   }
	   else {
	    pref.append("0");
	   }
      	 }
	code_in_string[quaternary_convertion(i, blocks)] = pref; 
	}
	 
  }

  
  Huff_tree::Huff_map& combine_in_blocks(Huff_tree::Huff_map &letter_probability, int blocks){
            
    int length = 0;
    Huff_tree::Huff_map &leafs = * new Huff_tree::Huff_map();
    int max = (int)pow(4.0, (double)blocks);
    
    for (int i =0 ; i < max; i++){
      length = 0;
      int copy_i = i;
      leafs[i] = 1;
      while (length < blocks+1){
	if (copy_i > 0){
	    leafs[i] *= letter_probability[(copy_i % 4)];
	    copy_i /= 4;
         }	
	else{
	 leafs[i] *= letter_probability[0];
   
	}
	length++;
     }
   }   
     
   return leafs;
    
 }

#endif 
  
  void string_encoding(std::string initial, std::map<DataType, std::vector<bool> >& code) {
   
  }
  
  void calculate_expectation (){
    typename std::vector<Node*>::iterator iter;
    double expectation = 0;   
  for ( iter=body.begin() ; iter != body.end(); iter++ ){
    if (!(*iter)->exist_child){
    expectation += (*iter)->frequency*(*iter)->encoded.size();
     }
  }
    cout << "Expectation is" << "\t" << expectation << "\n";
   
  }
  
 
 
 std::vector<int> quaternary_convertion(int i, int blocks) {
  
    std::vector<int> quaternary;
    
      int copy_i = i;
      while (quaternary.size() < blocks+1){
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

  
};


  void print2(std::map<std::string, std::string>& code_in_string) {
    std::map<std::string, std::string>::iterator it;
      for ( it=code_in_string.begin() ; it != code_in_string.end(); it++ ){
    cout << (*it).first << "\t" << (*it).second << "\n";
       }
    
  }
  
int main(int argc, char **argv) {
    std::cout << "Hello, world!" << std::endl;
    Huff_tree <int, double>* tree = new Huff_tree <int, double>;
    std::map<int, double> letter_probability;
   
    letter_probability[0] = 0.1;
    letter_probability[1] = 0.4;
    letter_probability[2] = 0.24;
    letter_probability[3] = 0.26;
   std::vector<bool> prefix;
   std::map<int, std::vector<bool> > code ;
   std::map<std::string, std::string> code_in_string;
   
  int blocks = 4;
  tree->make_leafs(letter_probability, blocks); 
  tree->construct_tree()->fill(prefix,code);
  tree->root->print();
  tree->calculate_expectation();
   
  //tree->construct_tree(tree->combine_in_blocks(letter_probability, blocks))->fill(prefix,code);
  // tree->root->print(code);
   //tree->map_convertion(code, code_in_string, blocks);
  // tree->calculate_expectation(tree->combine_in_blocks(letter_probability, blocks),code);
  //print2 (code_in_string);
  return 0;
    
}
