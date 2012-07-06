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
 class Node;  

 class Node {
    public:
    Frequency frequency;
    DataType data;
    Node* left_child;
    Node* right_child;
    bool exist_child;
    
    Node (Node* left, Node* right){
      left_child = left;
      right_child = right;
      exist_child = true;
      this->frequency = left_child->frequency + right_child->frequency;
      data = 0;         
    }
    
    Node (Frequency f, DataType d){
      frequency = f;
      data  = d;
      this->left_child = 0;
      this->right_child = 0;
      exist_child = false;
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
    
  
  
  
  std::map<DataType, std::vector<bool> >&  fill(std::vector<bool> prefix, std::map<DataType, std::vector<bool> >& code){
    
   //  std::cout << "fill tree!"  << std::endl;
      
    if (exist_child){
       //std::cout << "fill tree!"  << std::endl;
      prefix.push_back(0);
      left_child->fill(prefix, code);
      prefix.pop_back();
      prefix.push_back(1);
      right_child->fill(prefix, code);
	
    }
    else{
    code[data] = prefix;
    //std::cout <<  code[data] << std::endl;
    
      }
      
     // std::cout << "fill tree!" << std::endl;
      return code;
   
    
  }
  
    void print(std::map<DataType, std::vector<bool> >& code){
    
   //  std::cout << "fill tree!"  << std::endl;
          
    
    if (exist_child){
      
      left_child->print(code);
       right_child->print(code);
	
    }
    else{
       std::vector<bool>::iterator it;
       std::cout << data << "\t";
      for ( it=code[data].begin() ; it != code[data].end(); it++ ){
      
	std::cout << *it; 
      
      }
      
        std::cout <<"\n";  
    
      }
      
      //std::cout << "fill tree!" << std::endl;
      
   
    
  }
    
  };

    
  Node* root;
  
  typedef typename std::map< DataType,  Frequency>	Huff_map;
  typedef typename std::map< DataType,  Frequency>::iterator Huff_iterator;
     

  
  //void construct_tree(std::map< DataType, Frequency> &leafs);
  
    
    
Node* construct_tree(Huff_tree::Huff_map &leafs){
    priority_queue<Node*, vector<Node*>, Node> pqueue;
   // std::cout << pqueue.size() << "\n";     
    
     Huff_tree::Huff_iterator it;
     //map<char *, double>::iterator it;

  
  for ( it=leafs.begin() ; it != leafs.end(); it++ ){
    
      Node* dataNode = new Node((*it).second,(*it).first);
      pqueue.push(dataNode);
  }
  while (!pqueue.empty()){
    //std::cout << pqueue.size() << std::endl;
       Node* top = pqueue.top();
       pqueue.pop();
       if (pqueue.empty()){
          root = top;
	  
	  // top->fill;
       }
       else {
        Node* top2 = pqueue.top();
        pqueue.pop();
        pqueue.push(new Node(top, top2));
       }
   }
   std::cout << "find_root" << std::endl;
   int size = pqueue.size();
    
     return root; 
}
 
  
  Huff_tree::Huff_map& combine_in_blocks(Huff_tree::Huff_map &letter_probability, int blocks){
    Huff_tree::Huff_map &leafs = * new Huff_tree::Huff_map();
    int max = (int)pow(4.0, (double)blocks);
    for (int i =0 ; i < max; i++){
        if(i < 4){
	  leafs[i] = letter_probability[i]*pow(letter_probability[0], blocks - 1); 
	  
	}
	else{ 
	  int copy_i = i; 
	  leafs[i] = 1;
	  while (copy_i > 0){
	    leafs[i] *= letter_probability[(copy_i % 4)];
	    copy_i /= 4;
          }	
      
        }
    
        
    
      }   
    std::cout << leafs.size() << "\n";   
    return leafs;
    
  }
};


  
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
   
   tree->construct_tree(tree->combine_in_blocks(letter_probability, 2))->fill(prefix,code);
   tree->root->print(code);
 
  return 0;
    
}
