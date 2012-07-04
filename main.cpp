#include <iostream>



template <class DataType, class Frequency> class Huff_tree <DataType, Frequency>{
  class Node {
    Frequency frequency;
    DataType data;
    Node* left_child;
    Node* right_child;
    
    Node (Node* left, Node* right){
      this.left_child = left;
      this.right_child =right;
      this.frequency =
      
    }
    
    Node (Frequency f, DataType d){
      this.frequency = f;
      this.data  = d;
    }
    void fill(vector <bool> prefix, map <DataType, vector <bool>> code){
      if (left_child){
	prefix.push_back(0);
	left_child->fill(prefix, code);
	prefix.pop_back();
	prefix.push_back(1);
	right_child->fill(prefix, code);
	
      }
      else{
	code[data] = prefix;
	
      }
      
      
    }
  }
int main(int argc, char **argv) {
    std::cout << "Hello, world!" << std::endl;
    return 0;

  
    
}
