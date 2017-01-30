#include <iostream>
#include <fstream>
#include <string>
#include <vector>


// ROOT //
#include "TFile.h"
#include "TTree.h"



using namespace std;



class Rectangle{
	//members of the class
	int width, height;
	public:
	Rectangle(int,int); //Constructor
	int area(){return width*height;}
	
	};
	
	
	/*Scope Operator :: define a member of the class 
	 * outside 
	 */ 
	Rectangle::Rectangle(int x, int y){  //instead of void Rectangle :: set_vales
		width = x;
		height = y;
		}
	
	
	
	/*
	 * n order to avoid that, a class can include a special function called its constructor, which is automatically called whenever a new object of this class is created, allowing the class to initialize member variables or allocate storage.

This constructor function is declared just like a regular member function, but with a name that matches the class name and without any return type; not even void.

The Rectangle class above can easily be improved by implementing a constructor:
	 * 
	 * */
	
	
	
	
	 //obj of class like int a
	/*After declaraiton of the class, any of the public members of object rect can be accessed s if they were a normal function 
	 * .between object ame and member name: rect.set_values();
myarea = rect.area(); */




int main(){
	Rectangle rect (3,4);
	cout<<"area: "<<rect.area();
	return 0;
	
	
	
	}
