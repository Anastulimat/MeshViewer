#pragma once
#include <glm/glm.hpp>
#include "myFace.h"
#include "myVertex.h"

class myVertex;
class myFace;

class myHalfedge
{
public:
	myVertex *source; 
	myFace *adjacent_face; 
	myHalfedge *next;  
	myHalfedge *prev;  
	myHalfedge *twin;  

	size_t index; //Extra variable. Use as you wish!

	void copy(myHalfedge *);
	void reset();

	myHalfedge(myHalfedge *);
	myHalfedge();
	~myHalfedge();

	glm::vec3 edgePoint();
};