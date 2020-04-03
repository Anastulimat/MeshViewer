#include "myBVH.h"
#include "myAABBnode.h"
#include <algorithm>
#include <stack>
#include <iostream>
#include "errors.h"

using namespace std;

void myBVH::buildBVH(myMesh *m, const splittingMethod s)
{
	reset();

	for (size_t i=0; i < m->faces.size(); i++)
		sorted_allfaces.push_back(m->faces[i]);

	bvh_root = new myAABBnode();
	_build(0, sorted_allfaces.size()-1, bvh_root, 0, s);
}

void myBVH::_build(size_t left_index, size_t right_index, myAABBnode* root, size_t depth, const splittingMethod s)
{
	 
}

bool myBVH::intersect(mySegment ray, myFace * & picked_face, float & min_t, glm::mat4 model_matrix)
{
	return false;
}

bool myBVH::intersect(mySegment ray, myFace * & picked_face, float & min_t, std::tuple<int, int, int, int> & stats, glm::mat4 model_matrix)
{
	return false;
}

size_t myBVH::splitWithPivot(const myAABB & box, const size_t left_index, const size_t right_index, const splittingMethod s)
{
	return 0;
}


void myBVH::reset()
{
	 
}

myBVH::myBVH()
{
	bvh_root = nullptr;
}


myBVH::~myBVH()
{
	reset();
}
