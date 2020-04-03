#include "myVertex.h"
#include "myHalfedge.h"
#include "myFace.h"
#include <glm/glm.hpp>
#include "errors.h"


myVertex::myVertex(myVertex *v)
{
	copy(v);
}

myVertex::myVertex(glm::vec3 p)
{
	reset();
	point = p;
}

myVertex::~myVertex()
{
}

void myVertex::reset()
{
	originof = nullptr;
	index = 0;
}

void myVertex::computeNormal()
{
	myassert(originof != nullptr);
	
	/* MY CODE START HERE */
	normal = glm::vec3(0.0f, 0.0f, 0.0f);
	myHalfedge* e = originof;

	do
	{
		normal += e->adjacent_face->normal;
		e = e->prev->twin;
	}
	while(e != nullptr && e != originof);
	normal = glm::normalize(normal);
	/* MY CODE END HERE */
}

void myVertex::copy(myVertex *v)
{
	myassert(v != nullptr);

	point = v->point;
	originof = v->originof;
	index = v->index;
}

glm::vec3 myVertex::vertexPoint()
{
	glm::vec3 q = glm::vec3(0.0f, 0.0f, 0.0f);
	glm::vec3 r = glm::vec3(0.0f, 0.0f, 0.0f);

	myHalfedge* e = this->originof;
	float n = 0;
	do
	{
		q += e->adjacent_face->centroid();
		r += (e->source->point + e->twin->source->point) / 2.0f;
		e = e->twin->next;
		n++;
	} while (e != this->originof);

	q = q / n;
	r = r / n;
	glm::vec3 temp = (1.0f / n) * q + (2.0f / n) * r + ((n - 3.0f / n) * this->point);
	glm::vec3 tempV = this->point;
	
	glm::vec3 temp2 = (1.0f / n) * q;
	glm::vec3 temp3 = (2.0f / n) * r;
	glm::vec3 temp4 = (((n - 3.0f) / n) * this->point);

	return ((1.0f / n) * q) + ((2.0f / n) * r) + (((n - 3.0f) / n) * this->point);




}
