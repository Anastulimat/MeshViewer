#include "myMesh.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <utility>
#include <glm/gtc/type_ptr.hpp>    
#include "helper_functions.h"
#include <wx/math.h>
#include <algorithm>
#include "errors.h"
#include "default_constants.h"
#include <unordered_set>


using namespace std;

myMesh::myMesh()
{
}


myMesh::~myMesh()
{
	clear();
}

void myMesh::clear()
{
	for (size_t i = 0; i < vertices.size(); i++) if (vertices[i]) delete vertices[i];
	for (size_t i = 0; i < halfedges.size(); i++) if (halfedges[i]) delete halfedges[i];
	for (size_t i = 0; i < faces.size(); i++) if (faces[i]) delete faces[i];

	vertices.clear();
	halfedges.clear();
	faces.clear();
}

myMesh::myMesh(std::vector<myVertex*> & in_vertices)
{
	for (size_t i = 0;i < in_vertices.size();i++)
		vertices.push_back(new myVertex(in_vertices[i]->point));
}

bool myMesh::checkMesh()
{
	cout << "Checking all mesh errors for " << name << ".\n";

	bool result = _checkfaces_boundaryedges() | _checkhalfedges_nextprev() | _checkhalfedges_source() | _checkhalfedges_twins() | _checkvertices_fans();

	cout << "Finished checking.\n";

	return result;
}


bool myMesh::_checkhalfedges_twins( )
{
	return _checkhalfedges_twins(halfedges);
}

bool myMesh::_checkhalfedges_twins(vector<myHalfedge *> & HE) const
{
	bool error = true;
	bool skip_further_errors = false;

	cout << "\tChecking the halfedge twin variables: " << HE.size() <<  endl;
	for (size_t i = 0; i < HE.size(); i++)
	{
		myassert(HE[i] != nullptr);

		if (HE[i]->twin == nullptr)
		{
			if (!skip_further_errors)
				cout << "\t\tWarning, boundary edges or otherwise nullptr twins.\n";
			skip_further_errors = true;
		}
		else
			myassert(HE[i] == HE[i]->twin->twin);
	}
	cout << "\tEnded check.\n\n";

	return error;
}

bool myMesh::_checkhalfedges_nextprev()
{
	return _checkhalfedges_nextprev(halfedges);
}

bool myMesh::_checkhalfedges_nextprev(vector<myHalfedge *> & HE) const
{
	bool error = true;

	cout << "\tChecking the halfedge next/prev variables: " << HE.size() << endl;
	for (size_t i = 0; i < HE.size(); i++)
	{
		myassert(HE[i] != nullptr);
		myassert(HE[i]->next != nullptr);
		myassert(HE[i]->prev != nullptr);
		myassert(HE[i]->next->prev == HE[i]);
		myassert(HE[i]->prev->next == HE[i]);
	}
	cout << "\tEnded check.\n\n";

	return error;
}


bool myMesh::_checkhalfedges_source()
{
	return _checkhalfedges_source(halfedges);
}
	
bool myMesh::_checkhalfedges_source(vector<myHalfedge *> & HE) const
{
	bool error = true;

	cout << "\tChecking the halfedge source variables: " << HE.size()  << endl;
	for (size_t i = 0; i < HE.size(); i++)
	{
		myassert(HE[i] != nullptr);
		myassert(HE[i]->source != nullptr);
	}
	cout << "\tEnded check.\n\n";

	return error;
}


bool myMesh::_checkvertices_fans()
{
	return _checkvertices_fans(vertices);
}

bool myMesh::_checkvertices_fans(vector<myVertex *> & V) const
{
	bool error = true;
	
	map<myVertex *, size_t> vertex_outdegrees;

	for (size_t i = 0;i < V.size();++i)
		vertex_outdegrees[V[i]] = 0;

	for (size_t i = 0;i < halfedges.size();++i)
		if (vertex_outdegrees.count(halfedges[i]->source))
			vertex_outdegrees[halfedges[i]->source]++;

	cout << "\tChecking fans of each vertex: " << V.size() << endl;
	for (size_t i = 0; i<V.size(); i++)
	{
		myassert(V[i] != nullptr);
		myassert(V[i]->originof != nullptr);

		myHalfedge *e = V[i]->originof;

		size_t k = 0;
		do 
		{
			myassert(e->prev != nullptr);
			e = e->prev->twin;

			k++;
			myassert(k <= MAX_VERTEX_DEGREE);
		} while (e != nullptr && e != V[i]->originof);
		myassert(k == vertex_outdegrees[V[i]]);
	}
	cout << "\tEnded check.\n\n";

	return error;
}

bool myMesh::_checkfaces_boundaryedges()
{
	return _checkfaces_boundaryedges(faces);
}

bool myMesh::_checkfaces_boundaryedges(vector<myFace *> & F) const
{
	bool error = true;
	bool istriangular = true;

	cout << "\tChecking edges of each face: " << F.size() << endl;

	size_t num_incidentedgesoverallfaces = 0;
	for (size_t i = 0; i<F.size(); i++)
	{
		myassert(F[i] != nullptr);

		myHalfedge *e = F[i]->adjacent_halfedge;
		size_t k = 0;
		do 
		{
			myassert(e != nullptr);
			e = e->next;

			k++;
			myassert(k <= MAX_FACE_DEGREE);
		} while (e != F[i]->adjacent_halfedge);

		num_incidentedgesoverallfaces += k;
		if (k>3) istriangular = false;
	}
	cout << "\tAverage number of edges per face: " << static_cast<float>(num_incidentedgesoverallfaces) / faces.size() << endl;
	if (istriangular) cout << "\t\tThe mesh is triangular.\n";
	else cout << "\t\tThe mesh is not triangular.\n";
	
	cout << "\tEnded check.\n\n";

	return error;
}



bool myMesh::readFile(std::string filename)
{
	string s, t, u, x, y, z;

	ifstream fin(filename);
	if (!fin.is_open()) 
	{
		PRINT(ERROR_FILEOPEN);
		return false;
	}
	name = filename;

	vector<glm::vec3> _points;
	vector<vector<size_t>> _faces;

	while (getline(fin, s))
	{
		stringstream myline(s);
		myline >> t;
		if (t == "v")
		{
			myline >> x >> y >> z;
			_points.push_back(glm::vec3(static_cast<float>(atof((x.substr(0, x.find("/"))).c_str())), static_cast<float>(atof((y.substr(0, y.find("/"))).c_str())), static_cast<float>(atof((z.substr(0, z.find("/"))).c_str()))));
		}
		if (t == "f")
		{
			_faces.push_back(vector<size_t>());
			while (myline >> u) _faces.back().push_back(atoi((u.substr(0, u.find("/"))).c_str()) - 1);
		}
		t = "";
	}

	/* MY CODE START HERE */

	/* CREATE VERTICES FROM POINTS AND ADD THEME */
	for (int i = 0; i < _points.size(); i++)
	{
		vertices.push_back(new myVertex(_points[i]));
	}

	map<pair<int, int>, myHalfedge*> table;

	for (int i = 0; i < _faces.size(); i++)
	{
		vector<size_t> currentFace = _faces[i];
		myFace* f = new myFace();
		vector<myHalfedge*> hes(currentFace.size());

		for (int j = 0; j < currentFace.size(); j++)
		{
			hes[j] = new myHalfedge();
		}

		f->adjacent_halfedge = hes[0];

		for (int k = 0; k < currentFace.size(); k++)
		{
			int ipo = (k + 1) % currentFace.size();
			int imo = (k - 1 + currentFace.size()) % currentFace.size();

			table[make_pair(currentFace[k], currentFace[ipo])] = hes[k];
			map<pair<int, int>, myHalfedge*>::iterator it = table.find(make_pair(currentFace[ipo], currentFace[k]));

			hes[k]->next = hes[ipo];
			hes[k]->prev = hes[imo];
			//hes[k]->twin = nullptr;

			if (it != table.end())
			{
				hes[k]->twin = it->second;
				it->second->twin = hes[k];
			}

			hes[k]->source = vertices[currentFace[k]];
			hes[k]->adjacent_face = f;

			vertices[currentFace[k]]->originof = hes[k];

			halfedges.push_back(hes[k]);
		}
		faces.push_back(f);
	}

	/* MY CODE END HERE */

	
 
	checkMesh();

	computeNormals();

	normalize();

	return true;
}

void myMesh::computeNormals()
{

	/* MY CODE START HERE */

	for (int i = 0; i < faces.size(); i++)
	{
		faces[i]->computeNormal();
	}

	for (int i = 0; i < vertices.size(); i++)
	{
		vertices[i]->computeNormal();
	}

	/* MY CODE END HERE */
}

void myMesh::normalize()
{
	if (vertices.size() == 0) return;

	enum { MIN, MAX };
	vector<glm::vec3> corner = { vertices[0]->point, vertices[0]->point };

	for (size_t i = 1; i < vertices.size(); i++) 
	{
		for (int coordinate = 0;coordinate < 3; coordinate++)
		{
			corner[MIN][coordinate] = std::min( corner[MIN][coordinate], vertices[i]->point[coordinate] );
			corner[MAX][coordinate] = std::max( corner[MAX][coordinate], vertices[i]->point[coordinate] );
		}
	}
	
	float scale_factor = std::max({ corner[MAX][0] - corner[MIN][0], corner[MAX][1] - corner[MIN][1], corner[MAX][2] - corner[MIN][2] } );

	for (size_t i = 0; i < vertices.size(); i++) 
	{
		for (int coordinate = 0;coordinate < 3; coordinate++)
			vertices[i]->point[coordinate] -= (corner[MAX][coordinate] + corner[MIN][coordinate]) / 2.0f;
		vertices[i]->point /= scale_factor;
	}
}

bool myMesh::writeFile(std::string filename)
{
	ofstream fout(filename);

	if (!fout.is_open())
	{
		PRINT(ERROR_FILEOPEN);
		return false;
	}

	for (size_t i = 0;i < vertices.size(); i++)
	{
		fout << "v " << vertices[i]->point.x << " " << vertices[i]->point.y << " " << vertices[i]->point.z << endl;
		vertices[i]->index = i + 1;
	}

	for (size_t i = 0;i < faces.size(); i++)
	{
		fout << "f";
		myHalfedge *e = faces[i]->adjacent_halfedge;
		do
		{
			fout << " " << e->source->index;
			e = e->next;
		} while (e != faces[i]->adjacent_halfedge);
		fout << endl;
	}
	fout.close();
	return true;
}


vector<glm::vec3> myMesh::voronoiReconstruction( )
{
	 
	return vector<glm::vec3>( );
} 

void myMesh::splitFaceTRIS(myFace *f, glm::vec3 p)
{
	if (f == nullptr) return;

	/**************TODO*************************/

}


bool myMesh::splitEdge(myHalfedge *e1, glm::vec3 p)
{
	if (e1 == nullptr)
	{
		PRINT(ERROR_NULLHE);
		return false;
	}

	/*MY CODE START HERE*/
	myVertex* v = new myVertex();
	myHalfedge* e2 = new myHalfedge();
	myHalfedge* e2_twin = new myHalfedge();

	/* Définir les attributs de la nouvelle vertice */
	v->point = p;
	v->originof = e2;

	/* Définir les attrbuts du nouvel edge e2 */
	e2->source = v;
	e2->prev = e1;
	e2->next = e1->next;
	e2->twin = e2_twin;
	e2->adjacent_face = e1->adjacent_face;

	/* Définir les attrbuts du nouvel edge e2_twin */
	e2_twin->source = e1->twin->source;
	e2_twin->next = e1->twin;
	e2_twin->prev = e1->twin->prev;
	e2_twin->twin = e2;
	e2_twin->adjacent_face = e1->twin->adjacent_face;

	/* Ajuster le reste des élément en fonction des changement */
	e1->twin->source = v;
	e1->next->prev = e2;
	e1->next = e2;

	e1->twin->prev->next = e2_twin;
	e1->twin->prev = e2_twin;

	vertices.push_back(v);
	halfedges.push_back(e2);
	halfedges.push_back(e2_twin);

	/*MY CODE END HERE*/


	return false;
}


void myMesh::splitFaceQUADS(myFace *f, glm::vec3 p)
{
	
	int nbEdges = f->size()/2;
	int nbNewFace = nbEdges - 1;


	myVertex* v = new myVertex();
	v->point = p;

	vector<myHalfedge*> in = vector<myHalfedge*>();
	vector<myHalfedge*> out = vector<myHalfedge*>();

	vector<myHalfedge*> E = vector<myHalfedge*>();
	vector<myHalfedge*> nextE = vector<myHalfedge*>();

	vector<myFace*> nf = vector<myFace*>();

	for (int i = 0; i < nbEdges; i++)
	{
		in.push_back(new myHalfedge());
		out.push_back(new myHalfedge());

	}

	myHalfedge* e = f->adjacent_halfedge;
	for (int i = 0; i < nbEdges; i++)
	{
		E.push_back(e);
		e = e->next;
		nextE.push_back(e);
		e = e->next;
	}

	nf.push_back(f);
	for (int i = 0; i < nbNewFace; i++)
	{
		nf.push_back(new myFace());
	}
	for (int i = 0; i < nbEdges; i++)
	{
		int ipo = (i + 1) % nbEdges;
		int imo = (i - 1 + nbEdges) % nbEdges;
			
		E[i]->adjacent_face = nf[i];
		E[i]->prev = out[i];

		out[i]->source = v;
		out[i]->next = E[i];
		out[i]->prev = in[i];
		out[i]->adjacent_face = nf[i];
		out[i]->twin = in[ipo];

		in[i]->source = nextE[i]->next->source;
		in[i]->next = out[i];
		in[i]->prev = nextE[i];
		in[i]->adjacent_face = nf[i];
		in[i]->twin = out[ipo];

		nextE[i]->adjacent_face = nf[i];
		nextE[i]->next = in[i];

		nf[i]->adjacent_halfedge = nextE[i];
		v->originof = out[i];
	}

	vertices.push_back(v);
	for (int i = 0; i < nbNewFace; i++)
	{
		faces.push_back(nf[1 + i]);
	}

	for (int i = 0; i < nbEdges; i++)
	{
		halfedges.push_back(in[i]);
		halfedges.push_back(out[i]);

	}
}


void myMesh::splitFace_size6(myFace *f, myHalfedge *starting_edge)
{
	/**************TODO*************************/
	
}



void myMesh::fractalize()
{
	/**************TODO*************************/

}


void myMesh::subdivisionCatmullClark_createNewPoints(std::vector<glm::vec3> & facepoints,
													 std::vector<glm::vec3> & edgepoints,
													 std::vector<glm::vec3> & newvertexpoints)
{
	facepoints.clear();
	edgepoints.clear();
	newvertexpoints.clear();
	/*MY CODE START HERE*/
	for (int i = 0; i < faces.size(); i++)
	{
		facepoints.push_back(faces[i]->centroid());
	}

	for (int i = 0; i < halfedges.size(); i++)
	{
		halfedges[i]->index = i+1;
		halfedges[i]->twin->index = i + 1;
		edgepoints.push_back(halfedges[i]->edgePoint());
	}

	for (int i = 0; i < vertices.size(); i++)
	{
		newvertexpoints.push_back(vertices[i]->vertexPoint());
	}

	/*MY CODE END HERE*/

}

void myMesh::subdivisionCatmullClark()
{
	/*MY CODE START HERE*/
	std::vector<glm::vec3> facepoints;
	std::vector<glm::vec3> edgepoints;
	std::vector<glm::vec3> newvertexpoints;

	subdivisionCatmullClark_createNewPoints(facepoints, edgepoints, newvertexpoints);

	// Update position of the old vertices
	for (size_t i = 0; i < vertices.size(); i++) {
		vertices[i]->point = newvertexpoints[i];
	}

	// Update position of the old vertices
	size_t num_hedges = halfedges.size();
	for (size_t i = 0; i < num_hedges; i++)
	{
		if (halfedges[i]->index == halfedges[i]->twin->index)
		{
			splitEdge(halfedges[i], edgepoints[i]);
			halfedges[i]->index = -1;
		}
		
	}

	// Add new vertices for each face
	size_t num_faces = faces.size();
	for (size_t i = 0; i < num_faces; i++) {
		splitFaceQUADS(faces[i], facepoints[i]);
	}
	/*MY CODE END HERE*/

}

void myMesh::setIndices()
{
	for (size_t i = 0; i < faces.size(); i++) 
		faces[i]->index = i;
	for (size_t i = 0; i < halfedges.size(); i++) 
		halfedges[i]->index = i;
	for (size_t i = 0; i < vertices.size(); i++)	
		vertices[i]->index = i;
}

void myMesh::triangulate()
{
	for (size_t i = 0; i < faces.size(); i++)
		triangulate(faces[i]);
}

//return false if already triangle, true othewise.
//new edge and face array stores an extra index  for shorter code.
bool myMesh::triangulate(myFace *f)
{
	myassert(f != nullptr);

	size_t faceSize = f->size();

	/*MY CODE START HERE*/
	if (faceSize == 3) {
		return false;
	}

	/* Déclaration les nouveaux élements à ajouter à nos tableaux 
	   Et qui vont nous servir à découper la face
	*/
	myFace* newFace = new myFace();
	myHalfedge* newEdge = new myHalfedge();
	myHalfedge* newEdgeTwin = new myHalfedge();

	while (f->size() > 3)
	{
		/* Définir les propriétés du nouveau Edge */
		newEdge->source = f->adjacent_halfedge->next->next->source;
		newEdge->next = f->adjacent_halfedge;
		newEdge->prev = f->adjacent_halfedge->next;
		newEdge->adjacent_face = newFace;
		newEdge->twin = newEdgeTwin;

		/* Définir les propriétés du nouveau EdgeTwin */
		newEdgeTwin->source = f->adjacent_halfedge->source;
		newEdgeTwin->next = f->adjacent_halfedge->next->next;
		newEdgeTwin->prev = f->adjacent_halfedge->prev;
		newEdgeTwin->adjacent_face = f;
		newEdgeTwin->twin = newEdge;

		/* Ajuster le reste des éléments en fonction du changement */
		f->adjacent_halfedge->next->next = newEdge;
		f->adjacent_halfedge->prev->next = newEdgeTwin;
		f->adjacent_halfedge->next->next->prev = newEdgeTwin;
		f->adjacent_halfedge->prev = newEdge;
		f->adjacent_halfedge->adjacent_face = newFace;

		newFace->adjacent_halfedge = f->adjacent_halfedge;

		f->adjacent_halfedge = newEdgeTwin;

		faces.push_back(newFace);
		halfedges.push_back(newEdge);
		halfedges.push_back(newEdgeTwin);
	}
	/*MY CODE END HERE*/

	return true;
}

void myMesh::copy(myMesh *m)
{
	for (size_t i = 0;i < m->vertices.size();i++)
	{
		m->vertices[i]->index = i;
		vertices.push_back(new myVertex());
	}

	for (size_t i = 0; i < m->halfedges.size();i++)
	{
		m->halfedges[i]->index = i;
		halfedges.push_back(new myHalfedge());
	}

	for (size_t i = 0;i < m->faces.size();i++)
	{
		m->faces[i]->index = i;
		faces.push_back(new myFace());
	}

	for (size_t i = 0;i < m->vertices.size(); i++)
	{
		vertices[i]->point = glm::vec3(m->vertices[i]->point);
		vertices[i]->normal = glm::vec3(m->vertices[i]->normal);
		vertices[i]->originof = halfedges[m->vertices[i]->originof->index];
	}

	for (size_t i = 0; i < halfedges.size();i++)
	{
		halfedges[i]->source = vertices[m->halfedges[i]->source->index];
		halfedges[i]->adjacent_face = faces[m->halfedges[i]->adjacent_face->index];
		halfedges[i]->next = halfedges[m->halfedges[i]->next->index];
		halfedges[i]->prev = halfedges[m->halfedges[i]->prev->index];
		if ( m->halfedges[i]->twin != nullptr && m->halfedges[i]->twin->index < m->halfedges.size() ) 
			halfedges[i]->twin = halfedges[m->halfedges[i]->twin->index];
	}

	for (size_t i = 0;i < m->faces.size();i++)
	{
		faces[i]->normal = glm::vec3(faces[i]->normal);
		faces[i]->adjacent_halfedge = halfedges[m->faces[i]->adjacent_halfedge->index];
	}
}

bool myMesh::intersect(mySegment ray, myFace * & picked_face, float & min_t, glm::mat4 model_matrix)
{
	std::tuple<int, int, int, int> tmp_stats;
	return intersect(ray, picked_face, min_t, tmp_stats, model_matrix);
}

bool myMesh::intersect(mySegment ray, myFace * & picked_face, float & min_t, std::tuple<int, int, int, int> & stats, glm::mat4 model_matrix)
{
	bool result = false;

	get<0>(stats) = 0;
	get<1>(stats) = 0;
	get<2>(stats) = 0;

	picked_face = nullptr;

	min_t = std::numeric_limits<float>::max();
	float tmp;

	for (size_t i = 0;i < faces.size();i++)
	{
		myassert(faces[i] != nullptr);
		if (faces[i]->intersect(ray, tmp, model_matrix) && tmp < min_t)
		{
			min_t = tmp;
			picked_face = faces[i];
			result = true;
			get<1>(stats)++;
		}
		get<0>(stats)++;
	}
	return result;
}

void myMesh::smoothen(float delta)
{
	/**************TODO*************************/

}

void myMesh::smoothenVertices(vector<myVertex *> V, float delta)
{
	/**************TODO*************************/


}

void myMesh::sharpen(float delta)
{
	/**************TODO*************************/

}

void myMesh::inflate(float delta)
{
	/*MY CODE START HERE*/
	for (myVertex* v : vertices)
	{
		v->point += v->normal * delta;
	}
	/*MY CODE START HERE*/

}

void myMesh::bumpen(myVertex *v, float delta, int size)
{
	/**************TODO*************************/

}

void myMesh::computeSilhouette(std::vector<myHalfedge *>& silhouette_edges, glm::vec3 camera_position)
{
	/* MY CODE START HERE */

	for (int i = 0; i < halfedges.size(); i++)
	{
		myHalfedge* e = halfedges[i];
		if (e->twin != nullptr)
		{
			myHalfedge* twin = e->twin;
			float x = glm::dot(e->adjacent_face->normal, camera_position - e->source->point);
			float y = glm::dot(twin->adjacent_face->normal, camera_position - twin->source->point);

			if (x <= 0 && y >= 0)
			{
				silhouette_edges.push_back(e);
			}
		}
	}

	/* MY CODE END HERE */

}
