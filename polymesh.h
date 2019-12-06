/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

#pragma once


#include "object.h"
#include "vertex.h"
#include "transform.h"
#include "ray.h"
#include "hit.h"
#include "vector.h"
#include <vector>
#include <string>
#include <cfloat>

typedef int TriangleIndex[3];

class PolyMesh : public Object {
public:
	int vertex_count;
	int triangle_count;
	Vertex *vertices;
	TriangleIndex *triangles;
    Vector *face_normal;
    Vector *vertex_normal;
	Object *bounding_sphere;

    void do_construct(char *file, Transform *transform);
    PolyMesh(char *file);
    PolyMesh(char *file, Transform *transform);

    PolyMesh(char *file, Transform *transform, Material *material);

    void find_bounding_sphere_values(Vertex &c, float &r);

    void intersection(Ray ray, Hit &hit) override;

    void compute_vertex_normals();
    void compute_face_normal(int which_triangle, Vector &normal);

    void compute_barycentric(Vertex p, Vertex a, Vertex b, Vertex c, float &u, float &v, float &w);
};
