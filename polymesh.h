/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

#pragma once

#include "vertex.h"
#include "transform.h"

typedef int TriangleIndex[3];

class PolyMesh {
public:
	int vertex_count;
	int triangle_count;
        Vertex *vertex;
	TriangleIndex *triangle;

	void do_construct(char *file, Transform *transform);
	
	PolyMesh(char *file);
	PolyMesh(char *file, Transform *transform);
};
