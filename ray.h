/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

// Ray is a class to store and maniplulate 3D rays.

#pragma once

#include "vertex.h"
#include "vector.h"

class Ray {
public:
	Vertex position;
	Vector direction;

    Ray() {}

	Ray(Vertex p, Vector d)
	{
		position = p;
		direction = d;
	}
};
