/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

#pragma once

#include <cmath>
#include "vector.h"

class Vertex {
public:
	float x;
	float y;
	float z;
	float w;

	Vertex()
	{
		x = 0.0;
		y = 0.0;
		z = 0.0;
		w = 1.0f;
	}

	Vertex(float px, float py, float pz, float pw)
	{
		x = px;
		y = py;
		z = pz;
		w = pw;
	}

	Vertex(float px, float py, float pz)
	{
		x = px;
		y = py;
		z = pz;
		w = 1.0f;
	}

	Vector toVector() {
	    return {x,y,z};
	}

    Vector operator - (Vertex v)
    {
        return {x - v.x,y - v.y,z - v.z};
    }

    Vertex operator + (Vertex v)
    {
        return {x + v.x,y + v.y,z + v.z};
    }

    Vertex operator + (Vector v)
    {
        return {x + v.x,y + v.y,z + v.z};
    }

    Vertex operator / (float s)
    {
        return {x / s,y / s,z / s};
    }

    float magnitude()
    {
        return sqrt(x * x + y * y + z * z);
    }
};
