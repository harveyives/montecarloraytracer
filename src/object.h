/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

// Object is the base class for objects.
#ifndef _OBJECT_H_
#define _OBJECT_H_

#include <cfloat>
#include "ray.h"
#include "hit.h"
#include "vector.h"
#include "material.h"

class Object {
public:
    Material *material;
    // used for bounding boxes
    Vertex centre;
    float radius;
    long double min_limit_x = MAXFLOAT;
    long double min_limit_y = MAXFLOAT;
    long double min_limit_z = MAXFLOAT;
    long double max_limit_x = -MAXFLOAT;
    long double max_limit_y = -MAXFLOAT;
    long double max_limit_z = -MAXFLOAT;

	Object() = default;

    Object(Material *m) {
	    material = m;
	}
	
	virtual void intersection(Ray ray, Hit &hit)
	{

	}
};

#endif
