/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

// Object is the base class for objects.
#ifndef _OBJECT_H_
#define _OBJECT_H_

#include "ray.h"
#include "hit.h"
#include "vector.h"
#include "material.h"

class Object {
public:
    Material material;

	Object() = default;

	Object(Material m) : material(Vector(), 0, 0, 0) {
	    material = m;
	}
	
	virtual void intersection(Ray ray, Hit &hit)
	{

	}
};

#endif
