/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

#pragma once

#include <math.h>

class Vector {
public:
	float x;
	float y;
	float z;

	Vector(float px, float py, float pz)
	{
		x = px;
		y = py;
		z = pz;
	}

	//Create vector from two points going point a - > b
    Vector(float ax, float ay, float az, float bx, float by, float bz)
    {
        x = (bx - ax);
        y = (by - ay);
        z = (bz - az);
    }

	Vector()
	{
		x = 0.0f;
		y = 0.0f;
		z = 0.0f;
	}

	void normalise()
	{
		float len = (float)sqrt((double)(x*x + y*y + z*z));
		x = x / len;
		y = y / len;
		z = z / len;
	}

	float dot(Vector &other)
	{
		return x*other.x + y*other.y + z*other.z;
	}


	void reflection(Vector initial, Vector &reflect)
	{
		float d;

		d = dot(initial);
		d = d * 2.0f;

		reflect.x = initial.x - d * x;
		reflect.y = initial.y - d * y;
		reflect.z = initial.z - d * z;
	}

	void negate()
	{
		x = -x;
		y = -y;
		z = -z;
	}

	void cross(Vector &other, Vector &result)
	{
	  result.x = y*other.z - z*other.y;
	  result.y = z*other.x - x*other.z;
	  result.z = x*other.y - y*other.x;
	}

    float magnitude()
    {
        return sqrt(x * x + y * y + z * z);
    }


    Vector operator * (float s)
    {
        return {s * x, s * y, s * z};
    }

    Vector operator - (Vector v)
    {
        return {x - v.x,y - v.y,z - v.z};
    }

    Vector operator + (Vector v)
    {
        return {x + v.x,y + v.y,z + v.z};
    }
};
