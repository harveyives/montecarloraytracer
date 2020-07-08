/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

// Transform is a class to store, manipulate and apply transforms.

#pragma once

#include "vertex.h"

class Transform {
public:
  
  float matrix[4][4];

  Transform()
  {
      matrix[0][0] = 1.0f;
      matrix[0][1] = 0.0f;
      matrix[0][2] = 0.0f;
      matrix[0][3] = 0.0f;
      matrix[1][0] = 0.0f;
      matrix[1][1] = 1.0f;
      matrix[1][2] = 0.0f;
      matrix[1][3] = 0.0f;
      matrix[2][0] = 0.0f;
      matrix[2][1] = 0.0f;
      matrix[2][2] = 1.0f;
      matrix[2][3] = 0.0f;
      matrix[3][0] = 0.0f;
      matrix[3][1] = 0.0f;
      matrix[3][2] = 0.0f;
      matrix[3][3] = 1.0f;
  };
  
  Transform(float a, float b, float c, float d,
	    float e, float f, float g, float h,
	    float i, float j, float k, float l,
	    float m, float n, float o, float p)
  {
      matrix[0][0] = a;
      matrix[0][1] = b;
      matrix[0][2] = c;
      matrix[0][3] = d;
      matrix[1][0] = e;
      matrix[1][1] = f;
      matrix[1][2] = g;
      matrix[1][3] = h;
      matrix[2][0] = i;
      matrix[2][1] = j;
      matrix[2][2] = k;
      matrix[2][3] = l;
      matrix[3][0] = m;
      matrix[3][1] = n;
      matrix[3][2] = o;
      matrix[3][3] = p;
  };

  void apply(Vertex &in, Vertex &out)
  {
    out.x = matrix[0][0] * in.x
          + matrix[0][1] * in.y
          + matrix[0][2] * in.z
          + matrix[0][3] * in.w;
    out.y = matrix[1][0] * in.x
          + matrix[1][1] * in.y
          + matrix[1][2] * in.z
          + matrix[1][3] * in.w;
    out.z = matrix[2][0] * in.x
          + matrix[2][1] * in.y
          + matrix[2][2] * in.z
          + matrix[2][3] * in.w;
    out.w = matrix[3][0] * in.x
          + matrix[3][1] * in.y
          + matrix[3][2] * in.z
          + matrix[3][3] * in.w;
  };
  
  void apply(Vertex &in)
  {
    float x,y,z,w;
    
    x = matrix[0][0] * in.x
      + matrix[0][1] * in.y
      + matrix[0][2] * in.z
      + matrix[0][3] * in.w;
    y = matrix[1][0] * in.x
      + matrix[1][1] * in.y
      + matrix[1][2] * in.z
      + matrix[1][3] * in.w;
    z = matrix[2][0] * in.x
      + matrix[2][1] * in.y
      + matrix[2][2] * in.z
      + matrix[2][3] * in.w;
    w = matrix[3][0] * in.x
      + matrix[3][1] * in.y
      + matrix[3][2] * in.z
      + matrix[3][3] * in.w;

    in.x = x;
    in.y = y;
    in.z = z;
    in.w = w;
  };
};
