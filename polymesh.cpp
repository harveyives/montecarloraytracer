/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <float.h>
#include <vector>
#include<algorithm>
#include <values.h>
#include "vector.h"

#include "polymesh.h"
#include "sphere.h"

using namespace std;

PolyMesh::PolyMesh(char *file) : PolyMesh(file, new Transform()){

}

PolyMesh::PolyMesh(char *file, Transform *transform) : PolyMesh(file, transform, Vector(255,255,255))
{
}

PolyMesh::PolyMesh(char *file, Transform *transform, Vector colour) : Object(colour){
    this->do_construct(file, transform);

    //post construct, form bounding sphere
    Vertex center;
    float radius;
    find_bounding_sphere_values(Vertex(), radius);
    bounding_sphere = new Sphere(center, radius);
}

void PolyMesh::do_construct(char *file, Transform *transform)
{
    ifstream  f;
    f.open(file, ios::in);
    string line;

    //kcply
    getline(f, line);

    //vertex
    getline(f, line);
    vector<string> vertexLine = this->split_string(line);
    vertex_count = stoi(vertexLine[2]);

    //face
    getline(f, line);
    vector<string> faceLine = this->split_string(line);
    triangle_count = stoi(faceLine[2]);

    Vertex *vertices = new Vertex[PolyMesh::vertex_count];
    for(int i = 0; i < vertex_count; i++) {
        getline(f, line);
        vector<string> temp = this->split_string(line);
        vertices[i] = *new Vertex(
                stof(temp[0]),
                stof(temp[1]),
                stof(temp[2])
        );
        transform->apply(vertices[i]);
    }
    PolyMesh::vertex = vertices;

    TriangleIndex *triangles = new TriangleIndex[PolyMesh::triangle_count];
    for(int i = 0; i < triangle_count; i++) {
        getline(f, line);
        vector<string> temp = this->split_string(line);
        triangles[i][0] = stoi(temp[1])-1;
        triangles[i][1] = stoi(temp[2])-1;
        triangles[i][2] = stoi(temp[3])-1;
    }
    PolyMesh::triangle = triangles;
}

vector<string> PolyMesh::split_string(string line) {
    istringstream ss(line);
    istream_iterator<string> begin(ss), end;
    vector<string> arr(begin, end);
    return arr;
}

bool PolyMesh::find_bounding_sphere_values(Vertex centre, float radius) {
    //Implementation of Jack Ritter's Bounding Sphere algorithm.

    //assume any point as min and max
    Vertex min_limit = vertex[0];
    Vertex max_limit = vertex[0];

    //find limits
    for(int i = 0; i < vertex_count; i++) {
        if(vertex[i].x > min_limit.x) min_limit.x = vertex[i].x;
        if(vertex[i].x < max_limit.x) max_limit.x = vertex[i].x;
        if(vertex[i].y > min_limit.y) min_limit.y = vertex[i].y;
        if(vertex[i].y < max_limit.y) max_limit.y = vertex[i].y;
        if(vertex[i].z > min_limit.z) min_limit.z = vertex[i].z;
        if(vertex[i].z < max_limit.z) max_limit.z = vertex[i].z;
    }

    float dx = abs(max_limit.x - min_limit.x);
    float dy = abs(max_limit.y - min_limit.y);
    float dz = abs(max_limit.z - min_limit.z);
    
    //form initial sphere
    radius = max({dx, dy, dz}) / 2;
    centre = (min_limit + max_limit) / 2;

    for(int i = 0; i < vertex_count; i++) {
        Vector point_vector = vertex[i] - centre;
        float distance = point_vector.magnitude();

        //test if point outside sphere
        if(distance > radius) {
            float difference = (distance - radius) / 2;
            radius = radius + difference;
            centre = centre + difference * point_vector;
        }
    }
}

void PolyMesh::intersection(Ray ray, Hit &hit) {
    //Before testing polymesh intersection, check bounding sphere;
    Hit bounding_hit = Hit();
    bounding_sphere->intersection(ray, bounding_hit);
    if(!bounding_hit.flag) return;


    hit.flag = false;
    hit.t = MAXFLOAT;
    float epsilon = 0.0000001;
    //for each triangle, check the intersections
    for(int i = 0; i < triangle_count; i++) {
        //Moller Trombore Algorithm for triangle intersection
        // barycentric:     p = wA + uB + vC
        // ray:             p = o + tD
        // convert to tuv space and solve
        Vertex a = vertex[triangle[i][0]];
        Vertex b = vertex[triangle[i][1]];
        Vertex c = vertex[triangle[i][2]];

        Vector ac = c - a;
        Vector ab = b - a;
        Vector pvec;
        ray.direction.cross(ab, pvec);
        float determinant = ac.dot(pvec);

        //if very close to zero, will not hit the triangle
        if(fabs(determinant) < epsilon)
            continue;

        float inverseDeterminant = 1 / determinant;

        //converting to tuv space
        Vector tvec = Vertex(0,0,0) - a;

        float u = tvec.dot(pvec) * inverseDeterminant;
        if (u < 0 || u > 1)
            continue;

        Vector qvec;
        tvec.cross(ac, qvec);
        float v = ray.direction.dot(qvec) * inverseDeterminant;
        if (v < 0 || u + v > 1)
            continue;

        //thus intersection has been found, so calculate distance
        float t = ab.dot(qvec) * inverseDeterminant;
        if(t < hit.t){
            hit.t = t;
            hit.position.x = ray.position.x + t * ray.direction.x;
            hit.position.y = ray.position.y + t * ray.direction.y;
            hit.position.z = ray.position.z + t * ray.direction.z;
            //triangle normal
            ab.cross(ac, hit.normal);
            hit.normal.normalise();
        }
        hit.what = this;
        hit.flag = true;
    }
}
