/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */
#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include "vector.h"
#include "polymesh.h"
#include "sphere.h"
#include "utils.h"

using namespace std;

// Class to model objects constructed from triangles

PolyMesh::PolyMesh(char *file, Transform *transform) : PolyMesh(file, transform, new Phong()) {}

PolyMesh::PolyMesh(char *file, Transform *transform, Material *material) : Object(material) {
    this->do_construct(file, transform);

    //post construct, form bounding sphere
    Vertex c;
    float r;
    find_bounding_sphere_values(c, r);
    bounding_sphere = new Sphere(c, r);
    centre = c;
    radius = r;
}

void PolyMesh::do_construct(char *file, Transform *transform) {
    ifstream  f;
    f.open(file, ios::in);
    string line;

    //kcply
    getline(f, line);

    //vertex
    getline(f, line);
    vector<string> vertexLine = Utils::split_string(line);
    vertex_count = stoi(vertexLine[2]);

    //face
    getline(f, line);
    vector<string> faceLine = Utils::split_string(line);
    triangle_count = stoi(faceLine[2]);

    vertices = new Vertex[vertex_count];
    triangles = new TriangleIndex[triangle_count];
    face_normal = new Vector[triangle_count];
    vertex_normal = new Vector[vertex_count];

    // build vertices
    for(int i = 0; i < vertex_count; i++) {
        getline(f, line);
        vector<string> temp = Utils::split_string(line);
        vertices[i] = *new Vertex(
                stof(temp[0]),
                stof(temp[1]),
                stof(temp[2])
        );
        transform->apply(vertices[i]);
    }

    // build triangles and normals
    for(int i = 0; i < triangle_count; i++) {
        getline(f, line);
        vector<string> temp = Utils::split_string(line);
        triangles[i][0] = stoi(temp[1]);
        triangles[i][1] = stoi(temp[2]);
        triangles[i][2] = stoi(temp[3]);
        compute_face_normal(i, face_normal[i]);
    }

    compute_vertex_normals();
}

// Implementation of Jack Ritter's Bounding Sphere algorithm.
void PolyMesh::find_bounding_sphere_values(Vertex &c, float &r) {
    //find limits
    for(int i = 0; i < vertex_count; i++) {
        if(vertices[i].x < min_limit_x) min_limit_x = vertices[i].x;
        if(vertices[i].x > max_limit_x) max_limit_x = vertices[i].x;
        if(vertices[i].y < min_limit_y) min_limit_y = vertices[i].y;
        if(vertices[i].y > max_limit_y) max_limit_y = vertices[i].y;
        if(vertices[i].z < min_limit_z) min_limit_z = vertices[i].z;
        if(vertices[i].z > max_limit_z) max_limit_z = vertices[i].z;
    }

    float dx = abs(max_limit_x - min_limit_x);
    float dy = abs(max_limit_y - min_limit_y);
    float dz = abs(max_limit_z - min_limit_z);

    // max points used in texture generation
    Vertex max_point = {static_cast<float>(max_limit_x), static_cast<float>(max_limit_y), static_cast<float>(max_limit_z)};
    Vertex min_point = {static_cast<float>(min_limit_x), static_cast<float>(min_limit_y), static_cast<float>(min_limit_z)};

    //form initial sphere
    r = max({dx, dy, dz}) / 2;
    c = (max_point + min_point) / 2;

    for(int i = 0; i < vertex_count; i++) {
        Vector point_vector = vertices[i] - c;
        float distance = point_vector.magnitude();

        //test if point outside sphere
        if (distance > r) {
            float difference = (distance - r) / 2;
            r = r + difference;
            c = c + difference * point_vector;
        }
    }
}

void PolyMesh::intersection(Ray ray, Hit &hit) {
    //Before testing polymesh intersection, check bounding sphere;
    Hit bounding_hit = Hit();
    bounding_sphere->intersection(ray, bounding_hit);
    if (!bounding_hit.flag) return;

    hit.flag = false;
    float epsilon = 0.00000001;
    //for each triangle, check the intersections
    for(int i = 0; i < triangle_count; i++) {
        //Moller Trombore Algorithm for triangle intersection
        // barycentric:     p = wA + uB + vC
        // ray:             p = o + tD
        // convert to tuv space and solve
        Vertex a = vertices[triangles[i][0]];
        Vertex b = vertices[triangles[i][1]];
        Vertex c = vertices[triangles[i][2]];

        Vector ac = c - a;
        Vector ab = b - a;
        Vector pvec;
        ray.direction.cross(ab, pvec);
        float determinant = ac.dot(pvec);

        //if very close to zero, will not hit the triangle
        if(fabs(determinant) < epsilon)
            continue;

        double inverseDeterminant = 1.0 / determinant;

        //converting to tuv space
        Vector tvec = ray.position - a;
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
        if(t < hit.t && t > 0){
            hit.t = t;
            hit.position = ray.position + t * ray.direction;

            float u,v,w;
            // interpolate normals based on barycentric coordinates
            compute_barycentric(hit.position, a, b, c, u, v, w);
            hit.normal = vertex_normal[triangles[i][0]] * u +
                         vertex_normal[triangles[i][1]] * v +
                         vertex_normal[triangles[i][2]] * w;
            hit.normal.negate();
            hit.normal.normalise();

            hit.flag = true;
            hit.what = this;
        }
    }
}

void PolyMesh::compute_vertex_normals() {
    // collect neighbouring normals on faces and normalise to find vertex normal
    for (int i = 0; i < triangle_count; i++) {
        for (int j = 0; j < 3; j++) {
            vertex_normal[triangles[i][j]] = vertex_normal[triangles[i][j]] + face_normal[i];
        }
    }

    for (int i = 0; i < vertex_count; i++) {
        vertex_normal[i].normalise();
    }
}

void PolyMesh::compute_face_normal(int triangle, Vector &normal) {
    // cross two triangle edges to find normals
    Vector u = vertices[triangles[triangle][1]] - vertices[triangles[triangle][0]];
    u.normalise();

    Vector v = vertices[triangles[triangle][2]] - vertices[triangles[triangle][0]];
    v.normalise();

    u.cross(v, normal);
    normal.normalise();
}

// Transcribed from Real-Time Collision Detection by Christer Ericson
void PolyMesh::compute_barycentric(Vertex p, Vertex a, Vertex b, Vertex c, float &u, float &v, float &w) {
    Vector v0 = b - a, v1 = c - a, v2 = p - a;
    float d00 = v0.dot(v0);
    float d01 = v0.dot(v1);
    float d11 = v1.dot(v1);
    float d20 = v2.dot(v0);
    float d21 = v2.dot(v1);
    float denom = d00 * d11 - d01 * d01;
    v = (d11 * d20 - d01 * d21) / denom;
    w = (d00 * d21 - d01 * d20) / denom;
    u = 1.0f - v - w;
}