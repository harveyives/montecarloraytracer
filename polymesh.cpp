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
#include <vector>

#include "polymesh.h"

using namespace std;

PolyMesh::PolyMesh(char *file)
{
  Transform *transform = new Transform();

  this->do_construct(file, transform);
}

PolyMesh::PolyMesh(char *file, Transform *transform)
{
  this->do_construct(file, transform);
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

    //putting all the tokens in the vector
    vector<string> array(begin, end);
    return array;
}