#include "clippolygone.h"
#include <iostream>  // header in standard library
#include <cmath>     // header for math functions
#include <random>
#include <string>
#include <algorithm>
#include <list>
#include <stdio.h>
#include <vector>
#include "liblbfgs-master/include/lbfgs.h"

/*
----------------------------VECTOR CLASS-----------------------------------------
*/
Vector::Vector(double x, double y){
    coords[0] = x;
    coords[1] = y;
}

bool &Vector::operator==(const Vector &b){
    bool cond1 = coords[0] == b[0];
    bool cond2 = coords[1] == b[1];
    static bool res = cond1 && cond2;
    return res;
}

bool &Vector::operator!=(const Vector &b){
    bool cond1 = coords[0] != b[0];
    bool cond2 = coords[1] != b[1];
    static bool res = cond1 || cond2;
    return res;
}

Vector &Vector::operator+=(const Vector &b){
    coords[0] += b[0];
    coords[1] += b[1];
    return *this;
}

Vector &Vector::operator*=(const Vector &b){
    coords[0] *= b[0];
    coords[1] *= b[1];
    return *this;
}

Vector &Vector::operator/=(const Vector &b){
    coords[0] /= b[0];
    coords[1] /= b[1];
    return *this;
}

Vector &Vector::operator-=(const Vector &b)
{
    coords[0] -= b[0];
    coords[1] -= b[1];
    return *this;
}

const double &Vector::operator[](int i) const { return coords[i]; }
double &Vector::operator[](int i) { return coords[i]; }

Vector Vector::operator+(const Vector &a){
    return Vector(a[0] + coords[0], a[1] + coords[1]);
}

Vector Vector::operator+(const double a){
    return Vector(a + coords[0], a + coords[1]);
}

Vector Vector::operator-(const Vector &a){
    return Vector(coords[0] - a[0], coords[1] - a[1]);
}

Vector Vector::operator-(const double a){
    return Vector(coords[0] - a, coords[1] - a);
}

Vector Vector::operator*(const Vector &a){
    return Vector(a[0] * coords[0], a[1] * coords[1]);
}

Vector Vector::operator*(const double a){
    return Vector(a * coords[0], a * coords[1]);
}

Vector Vector::operator/(const Vector &a){
    return Vector(coords[0] / a[0], coords[1] / a[1]);
}

Vector Vector::operator/(const double a){
    return Vector(coords[0] / a, coords[1] / a);
}

double Vector::dot(const Vector &a){
    return a[0] * coords[0] + a[1] * coords[1];
}


/*
----------------------------EDGE CLASS-----------------------------------------
*/
Edge::Edge(const Vector &a, const Vector &b){
        this->point_a = a;
    this->point_b = b;
}


/*
----------------------------POLYGON CLASS-----------------------------------------
*/
Polygon::Polygon(){
}

/*
---------------------------- HELPERS -----------------------------------------
*/
Vector intersect(Vector &prevVertex, Vector &curVertex, Edge &clipEdge){
    Vector N = Vector(clipEdge.point_b[1] - clipEdge.point_a[1], clipEdge.point_a[0] - clipEdge.point_b[0]);
    Vector ab = curVertex - prevVertex;
    double t = N.dot(clipEdge.point_a - prevVertex) / N.dot(ab);
    Vector P = prevVertex + ab * t;
    if (t < 0 || t > 1)
        return Vector(0., 0.);
    return P;
}


bool is_inside(Vector &vertex, Edge &clipEdge)
{
    Vector N = Vector(clipEdge.point_b[1] - clipEdge.point_a[1], clipEdge.point_a[0] - clipEdge.point_b[0]) * (-1);
    bool test = N.dot(vertex - clipEdge.point_a) <= 0;
    return test;
}


/*
---------------------------- ALGORITHM CLIPPOLYGON -----------------------------------------
Sutherland-Hodgman by p.84 of lecture notes
*/
Polygon clipPolygonfunction(Polygon &subjectPolygon, Polygon &clipPolygone)
{
    int prevIndex;
    for (int i = 0; i < clipPolygone.edges.size(); i++)
    {
        Edge clipEdge = clipPolygone.edges[i];
        Polygon outPolygone = Polygon();
        for (int j = 0; j < subjectPolygon.vertices.size(); j++)
        {
            Vector curVertex = subjectPolygon.vertices[j];
            if (j > 0)
                prevIndex = j - 1;
            else
                prevIndex = subjectPolygon.vertices.size() - 1;
            Vector prevVertex = subjectPolygon.vertices[prevIndex];
            Vector intersection = intersect(prevVertex, curVertex, clipEdge);
            if (is_inside(curVertex, clipEdge))
            {
                if (!is_inside(prevVertex, clipEdge))
                {
                    outPolygone.vertices.push_back(intersection);
                }
                outPolygone.vertices.push_back(curVertex);
            }
            else if (is_inside(prevVertex, clipEdge))
            {
                outPolygone.vertices.push_back(intersection);
            }
        }

        subjectPolygon = outPolygone;
    }
    return subjectPolygon;
}

 
// saves a static svg file. The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
void save_svg(const std::vector<Polygon> &polygons, std::string filename, std::string fillcol = "none") {
        FILE* f = fopen(filename.c_str(), "w+"); 
        fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
        for (int i=0; i<polygons.size(); i++) {
            fprintf(f, "<g>\n");
            fprintf(f, "<polygon points = \""); 
            for (int j = 0; j < polygons[i].vertices.size(); j++) {
                fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
            }
            fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
            fprintf(f, "</g>\n");
        }
        fprintf(f, "</svg>\n");
        fclose(f);
}