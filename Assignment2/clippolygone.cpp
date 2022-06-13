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

double Vector::norm(){
    return sqrt(dot(*this));
}

Vector Vector::normalize(){
    return *this / norm();
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
Vector intersect(Vector &prevVertex, Vector &tempVector, Edge &clipEdge){
    Vector N = Vector(clipEdge.point_b[1] - clipEdge.point_a[1], clipEdge.point_a[0] - clipEdge.point_b[0]);
    Vector diff = tempVector - prevVertex;
    double t = N.dot(clipEdge.point_a - prevVertex) / N.dot(diff);
    Vector P = prevVertex + diff * t;
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
    int formerIndex;
    for (int i = 0; i < clipPolygone.edges.size(); i++)
    {
        Edge clipEdge = clipPolygone.edges[i];
        Polygon outPolygone = Polygon();
        for (int j = 0; j < subjectPolygon.vertices.size(); j++)
        {
            Vector tempVector = subjectPolygon.vertices[j];
            if (j > 0)
                formerIndex = j - 1;
            else
                formerIndex = subjectPolygon.vertices.size() - 1;
            Vector prevVertex = subjectPolygon.vertices[formerIndex];
            Vector intersection = intersect(prevVertex, tempVector, clipEdge);
            if (is_inside(tempVector, clipEdge))
            {
                if (!is_inside(prevVertex, clipEdge))
                {
                    outPolygone.vertices.push_back(intersection);
                }
                outPolygone.vertices.push_back(tempVector);
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


/*
---------------------------- VORONOI DIAGRAM -----------------------------------------
*/

Polygon clipPolygonLine(Polygon &subjectPolygon, Vector M, Vector vectorIJ)
{
    int formerIndex;
    Polygon resultPoly = Polygon();
    Vector vectorJI = vectorIJ * (-1);
#define inside(X) (X - M).dot(vectorIJ)

    for (int j = 0; j < subjectPolygon.vertices.size(); j++)
    {
        if (j<=0){
            formerIndex = subjectPolygon.vertices.size() - 1;
        }
        else{
            formerIndex = j - 1;
        }
        Vector tempVector = subjectPolygon.vertices[j];
        Vector prevVertex = subjectPolygon.vertices[formerIndex];
        Vector diff = tempVector - prevVertex;
        double t = vectorJI.dot(M - prevVertex) / vectorJI.dot(diff);
        Vector intersection = (t >= 0 && t <= 1) ? prevVertex + diff * t : Vector(0., 0.);
        if (inside(tempVector) < 0)
        {
            if (!(inside(prevVertex) < 0))
            {
                resultPoly.vertices.push_back(intersection);
            }
            resultPoly.vertices.push_back(tempVector);
        }
        else if (inside(prevVertex) < 0)
        {

            resultPoly.vertices.push_back(intersection);
        }
    }
    return resultPoly;
}

std::vector<Polygon> voronoi(Polygon &clipPolygon, std::vector<Vector> &points, const double *weights)
{
    Vector firstPoint, secondPoint, middlePoint, vectorIJ;
    double firstWeight, secondWeight;
    std::vector<Polygon> res;
    for (int i = 0; i < points.size(); i++)
    {
        firstWeight = weights[i];
        firstPoint = points[i];
        Polygon resultPoly = clipPolygon;
        for (int j = 0; j < points.size(); j++)
        {
            secondPoint = points[j];
            secondWeight = weights[j];
            if (i == j){
                continue;
            }
            middlePoint = (firstPoint + secondPoint) * 0.5;
            vectorIJ = (secondPoint - firstPoint);
            middlePoint = middlePoint + vectorIJ * (firstWeight - secondWeight) / (2 * pow(vectorIJ.norm(), 2.));
            resultPoly = clipPolygonLine(resultPoly, middlePoint, vectorIJ);
        }
        res.push_back(resultPoly);
    }
    return res;
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