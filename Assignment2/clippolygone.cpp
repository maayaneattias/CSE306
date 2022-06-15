#include "clippolygone.h"
#include <iostream>  // header in standard library
#include <cmath>     // header for math functions
#include <random>
#include <string>
#include <algorithm>
#include <list>
#include <stdio.h>
#include <vector>


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

double Polygon::area()
{
    double result = 0;
    int n = vertices.size();
    if (n == 0)
    {
        return 0;
    }
    for (int i = 0; i < n; i++)
    {
        Vector point1 = vertices[i];
        Vector point2 = (i < n - 1) ? vertices[i + 1] : vertices[0];
        result += point1[0] * point2[1] - point1[1] * point2[0];
    }
    return std::abs(0.5 * result);
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

double polygon_area(std::vector<Vector> vertices)
{
    double result = 0;
    int n = vertices.size();
    if (n == 0)
    {
        return 0;
    }
    for (int i = 0; i < n; i++)
    {
        Vector point1 = vertices[i];
        Vector point2 = (i < n - 1) ? vertices[i + 1] : vertices[0];
        result += point1[0] * point2[1] - point1[1] * point2[0];
    }
    return std::abs(0.5 * result);
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

Polygon clipPolygonBis(Polygon &subjectPolygon, Vector M, Vector vectorIJ)
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
            resultPoly = clipPolygonBis(resultPoly, middlePoint, vectorIJ);
        }
        res.push_back(resultPoly);
    }
    return res;
}

/*
----------------------------  FOR OPTIMAL TRANSPORT -----------------------------------------
*/
double triangulate(std::vector<Vector> vertices, Vector point)
{
    if (vertices.size() == 0){
        return 0;
    }
    Vector initial_vertex = vertices[0];
    double res = 0;
    for (int i = 0; i < vertices.size() - 2; i++)
    {
        Vector vertex1 = vertices[i + 1];
        Vector vertex2 = vertices[i + 2];
        std::vector<Vector> triangle;
        triangle.push_back(initial_vertex);
        triangle.push_back(vertex1);
        triangle.push_back(vertex2);
        double area = std::abs(polygon_area(triangle));
        res += (area / 6.) * ((initial_vertex - point).dot(initial_vertex - point) + (initial_vertex - point).dot(vertex1 - point) + (initial_vertex - point).dot(vertex2 - point) +
                              (vertex1 - point).dot(vertex1 - point) + (vertex1 - point).dot(vertex2 - point) +
                              (vertex2 - point).dot(vertex2 - point));
    }
    return std::abs(res);
}

void gradient_descent(Polygon &clipPolygon, std::vector<Vector> &points, const double *lambdas, double eps, double step, double *weights)
{
    double fx, normOfGrad;
    double grad[points.size()], area[points.size()];
    std::vector<Polygon> polygons = voronoi(clipPolygon, points, lambdas);
    save_svg(polygons, "beforeOptimisation.svg");
    double error = 1;
    int k = 1;
    while (error > eps)
    {
        error = 0;
        fx = 0;
        normOfGrad = 0;
        polygons = voronoi(clipPolygon, points, weights);

        for (int i = 0; i < points.size(); i++)
        {
            area[i] = polygons[i].area();
            grad[i] = (lambdas[i] - area[i]);
            normOfGrad += grad[i] * grad[i];
        }
        for (int i = 0; i < points.size(); i++)
        {
            Vector point;
            std::vector<Vector> vertices;
            double temp;
            point = points[i];
            vertices = polygons[i].vertices;
            temp = (step / sqrt(normOfGrad)) * grad[i];
            error += temp * temp;
            weights[i] += temp;
            fx += (triangulate(vertices, point) - weights[i] * area[i] + lambdas[i] * weights[i]);
        }
        error = sqrt(error);
        k += 1;
        if (k % 50 == 0){
            std::string filename = "opti_" + std::to_string(k) + ".svg";
            save_svg(polygons, filename);
        }
    }
}

 
/*
----------------------------OPTIMAL TRANSPORT CLASS-----------------------------------------
*/
OptimalTransport::OptimalTransport(Polygon &clipPolygon, std::vector<Vector> &points, int maximumterations){
    this->clipPolygon = clipPolygon;
    this->points = points;
    this->max_iterations = maximumterations;
}

OptimalTransport::~OptimalTransport(){
    if (m_x != NULL)
        {
            lbfgs_free(m_x);
            m_x = NULL;
        }
}

int OptimalTransport::execute(){
        lbfgsfloatval_t fx;
        int N = this->points.size();
        this->lambdas = (double *)malloc((N) * sizeof(double));

        lbfgsfloatval_t *m_x = lbfgs_malloc(N);
        if (m_x == NULL)
        {
            printf("ERROR: No memory allocated\n");
            return 1;
        }

        Vector C = Vector(0.5, 0.5);
        Vector diff;
        double tot = 0;
        double maximum = 0;
        for (int i = 0; i < N; i++)
        {
            lambdas[i] = std::exp(-pow(diff.norm(), 2.) / 0.02);
            tot += lambdas[i];
            diff = C - this->points[i];
        }
        for (int i = 0; i < N; i++)
        {
            m_x[i] = 1;
            lambdas[i] /= tot;
            if (lambdas[i] > maximum){
                maximum = lambdas[i];
            }
        }
        std::cout << "dirac:" << maximum << std::endl;
        this->polygons = voronoi(clipPolygon, this->points, lambdas);
        save_svg(polygons, "beforeOptimisation.svg");
        lbfgs_parameter_t param;
        lbfgs_parameter_init(&param);
        param.max_iterations = this->max_iterations;

        int ret = lbfgs(N, m_x, &fx, _evaluate, _progress, this, &param);
        this->polygons = voronoi(clipPolygon, this->points, m_x);

        printf("L-BFGS optimization terminated. Status code: %d\n", ret);
        free(lambdas);
        return ret;
}

lbfgsfloatval_t OptimalTransport::_evaluate(void *instance,const lbfgsfloatval_t *x,lbfgsfloatval_t *g,
                                    const int n, const lbfgsfloatval_t step){
        return reinterpret_cast<OptimalTransport *>(instance)->evaluate(x, g, n, step);
}

lbfgsfloatval_t OptimalTransport::evaluate(const lbfgsfloatval_t *x,lbfgsfloatval_t *g,const int n,const lbfgsfloatval_t step){
        lbfgsfloatval_t fx = 0.0;
        this->polygons = voronoi(clipPolygon, this->points, x);
        for (int i = 0; i < n; i++){
            std::vector<Vector> vertices = this->polygons[i].vertices;
            Vector point = this->points[i];
            double area = this->polygons[i].area();
            double temp = triangulate(vertices, point);
            fx += temp - x[i] * area + this->lambdas[i] * x[i];
            g[i] = area - this->lambdas[i];
        }
        fx = fx * -1.;
        return fx;
}

int OptimalTransport::_progress(void *instance,const lbfgsfloatval_t *x,const lbfgsfloatval_t *g,const lbfgsfloatval_t fx,
                        const lbfgsfloatval_t xnorm,const lbfgsfloatval_t gnorm,const lbfgsfloatval_t step,
                        int n,int k,int ls){
        return reinterpret_cast<OptimalTransport *>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);

}

int OptimalTransport::progress(const lbfgsfloatval_t *x,const lbfgsfloatval_t *g,const lbfgsfloatval_t fx,const lbfgsfloatval_t xnorm,
                const lbfgsfloatval_t gnorm,const lbfgsfloatval_t step,int n,int k,int ls){
        printf("Number of iteration is %d\n", k);
        printf(" Progress: fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
        printf("\n");
        return 0;   
}
