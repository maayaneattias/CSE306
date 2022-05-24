#include "auxiliary.h" 
#include <iostream> 
#include <cmath>   
#include <random>
#include <algorithm>
#include <list>
#include <stdio.h>
#include <string.h>

static std::default_random_engine engine(10); // random seed = 10
static std::uniform_real_distribution<double> uniform(0., 1.);

/*-----------------------------------------------------------------------------------------------*/

/*                          Vector class                                                         */

/*-----------------------------------------------------------------------------------------------*/
Vector::Vector(double x, double y, double z)
{
    coords[0] = x;
    coords[1] = y;
    coords[2] = z;
}

bool &Vector::operator==(const Vector &b)
{
    bool cond1 = coords[0] == b[0];
    bool cond2 = coords[1] == b[1];
    bool cond3 = coords[2] == b[2];
    static bool res = cond1 && cond2 && cond3;
    return res;
}

bool &Vector::operator!=(const Vector &b)
{
    bool cond1 = coords[0] != b[0];
    bool cond2 = coords[1] != b[1];
    bool cond3 = coords[2] != b[2];
    static bool res = cond1 || cond2 || cond3;
    return res;
}

Vector &Vector::operator+=(const Vector &b)
{
    coords[0] += b[0];
    coords[1] += b[1];
    coords[2] += b[2];
    return *this;
}

Vector &Vector::operator*=(const Vector &b)
{
    coords[0] *= b[0];
    coords[1] *= b[1];
    coords[2] *= b[2];
    return *this;
}

Vector &Vector::operator/=(const Vector &b)
{
    coords[0] /= b[0];
    coords[1] /= b[1];
    coords[2] /= b[2];
    return *this;
}

Vector &Vector::operator-=(const Vector &b)
{
    coords[0] -= b[0];
    coords[1] -= b[1];
    coords[2] -= b[2];
    return *this;
}

const double &Vector::operator[](int i) const { return coords[i]; }
double &Vector::operator[](int i) { return coords[i]; }

Vector Vector::operator+(const Vector &a)
{
    return Vector(a[0] + coords[0], a[1] + coords[1], a[2] + coords[2]);
}

Vector Vector::operator+(const double a)
{
    return Vector(a + coords[0], a + coords[1], a + coords[2]);
}

Vector Vector::operator-(const Vector &a)
{
    return Vector(coords[0] - a[0], coords[1] - a[1], coords[2] - a[2]);
}

Vector Vector::operator-(const double a)
{
    return Vector(coords[0] - a, coords[1] - a, coords[2] - a);
}

Vector Vector::operator*(const Vector &a)
{
    return Vector(a[0] * coords[0], a[1] * coords[1], a[2] * coords[2]);
}

Vector Vector::operator*(const double a)
{
    return Vector(a * coords[0], a * coords[1], a * coords[2]);
}

Vector Vector::operator/(const Vector &a)
{
    return Vector(coords[0] / a[0], coords[1] / a[1], coords[2] / a[2]);
}

Vector Vector::operator/(const double a)
{
    return Vector(coords[0] / a, coords[1] / a, coords[2] / a);
}

double Vector::dot(const Vector &a)
{
    return a[0] * coords[0] + a[1] * coords[1] + a[2] * coords[2];
}

Vector Vector::min(double s)
{
    return Vector(std::min(coords[0], s), std::min(coords[1], s), std::min(coords[2], s));
}

Vector Vector::max(double s)
{
    return Vector(std::max(coords[0], s), std::max(coords[1], s), std::max(coords[2], s));
}

Vector Vector::pow(double s)
{
    return Vector(std::pow(coords[0], s), std::pow(coords[1], s), std::pow(coords[2], s));
}

Vector Vector::cross_product(const Vector &a)
{

    return Vector(coords[1] * a[2] - coords[2] * a[1], coords[2] * a[0] - coords[0] * a[2], coords[0] * a[1] - coords[1] * a[0]);
}



int Vector::argmin()
{
    int res = 0;
    double mini = abs(coords[0]);
    int i = 1;
    while (i < 3)
    {
        if (abs(coords[i]) < mini)
        {
            res = i;
            mini = abs(coords[i]);
        }
        i += 1;
    }
    return res;
}

int Vector::get_longest()
{
    int res = 0;
    double maxi = coords[0];
    int i = 1;
    while (i < 3)
    {
        if (coords[i] > maxi)
        {
            res = i;
            maxi = coords[i];
        }
        i += 1;
    }
    return res;
}

double Vector::norm()
{
    return sqrt(dot(*this));
}

Vector Vector::normalize()
{
    return *this / norm();
}

/*-----------------------------------------------------------------------------------------------*/

/*                          Sphere class                                                         */

/*-----------------------------------------------------------------------------------------------*/
Sphere::Sphere(Vector C, Vector albedo, double R, bool mirror, bool transparent, double refraction_index)
{
    this->C = C;
    this->albedo = albedo;
    this->R = R;
    this->mirror = mirror;
    this->transparent = transparent;
    this->refraction_index = refraction_index;
}

Intersection Sphere::intersect(Ray &r)
{
    Intersection intersection;
    Vector c = r.O - this->C;
    double dotproductUC = r.u.dot(c);
    double norm_c = c.norm();
    double delta = dotproductUC * dotproductUC - norm_c * norm_c + R * R;
    if (delta < 0)
    {
        intersection.is_intersect = false;
        return intersection;
    }
    else
    {
        double temp = r.u.dot(c) * (-1) + sqrt(delta);
        if (temp < 0)
        {
            intersection.is_intersect = false;
        }
        else
        {
            intersection.is_intersect = true;
            double tempbis = temp - 2 * sqrt(delta);
            if (tempbis >= 0)
            {
                intersection.t = tempbis;
            }
            else
            {
                intersection.t = temp;
            };
        };
    };
    intersection.P = r.O + r.u * intersection.t;
    Vector C_P = intersection.P - this->C;
    intersection.N = C_P.normalize();
    return intersection;
}


/*-----------------------------------------------------------------------------------------------*/

/*                          Ray class                                                         */

/*-----------------------------------------------------------------------------------------------*/
Ray::Ray(Vector O, Vector u, double refraction_index)
{
    this->O = O;
    this->u = u;
    this->refraction_index = refraction_index;
}


/*-----------------------------------------------------------------------------------------------*/

/*                          HELPERS                                                            */

/*-----------------------------------------------------------------------------------------------*/

Vector swap(const Vector &N, int first, int second)
{
    Vector res = Vector(0, 0, 0);
    res[first] = N[second];
    res[second] = -N[first];
    return res;
}

Vector randomize_cross_product(Vector &N)
{
    //to modify
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    double x = cos(2 * M_PI * r1) * sqrt(1 - r2);
    double y = sin(2 * M_PI * r1) * sqrt(1 - r2);
    double z = sqrt(r2);

    Vector T1 = Vector(0, 0, 0);
    int min_idx = N.argmin();
    T1[min_idx] = 0;
    if (min_idx == 0)
        T1 = swap(N, 1, 2);
    if (min_idx == 1)
        T1 = swap(N, 0, 2);
    if (min_idx == 2)
        T1 = swap(N, 0, 1);
    T1 = T1.normalize();
    Vector T2 = N.cross_product(T1);
    return (T1 * x + T2 * y + N * z).normalize();
}

void boxMuller(double stdev, double &x, double &y)
{
    double r1 = uniform(engine);
    double r2 = uniform(engine);

    x = sqrt(-2 * log(r1)) * cos(2 * M_PI * r2) * stdev;
    y = sqrt(-2 * log(r1)) * sin(2 * M_PI * r2) * stdev;
}

Ray ray_centering(Vector &camera_pos, int W, int H, double fov, double refraction_index, int i, int j, bool random)
{   
    Vector r;
    r[1] = camera_pos[1] - i - 0.5 + H / 2;
    r[2] = camera_pos[2] - W / (2 * tan(fov / 2));
    r[0] = camera_pos[0] + j + 0.5 - W / 2;
    if (random){
        double x, y;
        boxMuller(0.5, x, y);
        r[0] += x;
        r[1] += y;
    }
    return Ray(camera_pos, (r - camera_pos).normalize(), refraction_index);
}


/*-----------------------------------------------------------------------------------------------*/

/*                          Bbox                                                                 */

/*-----------------------------------------------------------------------------------------------*/

Bbox::Bbox()
{
    this->b_min = Vector(0, 0, 0);
    this->b_max = Vector(0, 0, 0);
}
Bbox::Bbox(Vector &b_min, Vector &b_max)
{
    this->b_min = b_min;
    this->b_max = b_max;
}

bool Bbox::intersect(Ray &r, double &distance_intersection)
{

    double txmin, tymin, tzmin;
    double txmax, tymax, tzmax;
    double temp, tempbis;

    temp = (b_min[1] - r.O[1]) / r.u[1];
    tempbis = (b_max[1] - r.O[1]) / r.u[1];
    tymin = std::min(temp, tempbis);
    tymax = std::max(temp, tempbis);
    temp = (b_min[2] - r.O[2]) / r.u[2];
    tempbis = (b_max[2] - r.O[2]) / r.u[2];
    tzmin = std::min(temp, tempbis);
    tzmax = std::max(temp, tempbis);
    temp = (b_min[0] - r.O[0]) / r.u[0];
    tempbis = (b_max[0] - r.O[0]) / r.u[0];
    txmin = std::min(temp, tempbis);
    txmax = std::max(temp, tempbis);

    double max = std::max(txmin, std::max(tymin, tzmin));
    double min = std::min(txmax, std::min(tymax, tzmax));
    if (min > max && max > 0)
    {
        distance_intersection = max;
        return true;
    }
    return false;
}

Vector Bbox::compute_diag()
{
    return b_max - b_min;
}

/*-----------------------------------------------------------------------------------------------*/

/*                          Triangle Mesh                                                        */

/*-----------------------------------------------------------------------------------------------*/
TriangleMesh::TriangleMesh(Vector albedo, double scaling_factor, Vector translation)
{
    this->albedo = albedo;
    this->scaling_factor = scaling_factor;
    this->translation = translation;
    this->root = new Node();
    this->mirror = false;
    this->transparent = false;
}

Intersection TriangleMesh::intersect(Ray &r)
{
    Intersection intersection;
    intersection.is_intersect = false;
    double tempdist;
    if (!root->bbox.intersect(r, tempdist)){
        return intersection;
    }

    std::list<Node *> listofvisits;
    listofvisits.push_front(root);
    double maxdist = std::numeric_limits<double>::max();
    Vector vertexI, vertexJ, C, N;
    Vector elements, elementsbis;
    while (!listofvisits.empty())
    {
        Node *curNode = listofvisits.back();
        listofvisits.pop_back();
        if (!curNode->has_children)
        {
            if (curNode->child_left->bbox.intersect(r, tempdist))
            {
                if (tempdist < maxdist)
                {
                    listofvisits.push_back(curNode->child_left);
                }
            }
            if (curNode->child_right->bbox.intersect(r, tempdist))
            {
                if (tempdist < maxdist)
                {

                    listofvisits.push_back(curNode->child_right);
                }
            }
        }
        else
        {
            for (int i = curNode->starting_triangle; i < curNode->ending_triangle; i++)
            {
                TriangleIndices triangle = indices[i];
                vertexI = vertices[triangle.vtxi];
                vertexJ = vertices[triangle.vtxj];
                C = vertices[triangle.vtxk];
                elements = vertexJ - vertexI;
                elementsbis = C - vertexI;
                N = elements.cross_product(elementsbis);
                tempdist = N.dot(vertexI - r.O) / N.dot(r.u);
                if (0 < tempdist && tempdist < maxdist)
                {
                    double un = 1 / N.dot(r.u);
                    double beta = elementsbis.dot((vertexI - r.O).cross_product(r.u)) * un;
                    double gamma = -elements.dot((vertexI - r.O).cross_product(r.u)) * un;
                    double alpha = 1 - beta - gamma;
                    if (beta > 0 && gamma > 0 && alpha > 0)
                    {
                        maxdist = tempdist;
                        intersection.t = maxdist;
                        intersection.N = N.normalize();
                        intersection.P = vertexI + elements * beta + elementsbis * gamma;
                        intersection.is_intersect = true;
                    }
                }
            }
        }
    }
    return intersection;
}

Bbox TriangleMesh::compute_bbox(int starting_triangle, int ending_triangle)
{
    Bbox bbox;
    double bmax_x, bmax_y, bmax_z = -std::numeric_limits<double>::max();
    double bmin_x, bmin_y, bmin_z = std::numeric_limits<double>::max();
    for (int i = starting_triangle; i < ending_triangle; i++)
    {
        TriangleIndices triangle = indices[i];
        std::vector<Vector> points;
        points.push_back(vertices[triangle.vtxi]);
        points.push_back(vertices[triangle.vtxj]);
        points.push_back(vertices[triangle.vtxk]);
        for (int i = 0; i < 3; i++)
        {
            if (points[i][0] > bmax_x)
            {
                bmax_x = points[i][0];
            }
            if (points[i][0] < bmin_x)
            {
                bmin_x = points[i][0];
            }
            if (points[i][1] > bmax_y)
            {
                bmax_y = points[i][1];
            }
            if (points[i][1] < bmin_y)
            {
                bmin_y = points[i][1];
            }
            if (points[i][2] > bmax_z)
            {
                bmax_z = points[i][2];
            }
            if (points[i][2] < bmin_z)
            {
                bmin_z = points[i][2];
            }
        }
    }
    bbox.b_min = Vector(bmin_x, bmin_y, bmin_z);
    bbox.b_max = Vector(bmax_x, bmax_y, bmax_z);
    return bbox;
}

Vector TriangleMesh::find_barycenter(int triangle_indice)
{
    TriangleIndices triangle = indices[triangle_indice];
    Vector res = (vertices[triangle.vtxi] + vertices[triangle.vtxj] + vertices[triangle.vtxk]) / 3;
    return res;
}

void TriangleMesh::performingBVH(Node *node, int starting_triangle, int ending_triangle)
{
    Bbox bbox = compute_bbox(starting_triangle, ending_triangle); 
    node->bbox = bbox;
    node->starting_triangle = starting_triangle;
    node->ending_triangle = ending_triangle;
    Vector diag = node->bbox.compute_diag();
    Vector middle_diag = node->bbox.b_min + diag * 0.5;
    int axis = diag.get_longest();
    int pivoting = starting_triangle;
    for (int i = starting_triangle; i < ending_triangle; i++)
    {
        Vector barycenter = this->find_barycenter(i);
        if (barycenter[axis] < middle_diag[axis])
        {
            std::swap(indices[i], indices[pivoting]);
            pivoting++;
        }
    }
    if (pivoting <= starting_triangle || pivoting >= ending_triangle - 1 || ending_triangle - starting_triangle < 5)
    {
        node->has_children = true;
        return;
    }

    node->has_children = false;
    node->child_left = new Node();
    node->child_right = new Node();
    performingBVH(node->child_left, starting_triangle, pivoting);
    performingBVH(node->child_right, pivoting, ending_triangle);
}

void TriangleMesh::readOBJ(const char *obj)
{

    //NOTE: I used the implementation in the link given
    char matfile[255];
    char grp[255];

    FILE *f;
    f = fopen(obj, "r");
    int curGroup = -1;
    while (!feof(f))
    {
        char line[255];
        if (!fgets(line, 255, f))
            break;

        std::string linetrim(line);
        linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
        strcpy(line, linetrim.c_str());

        if (line[0] == 'u' && line[1] == 's')
        {
            sscanf(line, "usemtl %[^\n]\n", grp);
            curGroup++;
        }

        if (line[0] == 'v' && line[1] == ' ')
        {
            Vector vec;

            Vector col;
            if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6)
            {
                col[0] = std::min(1., std::max(0., col[0]));
                col[1] = std::min(1., std::max(0., col[1]));
                col[2] = std::min(1., std::max(0., col[2]));

                vertices.push_back(vec);
                vertexcolors.push_back(col);
            }
            else
            {
                sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                vertices.push_back(vec);
            }
        }
        if (line[0] == 'v' && line[1] == 'n')
        {
            Vector vec;
            sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
            normals.push_back(vec);
        }
        if (line[0] == 'v' && line[1] == 't')
        {
            Vector vec;
            sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
            uvs.push_back(vec);
        }
        if (line[0] == 'f')
        {
            TriangleIndices t;
            int i0, i1, i2, i3;
            int j0, j1, j2, j3;
            int k0, k1, k2, k3;
            int nn;
            t.group = curGroup;

            char *consumedline = line + 1;
            int offset;

            nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
            if (nn == 9)
            {
                if (i0 < 0)
                    t.vtxi = vertices.size() + i0;
                else
                    t.vtxi = i0 - 1;
                if (i1 < 0)
                    t.vtxj = vertices.size() + i1;
                else
                    t.vtxj = i1 - 1;
                if (i2 < 0)
                    t.vtxk = vertices.size() + i2;
                else
                    t.vtxk = i2 - 1;
                if (j0 < 0)
                    t.uvi = uvs.size() + j0;
                else
                    t.uvi = j0 - 1;
                if (j1 < 0)
                    t.uvj = uvs.size() + j1;
                else
                    t.uvj = j1 - 1;
                if (j2 < 0)
                    t.uvk = uvs.size() + j2;
                else
                    t.uvk = j2 - 1;
                if (k0 < 0)
                    t.ni = normals.size() + k0;
                else
                    t.ni = k0 - 1;
                if (k1 < 0)
                    t.nj = normals.size() + k1;
                else
                    t.nj = k1 - 1;
                if (k2 < 0)
                    t.nk = normals.size() + k2;
                else
                    t.nk = k2 - 1;
                indices.push_back(t);
            }
            else
            {
                nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                if (nn == 6)
                {
                    if (i0 < 0)
                        t.vtxi = vertices.size() + i0;
                    else
                        t.vtxi = i0 - 1;
                    if (i1 < 0)
                        t.vtxj = vertices.size() + i1;
                    else
                        t.vtxj = i1 - 1;
                    if (i2 < 0)
                        t.vtxk = vertices.size() + i2;
                    else
                        t.vtxk = i2 - 1;
                    if (j0 < 0)
                        t.uvi = uvs.size() + j0;
                    else
                        t.uvi = j0 - 1;
                    if (j1 < 0)
                        t.uvj = uvs.size() + j1;
                    else
                        t.uvj = j1 - 1;
                    if (j2 < 0)
                        t.uvk = uvs.size() + j2;
                    else
                        t.uvk = j2 - 1;
                    indices.push_back(t);
                }
                else
                {
                    nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                    if (nn == 3)
                    {
                        if (i0 < 0)
                            t.vtxi = vertices.size() + i0;
                        else
                            t.vtxi = i0 - 1;
                        if (i1 < 0)
                            t.vtxj = vertices.size() + i1;
                        else
                            t.vtxj = i1 - 1;
                        if (i2 < 0)
                            t.vtxk = vertices.size() + i2;
                        else
                            t.vtxk = i2 - 1;
                        indices.push_back(t);
                    }
                    else
                    {
                        nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                        if (i0 < 0)
                            t.vtxi = vertices.size() + i0;
                        else
                            t.vtxi = i0 - 1;
                        if (i1 < 0)
                            t.vtxj = vertices.size() + i1;
                        else
                            t.vtxj = i1 - 1;
                        if (i2 < 0)
                            t.vtxk = vertices.size() + i2;
                        else
                            t.vtxk = i2 - 1;
                        if (k0 < 0)
                            t.ni = normals.size() + k0;
                        else
                            t.ni = k0 - 1;
                        if (k1 < 0)
                            t.nj = normals.size() + k1;
                        else
                            t.nj = k1 - 1;
                        if (k2 < 0)
                            t.nk = normals.size() + k2;
                        else
                            t.nk = k2 - 1;
                        indices.push_back(t);
                    }
                }
            }

            consumedline = consumedline + offset;

            while (true)
            {
                if (consumedline[0] == '\n')
                    break;
                if (consumedline[0] == '\0')
                    break;
                nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                TriangleIndices t2;
                t2.group = curGroup;
                if (nn == 3)
                {
                    if (i0 < 0)
                        t2.vtxi = vertices.size() + i0;
                    else
                        t2.vtxi = i0 - 1;
                    if (i2 < 0)
                        t2.vtxj = vertices.size() + i2;
                    else
                        t2.vtxj = i2 - 1;
                    if (i3 < 0)
                        t2.vtxk = vertices.size() + i3;
                    else
                        t2.vtxk = i3 - 1;
                    if (j0 < 0)
                        t2.uvi = uvs.size() + j0;
                    else
                        t2.uvi = j0 - 1;
                    if (j2 < 0)
                        t2.uvj = uvs.size() + j2;
                    else
                        t2.uvj = j2 - 1;
                    if (j3 < 0)
                        t2.uvk = uvs.size() + j3;
                    else
                        t2.uvk = j3 - 1;
                    if (k0 < 0)
                        t2.ni = normals.size() + k0;
                    else
                        t2.ni = k0 - 1;
                    if (k2 < 0)
                        t2.nj = normals.size() + k2;
                    else
                        t2.nj = k2 - 1;
                    if (k3 < 0)
                        t2.nk = normals.size() + k3;
                    else
                        t2.nk = k3 - 1;
                    indices.push_back(t2);
                    consumedline = consumedline + offset;
                    i2 = i3;
                    j2 = j3;
                    k2 = k3;
                }
                else
                {
                    nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                    if (nn == 2)
                    {
                        if (i0 < 0)
                            t2.vtxi = vertices.size() + i0;
                        else
                            t2.vtxi = i0 - 1;
                        if (i2 < 0)
                            t2.vtxj = vertices.size() + i2;
                        else
                            t2.vtxj = i2 - 1;
                        if (i3 < 0)
                            t2.vtxk = vertices.size() + i3;
                        else
                            t2.vtxk = i3 - 1;
                        if (j0 < 0)
                            t2.uvi = uvs.size() + j0;
                        else
                            t2.uvi = j0 - 1;
                        if (j2 < 0)
                            t2.uvj = uvs.size() + j2;
                        else
                            t2.uvj = j2 - 1;
                        if (j3 < 0)
                            t2.uvk = uvs.size() + j3;
                        else
                            t2.uvk = j3 - 1;
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        indices.push_back(t2);
                    }
                    else
                    {
                        nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                        if (nn == 2)
                        {
                            if (i0 < 0)
                                t2.vtxi = vertices.size() + i0;
                            else
                                t2.vtxi = i0 - 1;
                            if (i2 < 0)
                                t2.vtxj = vertices.size() + i2;
                            else
                                t2.vtxj = i2 - 1;
                            if (i3 < 0)
                                t2.vtxk = vertices.size() + i3;
                            else
                                t2.vtxk = i3 - 1;
                            if (k0 < 0)
                                t2.ni = normals.size() + k0;
                            else
                                t2.ni = k0 - 1;
                            if (k2 < 0)
                                t2.nj = normals.size() + k2;
                            else
                                t2.nj = k2 - 1;
                            if (k3 < 0)
                                t2.nk = normals.size() + k3;
                            else
                                t2.nk = k3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            k2 = k3;
                            indices.push_back(t2);
                        }
                        else
                        {
                            nn = sscanf(consumedline, "%u%n", &i3, &offset);
                            if (nn == 1)
                            {
                                if (i0 < 0)
                                    t2.vtxi = vertices.size() + i0;
                                else
                                    t2.vtxi = i0 - 1;
                                if (i2 < 0)
                                    t2.vtxj = vertices.size() + i2;
                                else
                                    t2.vtxj = i2 - 1;
                                if (i3 < 0)
                                    t2.vtxk = vertices.size() + i3;
                                else
                                    t2.vtxk = i3 - 1;
                                consumedline = consumedline + offset;
                                i2 = i3;
                                indices.push_back(t2);
                            }
                            else
                            {
                                consumedline = consumedline + 1;
                            }
                        }
                    }
                }
            }
        }
    }
    fclose(f);
}


