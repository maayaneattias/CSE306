#include <iostream>
#include <vector>

class Vector{
public:
    explicit Vector(double x = 0., double y = 0., double z = 0.);
    bool &operator==(const Vector &b);
    bool &operator!=(const Vector &b);
    Vector &operator+=(const Vector &b);
    Vector &operator-=(const Vector &b);
    Vector &operator*=(const Vector &b);
    Vector &operator/=(const Vector &b);
    Vector operator+(const Vector &a);
    Vector operator+(const double a);
    Vector operator-(const Vector &a);
    Vector operator-(const double a);
    Vector operator*(const Vector &a);
    Vector operator*(const double a);
    Vector operator/(const Vector &a);
    Vector operator/(const double a);
    const double &operator[](int i) const;
    double &operator[](int i);

    double dot(const Vector &a);
    Vector cross_product(const Vector &a);
    Vector pow(const double a);
    Vector max(const double a);
    Vector min(const double a);
    int argmin();
    int get_longest();
    double norm();
    Vector normalize();

private:
    double coords[3];
};

struct Intersection{
    bool is_intersect;
    double t;
    Vector P;
    Vector N;
    int id;
};

class Ray{
public:
    explicit Ray(Vector O, Vector u, double refraction_index);
    Vector O;
    Vector u;
    double refraction_index;
};

class Geometry{
public:
    virtual Intersection intersect(Ray &r) = 0;
    Vector albedo;
    bool mirror;
    bool transparent;
    double refraction_index;
};

class Sphere : public Geometry{
public:
    explicit Sphere(Vector C, Vector albedo, double R, bool mirror, bool transparent, double refraction_index = 1.5);
    virtual Intersection intersect(Ray &r);
    Vector C;
    double R;
};


class Bbox {
public:
    Bbox();
    Bbox(Vector &b_min,Vector &b_max);
    Vector b_min;
    Vector b_max;
    bool intersect(Ray &r,double &inter_distance);
    Vector compute_diag();
};

class Node{
public:
    int starting_triangle;
    int ending_triangle;
    Bbox bbox;
    Node* child_left;
    Node* child_right;
    bool has_children;
    Node(){}
};

class TriangleIndices{
//NOTE: I got inspired by the implementation in the link on Moodle
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group){};
    int vtxi, vtxj, vtxk; // for vertex coordinates 
    int uvi, uvj, uvk;    // for uv coordinates 
    int ni, nj, nk;       // for normals 
    int group;            // for face group
};

class TriangleMesh : public Geometry{
public:
    virtual Intersection intersect(Ray &r);
    ~TriangleMesh() {}
    TriangleMesh(Vector albedo, double scaling_factor, Vector translation);
    void readOBJ(const char *obj);
    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
    double scaling_factor;
    Vector translation;
    Bbox compute_bbox(int starting_triangle, int ending_triangle);
    Node* root;
    void performingBVH(Node* node,int starting_triangle, int ending_triangle);
    Vector find_barycenter(int triangle_indice);
};