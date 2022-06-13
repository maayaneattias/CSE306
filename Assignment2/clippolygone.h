#include <iostream>
#include <vector>

class Vector
{
public:
    explicit Vector(double x = 0., double y = 0.);
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
    double norm();
    Vector normalize();

private:
    double coords[3];
};

class Edge{
    public:
        explicit Edge(const Vector &a,const Vector &b);
        Vector point_a;
        Vector point_b;
};

class Polygon {
    public:
        explicit Polygon();
        std::vector<Edge> edges;
        std::vector<Vector> vertices;
};