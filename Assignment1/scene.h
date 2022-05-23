#include <iostream>
#include <vector>
#include "auxiliary.cpp"

class Scene{
public:
    Scene(double I, Vector &S, double refraction_index, bool fresnel,bool bvh, int scene_type);
    Intersection intersect(Ray &r);
    Vector get_color(Ray &ray, int ray_depth,bool indirect_lightning);
    double refraction_index;

private:
    std::vector<Geometry *> geometries; //changed in Lab 3 pointers to Geometry now
    double I; 
    Vector S; 
    bool fresnel;
};