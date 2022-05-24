#include "scene.h" 
#include <iostream> 
#include <cmath>    
#include <random>
#include <algorithm>
#include <list>

Scene::Scene(double I, Vector &S, double refraction_index, bool fresnel, bool bvh, int scene_type)
{
    this->I = I;
    this->S = S;
    this->refraction_index = refraction_index;
    this->fresnel = fresnel;

    Sphere *top_sphere = new Sphere(Vector(0, 1000, 0), Vector(1, 0, 0), 940, false, false);
    this->geometries.push_back(top_sphere);

    Sphere *right_sphere = new Sphere(Vector(0, 0, -1000), Vector(0, 1, 0), 940, false, false);
    this->geometries.push_back(right_sphere);

    Sphere *bottom_sphere = new Sphere(Vector(0, -1000, 0), Vector(0, 0, 1), 990, false, false);
    this->geometries.push_back(bottom_sphere);

    Sphere *left_sphere = new Sphere(Vector(0, 0, 1000), Vector(1, 0, 1), 940, false, false);
    this->geometries.push_back(left_sphere);

    Sphere *aqua_sphere = new Sphere(Vector(-1000, 0, 0), Vector(0, 1, 1), 940, false, false);
    this->geometries.push_back(aqua_sphere);

    Sphere *yellow_sphere = new Sphere(Vector(1000, 0, 0), Vector(1, 1, 0), 940, false, false);
    this->geometries.push_back(yellow_sphere);

    if (scene_type == 1)
    {
        Sphere *center_left_sphere = new Sphere(Vector(-20, 0, 0), Vector(0.5, 0.5, 0.5), 10, true, false);
        this->geometries.push_back(center_left_sphere);

        Sphere *center_sphere = new Sphere(Vector(0, 0, 0), Vector(0.5, 0.5, 0.5), 10, false, true);
        this->geometries.push_back(center_sphere);

        Sphere *center_right_sphere = new Sphere(Vector(20, 0, 0), Vector(0.5, 0.5, 0.5), 10, false, true);
        this->geometries.push_back(center_right_sphere);

        Sphere *center_right_sphere_hollow = new Sphere(Vector(20, 0, 0), Vector(0.5, 0.5, 0.5), 10 - 1, false, true, 1);
        this->geometries.push_back(center_right_sphere_hollow);
    }
    else if (scene_type == 2)
    {
        double scaling_factor = 0.6; //can be changed
        //Vector albedo = Vector(0.3, 0.2, 0.25);  
        Vector albedo =  Vector(1, 1, 1) ;
        Vector translation = Vector(0,-10,0);
        TriangleMesh *meshes = new TriangleMesh(albedo, scaling_factor, translation);
        meshes->readOBJ("cat/cat.obj");
        for (int i = 0; i < meshes->vertices.size(); i++)
        {
            meshes->vertices[i] = meshes->vertices[i] * scaling_factor;
            meshes->vertices[i] = meshes->vertices[i] + translation;
        }
        if (bvh)
            meshes->order_BVH(meshes->root, 0, meshes->indices.size());
        else //bounding box
        {
            meshes->root->bbox = meshes->compute_bbox(0, meshes->indices.size());
            meshes->root->has_children = true;
            meshes->root->starting_triangle = 0;
            meshes->root->ending_triangle = meshes->indices.size();
        }
        this->geometries.push_back(meshes);
    }
};

Intersection Scene::intersect(Ray &r)
{
    Intersection result;
    result.is_intersect = false;
    double last_t = std::numeric_limits<double>::max();

    for (std::vector<Geometry *>::size_type i = 0; i != geometries.size(); i++)
    {
        Intersection intersection = geometries[i]->intersect(r);
        if (intersection.is_intersect && intersection.t < last_t)
        {
            last_t = intersection.t;
            intersection.id = i;
            result = intersection;
        };
    };
    return result;
}

Vector Scene::get_color(Ray &ray, int ray_depth, bool indirect_lightning)
{
    if (ray_depth < 0)
    {
        return Vector(0, 0, 0);
    }
    Intersection intersection = this->intersect(ray);
    if (not intersection.is_intersect)
    {
        return Vector(0, 0, 0);
    }
    if (geometries[intersection.id]->mirror)
    {
        Ray reflected_ray = Ray(intersection.P + intersection.N * 1e-3, ray.u - intersection.N * (2 * intersection.N.dot(ray.u)), refraction_index);
        return get_color(reflected_ray, ray_depth - 1, indirect_lightning);
    }

    else if (geometries[intersection.id]->transparent)
    {
        bool is_exiting = false;
        Vector original_N = intersection.N;
        if (ray.u.dot(intersection.N) > 0)
        {
            intersection.N = intersection.N * (-1);
            is_exiting = true;
        }
        double n2;
        if (is_exiting)
        {
            n2 = refraction_index;
        }
        else
        {
            n2 = geometries[intersection.id]->refraction_index;
        }

        //fresnel
        double k0 = (ray.refraction_index - n2) * (ray.refraction_index - n2) / ((ray.refraction_index + n2) * (ray.refraction_index + n2));
        double N_dot_u = abs(intersection.N.dot(ray.u));
        double R = k0 + (1 - k0) * (1 - N_dot_u) * (1 - N_dot_u) * (1 - N_dot_u) * (1 - N_dot_u) * (1 - N_dot_u);
        double u = uniform(engine);
        if (fresnel && u < R)
        { //reflection
            Ray reflected_ray = Ray(intersection.P + intersection.N * 1e-3, ray.u - intersection.N * (2 * intersection.N.dot(ray.u)), refraction_index);
            return get_color(reflected_ray, ray_depth - 1, indirect_lightning);
        }

        else
        {
            double dot_u_n = ray.u.dot(intersection.N);
            double refraction = ray.refraction_index / n2;
            double x = 1 - refraction * refraction * (1 - dot_u_n * dot_u_n);
            if (x > 0)
            { //refraction
                Vector w_t = (ray.u - intersection.N * dot_u_n) * refraction;
                Vector w_n = intersection.N * (-1) * sqrt(x);
                Vector w = w_t + w_n;
                Ray transparent_ray = Ray(intersection.P - intersection.N * 1e-3, w, n2);
                return get_color(transparent_ray, ray_depth - 1, indirect_lightning);
            }
            else
            { //reflection
                Ray total_internal_refraction = Ray(intersection.P + intersection.N * 1e-3, ray.u - intersection.N * (2 * intersection.N.dot(ray.u)), n2);
                return get_color(total_internal_refraction, ray_depth - 1, indirect_lightning);
            }
        }
    }
    else
    {
        Vector p_s = this->S - intersection.P;
        double d = p_s.norm();
        Vector w = p_s / d;
        Ray r = Ray(intersection.P + intersection.N * 1e-3, w, refraction_index);
        Intersection visibility = this->intersect(r);
        double hidden = not(visibility.is_intersect && visibility.t < d);
        Vector L0 = geometries[intersection.id]->albedo * intersection.N.dot(w) * this->I * hidden / (4 * M_PI * d * M_PI * d);
        if (indirect_lightning)
        {
            Ray random_ray = Ray(intersection.P + intersection.N * 1e-3, random_cos(intersection.N), refraction_index);
            L0 += get_color(random_ray, ray_depth - 1, indirect_lightning) * geometries[intersection.id]->albedo;
        }
        return L0;
    }
}
