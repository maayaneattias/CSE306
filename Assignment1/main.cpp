#include <cmath>
#include <chrono>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"
#include "scene.cpp"


int main()
{
    using std::chrono::duration;
    using std::chrono::duration_cast;
    using std::chrono::high_resolution_clock;
    using std::chrono::milliseconds;
    
    // You can modifiy those parameters
    auto t1 = high_resolution_clock::now();
    int W = 512;
    int H = 512;
    int n_rays = 32;
    int scene_type = 2;
    bool bvh = true;
    bool antialiasing = false;
    bool fresnel = false;
    bool indirect_lightning = false;
    if (fresnel) n_rays = 32;
    double I = 2e10;
    int ray_depth = 5;
    double fov = 1.047; // 60 deg
    Vector camera_pos = Vector(0, 0, 55);
    Vector S = Vector(-10, 20, 40);
    //end of parameters setting

    std::vector<unsigned char> image(W * H * 3, 0);
    Scene scene = Scene(I, S, 1,fresnel,bvh,scene_type);
    #pragma omp parallel for schedule(dynamic, 1) //parrallelization

    for (int i = 0; i < H; i++)
    {
        for (int j = 0; j < W; j++)
        {
            Vector colors = Vector(0,0,0);
            for (int k = 0; k < n_rays;k++)
            {
                Ray ray = center_ray(camera_pos,W,H,fov, scene.refraction_index,i,j,antialiasing);
                colors += scene.get_color(ray, ray_depth,indirect_lightning);
            }
            colors = colors/n_rays;
            colors = colors.pow(1. / 2.2).min(255).max(0);
            image[(i * W + j) * 3] = colors[0];
            image[(i * W + j) * 3 + 1] = colors[1];
            image[(i * W + j) * 3 + 2] = colors[2];
        };
    };
    stbi_write_png("imagelast.png", W, H, 3, &image[0], 0); //save image
    auto t2 = high_resolution_clock::now();
    //auto ms_int = duration_cast<milliseconds>(t2 - t1); //time taken
    duration<double, std::milli> ms_double = t2 - t1;

    //std::cout << ms_int.count() << "ms\n";
    std::cout << ms_double.count() << "ms"; //time taken
    return 0;
};
