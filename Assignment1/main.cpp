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
    bool bvh = false;
    bool antialiasing = false;
    bool fresnel = false;
    bool indirect_lightning = false;
    if (fresnel == true) {
        n_rays = 32;
    }
    double I = 3e10;
    int ray_depth = 5;
    double field_of_view = 1.047; // 60 deg
    Vector position_camera = Vector(0, 0, 55);
    Vector S = Vector(-10, 20, 40);
    //end of parameters setting

    std::vector<unsigned char> result_grid(W * H * 3, 0);
    Scene scene = Scene(I, S, 1,fresnel,bvh,scene_type);
    #pragma omp parallel for schedule(dynamic, 1) //parrallelization

    for (int i = 0; i < H; i++)
    {
        for (int j = 0; j < W; j++)
        {
            Vector colors = Vector(0,0,0);
            for (int k = 0; k < n_rays;k++)
            {
                Ray ray = ray_centering(position_camera,W,H,field_of_view, scene.refraction_index,i,j,antialiasing);
                colors += scene.get_color(ray, ray_depth,indirect_lightning);
            }
            colors = colors/n_rays;
            colors = colors.pow(1. / 2.2).min(255).max(0);
            result_grid[(i * W + j) * 3] = colors[0];
            result_grid[(i * W + j) * 3 + 1] = colors[1];
            result_grid[(i * W + j) * 3 + 2] = colors[2];
        };
    };
    stbi_write_png("result.png", W, H, 3, &result_grid[0], 0); //save result_grid
    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;
    std::cout << ms_double.count() << "ms"; //time taken
    return 0;
};
