#include <iostream>
#include <random>
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../stb/stb_image_write.h"
#include "../stb/stb_image.h"

/*
NOTE: I implemented the second choice: image retargeting algorithm for shrinking images (lecture notes sec. 3.3) using 

Seam Carving method steps:
-Determine associated energy
-compute a cumulated energy map
-find the starting point of the seam
-progressively build the seam
-remove the seam
*/

std::vector<unsigned char> forenergy(std::vector<unsigned char> vect_intensity, int h, int w)
{
    int energy;
    int temp;
    int size = w*h;
    std::vector<unsigned char> vect_energy(size, 0);
    for (int i = 0; i < h; i++)
    {
        for (int j = 0; j < w; j++)
        {
            temp = 0;
            energy = 0;
            if (j + 1 < w)
            {
                temp += vect_intensity[i * w + j + 1];
            }
            if (j - 1 >= 0)
            {
                temp -= vect_intensity[i * w + j - 1];
            }
            energy += std::abs(temp);

            temp = 0;
            if (i + 1 < h)
            {
                temp += vect_intensity[(i + 1) * w + j];
            }
            if (i - 1 >= 0)
            {
                temp -= vect_intensity[(i - 1) * w + j];
            }
            energy += std::abs(temp);
            vect_energy[i * w + j] = energy;
        }
    }
    return vect_energy;
}

std::vector<unsigned char> forcumulate(std::vector<unsigned char> vect_energy, int h, int w)
{
    int energy;
    int size = w*h;
    std::vector<unsigned char> vect_cumulate(size, 0);
    for (int i = 0; i < h; i++)
    {
        for (int j = 0; j < w; j++)
        {
            vect_cumulate[i * w + j] = vect_energy[i * w + j];
            if (j + 1 < w && i - 1 >= 0 && j - 1 >= 0)
            {
                vect_cumulate[i * w + j] += std::min(std::min(vect_energy[(i - 1) * w + j + 1], vect_energy[(i - 1) * w + j]), vect_energy[(i - 1) * w + j - 1]);;
            }
        }
    }
    return vect_cumulate;
}

std::vector<unsigned char> forintensity(std::vector<unsigned char> imagee, int h, int w)
{
    int size = w*h;
    std::vector<unsigned char> vect_intensity(size,0);
    for (int i = 0; i < h; i++)
    {
        for (int j = 0; j < w; j++)
        {
            vect_intensity[i * w + j ] = imagee[(i * w + j) * 3]+imagee[(i * w + j) * 3 + 1]+imagee[(i * w + j) * 3+2];
        }
    }
    return vect_intensity;
}

void translation(std::vector<unsigned char> imagee, int w, int c, int r)
{

    for (int i = c; i < w - 1; i++)
    {
        imagee[(r * (w-1) + i) * 3] = imagee[(r * w + i) * 3 + 3];
        imagee[(r * (w-1) + i) * 3 + 1] = imagee[(r * w + i) * 3 + 4];
        imagee[(r * (w-1) + i) * 3 + 2] = imagee[(r * w + i) * 3 + 5];
    }
}

std::vector<unsigned char> cropping(std::vector<unsigned char> image, std::vector<int> secondaryvector, int h,int w){
    std::vector<unsigned char> cropped((w-1)*h*3,0);
    for (int i=0; i<h; i++){
        int value_path = secondaryvector[i];
        for (int j=0; j<w-1; j++){
            if (j>=value_path){
                cropped[(i*(w-1)+j)*3 + 0] = image[(i*(w)+j+1)*3 + 0]; 
                cropped[(i*(w-1)+j)*3 + 1] = image[(i*(w)+j+1)*3 + 1];
                cropped[(i*(w-1)+j)*3 + 2] = image[(i*(w)+j+1)*3 + 2];  
            }
            else{
                cropped[(i*(w-1)+j)*3 + 0] = image[(i*(w)+j)*3 + 0]; 
                cropped[(i*(w-1)+j)*3 + 1] = image[(i*(w)+j)*3 + 1];
                cropped[(i*(w-1)+j)*3 + 2] = image[(i*(w)+j)*3 + 2];  
            }
        }
    }
    return cropped;
}


int main()
{
    int w, h, forload;
    std::vector<unsigned char> vect_intensity,vect_energy,vect_cumulate;
    int N = 40;
    unsigned char *image = stbi_load("test.jpg", &w, &h, &forload, 0);
    std::vector<unsigned char> imagee(w*h*3,0);
    for (int i=0; i<h; i++){
        for (int j = 0; j <w; j++){
            imagee[(i*w + j) * 3] = image[(i*w + j) * 3];
            imagee[(i*w + j) * 3 + 1] = image[(i*w + j) * 3 + 1];
            imagee[(i*w + j) * 3 + 2] = image[(i*w + j) * 3 + 2];
        }
    }
    stbi_image_free(image);
    std::vector<int> secondaryvector(h);
    for (int k = 0; k < N; k++)
    {
        vect_intensity = forintensity(imagee, h, w);
        vect_energy = forenergy(vect_intensity, h, w);
        vect_cumulate = forcumulate(vect_energy, h, w);
        double min_set = std::numeric_limits<double>::max();
        int c;
        for (int i = 0; i < w; i++)
        {
            if (vect_cumulate[(h - 1) * w + i] < min_set)
            {
                min_set = vect_cumulate[(h - 1) * w + i];
                c = i;
            }
        }
        secondaryvector[h-1] = c;
        int temp,new_w;
        
        for (int i = h - 1; i > 0; i--)
        {
            min_set = std::numeric_limits<double>::max();
            if (c - 1 >= 0 && vect_cumulate[(i-1) * w + c - 1] < min_set)
            {
                temp = c - 1;
                min_set = vect_cumulate[(i-1) * w + c - 1];
            }
            if (vect_cumulate[(i-1) * w + c] < min_set)
            {
                temp = c;
                min_set = vect_cumulate[(i-1) * w + c];
            }
            if (c + 1 < w && vect_cumulate[(i-1) * w + c + 1] < min_set)
            {
                temp = c + 1;
                min_set = vect_cumulate[(i-1) * w + c + 1];
            }
            c = temp;
            secondaryvector[i-1] = c;
            translation(imagee, w, c, i - 1);
        }
        imagee = cropping(imagee, secondaryvector,h, w);
        w -= 1;
    }
    stbi_write_png("result.png", w, h, 3, &imagee[0], 0);
    
}


