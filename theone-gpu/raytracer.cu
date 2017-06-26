// [header]
// A very basic raytracer example.
// [/header]
// [compile]
// c++ -o raytracer -O3 -Wall raytracer.cpp
// [/compile]
// [ignore]
// Copyright (C) 2012  www.scratchapixel.com
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// [/ignore]
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>
#include <cstring>
#include <sys/time.h>

#if defined __linux__ || defined __APPLE__
// "Compiled for Linux
#else
// Windows doesn't define these values by default, Linux does
#define M_PI 3.141592653589793
#define INFINITY 1e8
#endif

#ifndef TILE_WIDTH
#define TILE_WIDTH 8
#endif

double rtclock()
{
    struct timezone Tzp;
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, &Tzp);
    if (stat != 0) printf("Error return from gettimeofday: %d",stat);
    return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}

class UnifiedMemoryClass
{
public:
    void* operator new(size_t len) {
        void *ptr;
        cudaMallocManaged(&ptr, len);
        return ptr;
    }

    void operator delete(void *ptr) {
        cudaFree(ptr);
    }

    void* operator new[](size_t len) {
        void *ptr;
        cudaMallocManaged(&ptr, len);
        return ptr;
    }
};

class Vec3 : public UnifiedMemoryClass
{
public:
    float x, y, z;
    __device__ __host__ Vec3() : x(float(0)), y(float(0)), z(float(0)) {}
    __device__ __host__ Vec3(float xx) : x(xx), y(xx), z(xx) {}
    __device__ __host__ Vec3(float xx, float yy, float zz) : x(xx), y(yy), z(zz) {}
    __device__ Vec3& normalize()
    {
        float nor2 = length2();
        if (nor2 > 0) {
            float invNor = 1 / sqrt(nor2);
            x *= invNor, y *= invNor, z *= invNor;
        }
        return *this;
    }
    __device__ Vec3 operator * (const float &f) const { return Vec3(x * f, y * f, z * f); }
    __device__ static Vec3 mult (const float &f, const Vec3 v) {
        return Vec3(v.x * f, v.y * f, v.z * f);
    }
    __device__ Vec3 operator * (const Vec3 &v) const { return Vec3(x * v.x, y * v.y, z * v.z); }
    __device__ static Vec3 mult (const Vec3 &v1, const Vec3 &v2) {
        return Vec3(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z);
    }
    __device__ float dot(const Vec3 &v) const { return x * v.x + y * v.y + z * v.z; }
    __device__ static float dot(const Vec3 &v1, const Vec3 &v2) {
        return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    }
    __device__ Vec3 operator - (const Vec3 &v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
    __device__ static Vec3 sub (const Vec3 &v1, const Vec3 &v2) {
        return Vec3(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
    }
    __device__ Vec3 operator + (const Vec3 &v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
    __device__ static Vec3 add (const Vec3 &v1, const Vec3 &v2) {
        return Vec3(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
    }
    __device__ Vec3& operator += (const Vec3 &v) { x += v.x, y += v.y, z += v.z; return *this; }
    __device__ Vec3& operator *= (const Vec3 &v) { x *= v.x, y *= v.y, z *= v.z; return *this; }
    __device__ Vec3 operator - () const { return Vec3(-x, -y, -z); }
    __device__ static Vec3 neg (const Vec3 &v) {
        return Vec3(-v.x, -v.y, -v.z);
    }
    __device__ float length2() const { return x * x + y * y + z * z; }
    __device__ static float length2(const Vec3 &v) {
        return v.x * v.x + v.y * v.y + v.z * v.z;
    }
    __device__ float length() const { return sqrt(length2()); }
    __device__ static float length(const Vec3 &v) {
        return sqrt(length2(v));
    }
    __host__ friend std::ostream & operator << (std::ostream &os, const Vec3 &v)
    {
        os << "[" << v.x << " " << v.y << " " << v.z << "]";
        return os;
    }
};

__device__ Vec3 sub (const Vec3 &v1, const Vec3 &v2) {
    return Vec3(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

typedef Vec3 Vec3f;
typedef Vec3 RGB;

class Sphere : public UnifiedMemoryClass
{
public:
    Vec3f center;                           /// position of the sphere
    float radius, radius2;                  /// sphere radius and radius^2
    Vec3f surfaceColor, emissionColor;      /// surface color and emission (light)
    float reflection;         /// surface transparency and reflectivity
    __device__ __host__ Sphere(
        const Vec3f &c,
        const float &r,
        const Vec3f &sc,
        const float &refl = 0,
        const Vec3f &ec = 0) :
        center(c), radius(r), radius2(r * r), surfaceColor(sc), emissionColor(ec),
        reflection(refl)
    { /* empty */ }
    __device__ __host__ Sphere(){}
    //[comment]
    // Compute a ray-sphere intersection using the geometric solution
    //[/comment]
    __device__ bool intersect(const Vec3f &rayorig, const Vec3f &raydir, float &t0, float &t1) const
    {
        Vec3f l = center - rayorig;
        float tca = l.dot(raydir);
        if (tca < 0) return false;
        float d2 = l.dot(l) - tca * tca;
        if (d2 > radius2) return false;
        float thc = sqrtf(radius2 - d2);
        t0 = tca - thc;
        t1 = tca + thc;

        return true;
    }
};

//[comment]
// This variable controls the maximum recursion depth
//[/comment]
#ifndef MAX_RAY_DEPTH
#define MAX_RAY_DEPTH 10
#endif

__device__ float mix(const float &a, const float &b, const float &mix)
{
    return b * mix + a * (1 - mix);
}

//[comment]
// This is the main trace function. It takes a ray as argument (defined by its origin
// and direction). We test if this ray intersects any of the geometry in the scene.
// If the ray intersects an object, we compute the intersection point, the normal
// at the intersection point, and shade this point using this information.
// Shading depends on the surface property (is it transparent, reflective, diffuse).
// The function returns a color for the ray. If the ray intersects an object that
// is the color of the object at the intersection point, otherwise it returns
// the background color.
//[/comment]
__global__ void trace(
    const Sphere* spheres,
    const unsigned n,
    Vec3f *image,
    unsigned width,
    unsigned height)
{
    unsigned x = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned y = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned tid = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned blockSize = blockDim.x * blockDim.y;

    extern __shared__ Sphere s_spheres[];

    if (tid < n) {
        s_spheres[tid] = spheres[tid];
    }

    if (blockSize < n && tid == 0) {
        for (unsigned i = blockSize; i < n; i++) {
            s_spheres[i] = spheres[i];
        }
    }

    __syncthreads();

    if (x < width && y < height) {
        float fov = 30;
        float aspectratio = width / float(height);
        float angle = tan(M_PI * 0.5 * fov / 180.);
        float xx = (2 * ((x + 0.5) / float(width)) - 1) * angle * aspectratio;
        float yy = (1 - 2 * ((y + 0.5) / float(height))) * angle;
        Vec3f rayorig(0);
        Vec3f raydir(xx, yy, -1);
        raydir.normalize();

        int depth = 0;
        Vec3f reflectionMultStack[MAX_RAY_DEPTH+1];
        Vec3f reflectionAddStack[MAX_RAY_DEPTH+1];
        Vec3f result;

        while (true) {
            //if (raydir.length() != 1) std::cerr << "Error " << raydir << std::endl;
            float tnear = INFINITY;
            const Sphere* sphere = NULL;
            // find intersection of this ray with the sphere in the scene
            for (unsigned i = 0; i < n; ++i) {
                float t0 = INFINITY, t1 = INFINITY;
                if (s_spheres[i].intersect(rayorig, raydir, t0, t1)) {
                    if (t0 < 0) t0 = t1;
                    if (t0 < tnear) {
                        tnear = t0;
                        sphere = &s_spheres[i];
                    }
                }
            }
            // if there's no intersection return black or background color
            if (!sphere) {
                //return Vec3f(2);
                result = Vec3f(2);
                break;
            }
            Vec3f surfaceColor = 0; // color of the ray/surfaceof the object intersected by the ray
            Vec3f phit = rayorig + raydir * tnear; // point of intersection
            Vec3f nhit = phit - sphere->center; // normal at the intersection point
            nhit.normalize(); // normalize normal direction
            // If the normal and the view direction are not opposite to each other
            // reverse the normal direction. That also means we are inside the sphere so set
            // the inside bool to true. Finally reverse the sign of IdotN which we want
            // positive.
            float bias = 1e-4; // add some bias to the point from which we will be tracing
            if (raydir.dot(nhit) > 0) nhit = -nhit;
            if (sphere->reflection > 0 && depth < MAX_RAY_DEPTH) {
                float facingratio = -raydir.dot(nhit);
                // change the mix value to tweak the effect
                float fresneleffect = mix(pow(1 - facingratio, 3), 1, 0.1);
                // compute reflection direction (not need to normalize because all vectors
                // are already normalized)
                Vec3f refldir = raydir - nhit * 2 * raydir.dot(nhit);
                refldir.normalize();
                reflectionMultStack[depth] = sphere->surfaceColor * fresneleffect;
                reflectionAddStack[depth] = sphere->emissionColor;
                rayorig = phit + nhit * bias;
                raydir = refldir;
                depth++;
                continue;
            }
            else {
                // it's a diffuse object, no need to raytrace any further
                for (unsigned i = 0; i < n; ++i) {
                    if (s_spheres[i].emissionColor.x > 0) {
                        // this is a light
                        Vec3f transmission = 1;
                        Vec3f lightDirection = s_spheres[i].center - phit;
                        lightDirection.normalize();
                        for (unsigned j = 0; j < n; ++j) {
                            if (i != j) {
                                float t0, t1;
                                if (s_spheres[j].intersect(phit + nhit * bias, lightDirection, t0, t1)) {
                                    transmission = 0;
                                    break;
                                }
                            }
                        }
                        surfaceColor += sphere->surfaceColor * transmission *
                        fmaxf(float(0), nhit.dot(lightDirection)) * s_spheres[i].emissionColor;
                        // return surfaceColor + sphere->emissionColor;
                    }
                }
                result = surfaceColor + sphere->emissionColor;
                break;
            }
        }

        for (depth = depth - 1; depth >= 0; depth--) {
            result = result * reflectionMultStack[depth] + reflectionAddStack[depth];
        }

        image[y*width+x] = result;
    }
}

//[comment]
// Main rendering function. We compute a camera ray for each pixel of the image
// trace it and return a color. If the ray hits a sphere, we return the color of the
// sphere at the intersection point, else we return the background color.
//[/comment]
void render(Vec3f* image, unsigned width, unsigned height, const Sphere* spheres, const unsigned n)
{
    dim3 dimGrid(ceil((float) width / TILE_WIDTH), ceil((float)height / TILE_WIDTH));
	dim3 dimBlock(TILE_WIDTH, TILE_WIDTH);

    trace<<<dimGrid, dimBlock, n * sizeof(Sphere)>>>(spheres, n, image, width, height);
    cudaDeviceSynchronize();
}

void save(const char* filename, Vec3f* image, unsigned width, unsigned height) {
    // Save result to a PPM image (keep these flags if you compile under Windows)
    std::ofstream ofs(filename, std::ios::out | std::ios::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (unsigned i = 0; i < width * height; ++i) {
        ofs << (unsigned char)(std::min(float(1), image[i].x) * 255) <<
               (unsigned char)(std::min(float(1), image[i].y) * 255) <<
               (unsigned char)(std::min(float(1), image[i].z) * 255);
    }
    ofs.close();
}

//[comment]
// In the main function, we will create the scene which is composed of some spheres
// and some light (which is also a sphere). Then, once the scene description is complete
// we render that scene, by calling the render() function.
//[/comment]
int main(int argc, char **argv)
{
    if (argc < 3) {
        std::cout << "No file detected" << '\n';
        return 1;
    }
    FILE* scene = NULL;
    scene = fopen(argv[1], "r");

    if (scene == NULL) {
        std::cout << "Error when reading file" << '\n';
    }
    unsigned width, height;
    unsigned s, l;
    fscanf(scene, "%u %u\n", &width, &height);
    fscanf(scene, "%u %u\n", &s, &l);
    Sphere* spheres = new Sphere[s+l];
    // spheres
    for(int i=0; i<s; i++) {
      float x, y, z, r, refl;
      fscanf(scene, "%f %f %f %f ", &x, &y, &z, &r);
      Vec3f position(x, y, z);
      fscanf(scene, "%f %f %f %f\n", &x, &y, &z, &refl);
      RGB color(x, y, z);
      // position, radius, surface color, reflectivity, transparency, emission color
      spheres[i] = Sphere(position, r, color, refl);
    }
    // lights
    for(int i=s; i<s+l; i++) {
      float x, y, z, r, refl, ec;
      fscanf(scene, "%f %f %f %f ", &x, &y, &z, &r);
      Vec3f position(x, y, z);
      fscanf(scene, "%f %f %f %f %f\n", &x, &y, &z, &refl, &ec);
      RGB color(x, y, z);
      // position, radius, surface color, reflectivity, transparency, emission color
      spheres[i] = Sphere(position, r, color, refl, ec);
    }
    Vec3f *image = new Vec3f[width * height];

    double runTime = rtclock();
    render(image, width, height, spheres, s+l);
    runTime = rtclock() - runTime;
    std::cout << "Run time: " << runTime << '\n';

    save(argv[2], image, width, height);

	delete image;
    delete spheres;

    return 0;
}
