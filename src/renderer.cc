//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Copyright (C) 2009  Sebastian Mach (*1983)
// * phresnel/at/gmail/dot/com
// * http://phresnel.org
// * http://picogen.org
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include <iostream>
#include <SDL/SDL.h>

#include "renderer.hh"
#include "john-ratcliff-perlin.def.hh"


Renderer::Renderer (unsigned int width_, unsigned int height_) 
: width (width_), height (height_)
, screen (0)
, heightfield_scale (10.0f)
, sunDirection (normalise (Vector3d (0.0,1.0,6.0)))
{
        // Initialize SDL video.
        if (SDL_Init (SDL_INIT_VIDEO) < 0) {
                std::cerr <<  "unable to init SDL:\n" 
                          << SDL_GetError() 
                          << std::endl;
        }
        atexit(SDL_Quit);

        // Create a new window.
        screen = SDL_SetVideoMode (
                        width, height, 
                        32,
                        SDL_HWSURFACE | SDL_DOUBLEBUF
                );

        if (0 == screen) {
                std::cerr << "Unable to set video to " 
                     << width << "x" << height
                     << ":\n" << SDL_GetError() 
                     << std::endl;
        }
        
        initHeightfield();
        initSkymap();
}



Renderer::~Renderer() {
}



void Renderer::initHeightfield () {
        
        heightfield = new Vertex1 [heightfield_size*heightfield_size];

        std::cout << "initializing heightfield" << std::endl;
        // syn: Perlin (float freq, float amp, int seed)
        struct Fun {
                Perlin<9> perlin;
                Fun () : perlin (0.025, 4.0, 32) {
                }
                
                inline float operator () (float fx, float fz) const {
                        //return 2.0 * sin (fx) * sin (fz);
                        return -3.0 + 4.0 * (1.0f - fabsf (perlin.Get (fx, fz)));
                }
        };
        Fun fun;
        
        for (int iz=0; iz<heightfield_size; ++iz) {
                for (int ix=0; ix<heightfield_size; ++ix) {
                        float fx = (float)ix / heightfield_scale;
                        float fz = (float)iz / heightfield_scale;
                        heightfield [makeHeightfieldIndex (ix, iz)].height = fun (fx, fz);                        
                }
        }
        
        std::cout << "kludging together normal map" << std::endl;
        for (int iz=0; iz<heightfield_size-1; ++iz) {
                for (int ix=0; ix<heightfield_size-1; ++ix) {
                        const double
                                h00 = heightfield [makeHeightfieldIndex (ix, iz)].height, 
                                h01 = heightfield [makeHeightfieldIndex (ix+1, iz)].height, 
                                h10 = heightfield [makeHeightfieldIndex (ix, iz+1)].height, 
                                //h11 = field [(iz+1)*size + ix+1].height,
                                hu = h00 - h01,
                                hv = h00 - h10,
                                du = 1.0 / heightfield_scale,
                                dv = 1.0 / heightfield_scale,
                                u_ [3] = {du, hu, 0.0f},
                                v_ [3] = {0.0f, hv, dv},
                                ilen_u_ = 1.0f / sqrtf (u_[0]*u_[0] + u_[1]*u_[1] + u_[2]*u_[2]),
                                ilen_v_ = 1.0f / sqrtf (v_[0]*v_[0] + v_[1]*v_[1] + v_[2]*v_[2]),
                                u [3] = {ilen_u_*u_[0], ilen_u_*u_[1], ilen_u_*u_[2] },
                                v [3] = {ilen_v_*v_[0], ilen_v_*v_[1], ilen_v_*v_[2] },
                                cross [3] = {
                                        u[1]*v[2] - u[2]*v[1], 
                                        u[2]*v[0] - u[0]*v[2], 
                                        u[0]*v[1] - u[1]*v[0], 
                                };
                        ;
                        const unsigned int index = makeHeightfieldIndex (ix, iz);
                        heightfield [index].normal[0] = cross [0];
                        heightfield [index].normal[1] = fabsf (cross [1]);
                        heightfield [index].normal[2] = cross [2];
                        
                        /*std::cout << "hxx="
                                  << h00 << ',' 
                                  << h01 << ',' 
                                  << h10 << ','                                           
                                  << h11 << '\n'
                                  << hu << ',' 
                                  << hv << '\n' ;
                        std::cout << "u="
                                  << u[0] << ',' 
                                  << u[1] << ',' 
                                  << u[2] << '\n';
                        std::cout << "v="
                                  << v[0] << ',' 
                                  << v[1] << ',' 
                                  << v[2] << '\n';
                        std::cout << "{" << field [iz*size + ix].normal[0]
                                    << "," << field [iz*size + ix].normal[1]
                                    << "," << field [iz*size + ix].normal[2] << "}" << std::endl;*/
                }
        }
        std::cout << "done." << std::endl;                
        //exit(0);
}




void Renderer::initSkymap () {
        skymap_fsize = static_cast<float>(skymap_size-1);
        
        skymap = new float [skymap_size*skymap_size][4];
        
        picogen::Preetham preetham;
        switch (1) {
        case 0:
                preetham.setColorFilter (picogen::Color (0.3, 0.2, 0.15));
                preetham.setSunDirection (picogen::Vector3d (0.4, 0.2, 0.4).normal());
                break;
        case 1:        
                preetham.setColorFilter (picogen::Color (0.3, 0.3, 0.3));
                preetham.setSunDirection (picogen::Vector3d (sunDirection.x[0], sunDirection.y[0], sunDirection.z[0]));
                break;
        };

        preetham.setSunColor (picogen::Color (500));                
        preetham.setTurbidity (2.5);
        preetham.invalidate();
        
        std::cout << "initializing skymap" << std::endl;
        for (size_t iz=0; iz<skymap_size; ++iz) {
                for (size_t ix=0; ix<skymap_size; ++ix) {
                        const float
                                x = -1.0 + 2.0 * (static_cast<float>(ix) / skymap_fsize),
                                z = -1.0 + 2.0 * (static_cast<float>(iz) / skymap_fsize),
                                h = f2_to_hemisphere (x,z)
                        ;
                        
                        picogen::Ray peek (
                                picogen::Vector3d (0.0, 0.0, 0.0),
                                picogen::Vector3d (x, h, z)
                        );
                        picogen::Color col;
                        preetham.shade (col, peek);
                        /*if (ix%10==0 || iz%10==0)
                                skymap [iz * size + ix] = picogen::Color(1.0,0.0,0.0);
                        else*/
                        
                        float r,g,b;
                        col.to_rgb (r,g,b);
                        r = r < 0.0 ? 0.0 : r > 1.0f ? 1.0f : r;
                        g = g < 0.0 ? 0.0 : g > 1.0f ? 1.0f : g;
                        b = b < 0.0 ? 0.0 : b > 1.0f ? 1.0f : b;                                        
                        const unsigned int index = makeSkymapIndex (ix, iz);
                        skymap [index][0] = r;
                        skymap [index][1] = g;
                        skymap [index][2] = b;
                }
        }
        std::cout << "done." << std::endl;
}



bool Renderer::isInitialized () const {
        return 0 != screen;
}



//
inline Renderer::Vertex1 Renderer::fun (float x, float z) const {
        
        struct { int operator () (float x) const {
                int const x_ = (int)x;
                return x_<0?0:x_>=heightfield_size?heightfield_size-1:x_;
                //return min (max(x_, 0), heightfield_size-1);
        } } wrap;
        
        struct { float operator () (float fac, float lhs, float rhs) const {
                return lhs + fac * (rhs-lhs);
        } } linear;
        
        struct { float operator () (float u, float v, float _00, float _01, float _10, float _11) const {
                float const a = (_00 + u * (_01-_00));
                float const b = (_10 + u * (_11-_10));
                return a + v * (b-a);
        } } bilinear;

        //
                
        x = x * heightfield_scale;
        z = z * heightfield_scale;
        
        int const 
                ix0 = wrap (x), 
                ix1 = wrap (x+1), 
                iz0 = wrap (z),
                iz1 = wrap (z+1)
        ;
        float const
                u = x - (float)ix0,         
                v = z - (float)iz0
        ;
        
        Vertex1 const
                _00 = heightfield [makeHeightfieldIndex (ix0, iz0)],
                _01 = heightfield [makeHeightfieldIndex (ix1, iz0)],
                _10 = heightfield [makeHeightfieldIndex (ix0, iz1)],
                _11 = heightfield [makeHeightfieldIndex (ix1, iz1)]
        ;
        
        const float
                n[3] = { 
                        bilinear (u, v, _00.normal[0], _01.normal[0], _10.normal[0], _11.normal[0]),
                        bilinear (u, v, _00.normal[1], _01.normal[1], _10.normal[1], _11.normal[1]),
                        bilinear (u, v, _00.normal[2], _01.normal[2], _10.normal[2], _11.normal[2])
                },
                ilen = 1.0f / sqrtf (n[0]*n[0] + n[1]*n[1] + n[2]*n[2]),
                normal[3] = { 
                        ilen * n[0],
                        ilen * n[1],
                        ilen * n[2]                        
                }                
        ;
        
        Vertex1 ret = {
                bilinear (u, v, _00.height, _01.height, _10.height, _11.height),
                {normal[0], normal[1], normal[2]}
        };
       
        return ret;
}



const float* Renderer::skylight (float nx, float ny, float nz) const {

        const float
                fx = 0.5 + 0.5 * nx,
                fz = 0.5 + 0.5 * nz
        ;
        const size_t
                ix = static_cast<size_t>(fx * skymap_fsize),
                iz = static_cast<size_t>(fz * skymap_fsize)
        ;
        return skymap [makeSkymapIndex (ix, iz)];
}


#ifdef SSE
#include <sse.hh>
using namespace grind::sse;


// TODO should carry a mask to see which are still valid
// TODO a) check if using a constant step yields better performance
// TODO b) check a), and check initializing the const-step based on the value d+d*x (might work as the marched distance is pretty small for each invocation of intersect())

inline Renderer::Intersection Renderer::intersect (Ray ray, Float depth, Mask active) const {

        const Float max (50.0f);        
        //const Float step (depth*0.005);
        const Float factor (0.005*1.0);

        Vertex vertex;
        
        for ( ; active.any();
                depth = depth + ternary (active, (depth*factor), Float(0.0))
                //depth = depth + ternary (active, step, Float(0.0))
        ) {
                const Vector3d position = ray (depth);
                vertex = fun (position.x, position.z);
                active = active & (position.y > vertex.height) & (depth < max);
        }
        return Intersection (ternary (depth < max, depth, Float(-1.0f)), vertex);
}



#if 1

inline Renderer::Vertex Renderer::fun (Float x, Float z) const {

        struct { Float operator () (Float x) const {                
                const Float s ((float)(heightfield_size-1));
                const Float _0 (0.0f);
                return trunc (saturate (x, _0, s));
        } } wrap;
        
        struct { float operator () (float fac, float lhs, float rhs) const {
                return lhs + fac * (rhs-lhs);
        } } linear;
        
        struct { 
                float operator () (float u, float v, float _00, float _01, float _10, float _11) const {
                        float const a = (_00 + u * (_01-_00));
                        float const b = (_10 + u * (_11-_10));
                        return a + v * (b-a);
                }
                
                Float operator () (Float u, Float v, Float _00, Float _01, Float _10, Float _11) const {
                        Float const a = (_00 + u * (_01-_00));
                        Float const b = (_10 + u * (_11-_10));
                        return a + v * (b-a);
                }
        } bilinear;
        
        //
        x = x * heightfield_scale;
        z = z * heightfield_scale;
        
        const Float s ((float)(heightfield_size-1));
        const Float _0 (0.0f);
        
        int s_ix0[4], s_ix1[4], s_iz0[4], s_iz1[4];

        for (int i=0; i<4; ++i) {
                s_ix0[i] = 0 + (int)x[i];
                s_ix1[i] = 1 + s_ix0[i];
                
                s_ix0[i] = s_ix0[i]>(heightfield_size-1)?(heightfield_size-1):s_ix0[i]<0?0:s_ix0[i];
                s_ix1[i] = s_ix1[i]>(heightfield_size-1)?(heightfield_size-1):s_ix1[i]<0?0:s_ix1[i];
                
                s_iz0[i] = 0 + (int)z[i];
                s_iz1[i] = 1 + s_iz0[i];
                
                s_iz0[i] = s_iz0[i]>(heightfield_size-1)?(heightfield_size-1):s_iz0[i]<0?0:s_iz0[i];
                s_iz1[i] = s_iz1[i]>(heightfield_size-1)?(heightfield_size-1):s_iz1[i]<0?0:s_iz1[i];
                
                //ms_iz0[i] = makeHeightfieldIndex.computeY(s_iz0[i]);
                //ms_iz1[i] = makeHeightfieldIndex.computeY(s_iz1[i]);
        }


        // get u/v
        Float const
                u = x - s_ix0,
                v = z - s_iz0
        ;
        
        Vertex1 *_00[4], *_01[4], *_10[4], *_11[4];
        
        for (size_t i=0; i<4; ++i) {
                #if 1
                const size_t x0 = makeHeightfieldIndex.computeX (s_ix0[i]);
                const size_t z0 = makeHeightfieldIndex.computeY (s_iz0[i]);
                const size_t x1 = makeHeightfieldIndex.computeX (s_ix1[i]);
                const size_t z1 = makeHeightfieldIndex.computeY (s_iz1[i]);
                _00[i] = heightfield + makeHeightfieldIndex.mergeXY (x0, z0);
                _01[i] = heightfield + makeHeightfieldIndex.mergeXY (x1, z0);
                _10[i] = heightfield + makeHeightfieldIndex.mergeXY (x0, z1);
                _11[i] = heightfield + makeHeightfieldIndex.mergeXY (x1, z1);
                #elif 0
                _00[i] = heightfield [makeHeightfieldIndex (s_ix0[i], s_iz0[i])];
                _01[i] = heightfield [makeHeightfieldIndex (s_ix1[i], s_iz0[i])];
                _10[i] = heightfield [makeHeightfieldIndex (s_ix0[i], s_iz1[i])];
                _11[i] = heightfield [makeHeightfieldIndex (s_ix1[i], s_iz1[i])];                
                #endif
        }
        _mm_empty();
        const Float height00 (_00[3]->height, _00[2]->height, _00[1]->height, _00[0]->height);
        const Float height01 (_01[3]->height, _01[2]->height, _01[1]->height, _01[0]->height);
        const Float height10 (_10[3]->height, _10[2]->height, _10[1]->height, _10[0]->height);
        const Float height11 (_11[3]->height, _11[2]->height, _11[1]->height, _11[0]->height);
        
        const Float normal00[3] = {
                Float(_00[3]->normal[0], _00[2]->normal[0], _00[1]->normal[0], _00[0]->normal[0]),
                Float(_00[3]->normal[1], _00[2]->normal[1], _00[1]->normal[1], _00[0]->normal[1]),
                Float(_00[3]->normal[2], _00[2]->normal[2], _00[1]->normal[2], _00[0]->normal[2])
        };
        const Float normal01[3] = {
                Float(_01[3]->normal[0], _01[2]->normal[0], _01[1]->normal[0], _01[0]->normal[0]),
                Float(_01[3]->normal[1], _01[2]->normal[1], _01[1]->normal[1], _01[0]->normal[1]),
                Float(_01[3]->normal[2], _01[2]->normal[2], _01[1]->normal[2], _01[0]->normal[2])
        };
        const Float normal10[3] = {
                Float(_10[3]->normal[0], _10[2]->normal[0], _10[1]->normal[0], _10[0]->normal[0]),
                Float(_10[3]->normal[1], _10[2]->normal[1], _10[1]->normal[1], _10[0]->normal[1]),
                Float(_10[3]->normal[2], _10[2]->normal[2], _10[1]->normal[2], _10[0]->normal[2])
        };
        const Float normal11[3] = {
                Float(_11[3]->normal[0], _11[2]->normal[0], _11[1]->normal[0], _11[0]->normal[0]),
                Float(_11[3]->normal[1], _11[2]->normal[1], _11[1]->normal[1], _11[0]->normal[1]),
                Float(_11[3]->normal[2], _11[2]->normal[2], _11[1]->normal[2], _11[0]->normal[2])
        };
        
       
        const Float n[3] = {
                #if 1
                bilinear (u, v, normal00[0], normal01[0], normal10[0], normal11[0]),
                bilinear (u, v, normal00[1], normal01[1], normal10[1], normal11[1]),
                bilinear (u, v, normal00[2], normal01[2], normal10[2], normal11[2])
                #elif 0
                bilinear (u, v, 
                        Float (_00[3].normal[0], _00[2].normal[0], _00[1].normal[0], _00[0].normal[0]), 
                        Float (_01[3].normal[0], _01[2].normal[0], _01[1].normal[0], _01[0].normal[0]), 
                        Float (_10[3].normal[0], _10[2].normal[0], _10[1].normal[0], _10[0].normal[0]), 
                        Float (_11[3].normal[0], _11[2].normal[0], _11[1].normal[0], _11[0].normal[0])
                ) ,
                bilinear (u, v, 
                        Float (_00[3].normal[1], _00[2].normal[1], _00[1].normal[1], _00[0].normal[1]), 
                        Float (_01[3].normal[1], _01[2].normal[1], _01[1].normal[1], _01[0].normal[1]), 
                        Float (_10[3].normal[1], _10[2].normal[1], _10[1].normal[1], _10[0].normal[1]), 
                        Float (_11[3].normal[1], _11[2].normal[1], _11[1].normal[1], _11[0].normal[1])
                ) ,
                bilinear (u, v, 
                        Float (_00[3].normal[2], _00[2].normal[2], _00[1].normal[2], _00[0].normal[2]), 
                        Float (_01[3].normal[2], _01[2].normal[2], _01[1].normal[2], _01[0].normal[2]), 
                        Float (_10[3].normal[2], _10[2].normal[2], _10[1].normal[2], _10[0].normal[2]), 
                        Float (_11[3].normal[2], _11[2].normal[2], _11[1].normal[2], _11[0].normal[2])
                )
                #elif 0
                // still fastest
                Float (
                        bilinear (u[3], v[3], _00[3].normal[0], _01[3].normal[0], _10[3].normal[0], _11[3].normal[0]),
                        bilinear (u[2], v[2], _00[2].normal[0], _01[2].normal[0], _10[2].normal[0], _11[2].normal[0]),
                        bilinear (u[1], v[1], _00[1].normal[0], _01[1].normal[0], _10[1].normal[0], _11[1].normal[0]),
                        bilinear (u[0], v[0], _00[0].normal[0], _01[0].normal[0], _10[0].normal[0], _11[0].normal[0])
                ),                                                                                    
                Float (                                                                               
                        bilinear (u[3], v[3], _00[3].normal[1], _01[3].normal[1], _10[3].normal[1], _11[3].normal[1]),
                        bilinear (u[2], v[2], _00[2].normal[1], _01[2].normal[1], _10[2].normal[1], _11[2].normal[1]),
                        bilinear (u[1], v[1], _00[1].normal[1], _01[1].normal[1], _10[1].normal[1], _11[1].normal[1]),
                        bilinear (u[0], v[0], _00[0].normal[1], _01[0].normal[1], _10[0].normal[1], _11[0].normal[1])
                ),                                                                                    
                Float (                                                                               
                        bilinear (u[3], v[3], _00[3].normal[2], _01[3].normal[2], _10[3].normal[2], _11[3].normal[2]),
                        bilinear (u[2], v[2], _00[2].normal[2], _01[2].normal[2], _10[2].normal[2], _11[2].normal[2]),
                        bilinear (u[1], v[1], _00[1].normal[2], _01[1].normal[2], _10[1].normal[2], _11[1].normal[2]),
                        bilinear (u[0], v[0], _00[0].normal[2], _01[0].normal[2], _10[0].normal[2], _11[0].normal[2])
                )
                #endif
        };

        Float ilen = rsqrt (n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
        Float normal[3] = { 
                        ilen * n[0],
                        ilen * n[1],
                        ilen * n[2]
                }
        ;
        
        //~ Vertex1 ret = {
                //~ bilinear (u, v, _00.height, _01.height, _10.height, _11.height),
                //~ {normal[0], normal[1], normal[2]}
        //~ };
       
        //~ return ret;
                
        return Vertex (
                bilinear (u, v, height00, height01, height10, height11)
                        /*(u,v,
                        Float (_00[3].height, _00[2].height, _00[1].height, _00[0].height), 
                        Float (_01[3].height, _01[2].height, _01[1].height, _01[0].height), 
                        Float (_10[3].height, _10[2].height, _10[1].height, _10[0].height), 
                        Float (_11[3].height, _11[2].height, _11[1].height, _11[0].height))*/
                ,
                Vector3d (normal[0], normal[1], normal[2])
        );
}

#else

inline Renderer::Vertex Renderer::fun (Float x, Float z) const {

        const Vertex1 v [] = {
                fun (x[3], z[3]),
                fun (x[2], z[2]),
                fun (x[1], z[1]),
                fun (x[0], z[0])
        };

        return Vertex (
                Float (v[0].height, v[1].height, v[2].height, v[3].height),
                Vector3d (
                        Float (v[0].normal[0], v[1].normal[0], v[2].normal[0], v[3].normal[0]),
                        Float (v[0].normal[1], v[1].normal[1], v[2].normal[1], v[3].normal[1]),
                        Float (v[0].normal[2], v[1].normal[2], v[2].normal[2], v[3].normal[2])
                )
        );
}

#endif


inline Ray Renderer::generateRay4 (int x, int y) const {

        return Ray (
                Vector3d (
                        Float (currentPosition[0]),
                        Float (currentPosition[1]),
                        Float (currentPosition[2])
                ),
                normalise (Vector3d (
                        Float ( (float)(3+x) / (float)width - 0.5f,
                                (float)(2+x) / (float)width - 0.5f,
                                (float)(1+x) / (float)width - 0.5f,
                                (float)(0+x) / (float)width - 0.5f),
                        Float (0.5 - (float)y / (float)height),
                        Float (0.35)
                ))
        );
}


inline void Renderer::generateRayDirection (float dir[3], int x, int y) const  {
        const float
                a = (float)x / (float)width - 0.5f,
                b = 0.5f - (float)y / (float)height,
                c = 0.35,
                rlen = 1.0f / sqrtf (a*a + b*b + c*c); 
        ;
        dir[0] = rlen * a;
        dir[1] = rlen * b;
        dir[2] = rlen * c;
}
        



#define UPWARDS 1
void Renderer::render (float px, float py, float pz) {
        if (SDL_MUSTLOCK(screen) && SDL_LockSurface(screen)<0)
                return;
        
        //const float clock = 1.0 * (float)::clock () / (float)CLOCKS_PER_SEC;
        //const float clock2 = 0.2 * (float)::clock () / (float)CLOCKS_PER_SEC;
        currentPosition [0] = px;//20.f;//30.0f * sin (clock * 0.01);
        currentPosition [1] = py;//5+5.0f * sin (clock * 0.1);
        currentPosition [2] = pz;//30.0f * sin (clock * 0.005);
        currentPosition [3] = 0.0;
        
        const Vector3d sunDirection = this->sunDirection; //sin(clock2), 0.5+sin(clock2*0.3), cos(clock2)));
        const Float sunColor[] = {
                1.0, 0.7, 0.5
        };
        
        
        static const Float
                _255(255.0f),
                _1(1.0f),
                _0(0.0f)
        ;
        
        
        Uint32 const pitch = screen->pitch>>2;
        
        const unsigned int x_step = 4;
        const unsigned int width_minus_step = width - x_step;
        
        /*if (0) {
                skylight (ray_t());
                for (unsigned int y=0; y<size; ++y) {
                for (unsigned int x=0; x<size; ++x) {
                        Uint32 *adr = ((Uint32*)screen->pixels + y*pitch + x);
                        float r,g,b;
                        skymap[y*size+x].saturate(0.0,1.0).to_rgb (r,g,b);
                        *adr = SDL_MapRGB(screen->format, (int)(r*255.0f),(int)(g*255.0f),(int)(b*255.0f));
                } }
        } else */
        #pragma omp parallel for
        for (unsigned int x=0; x<width_minus_step; x+=x_step) {

                Intersection I;
                I.depth = 3.0;
               
                unsigned int y = height-1;
                
                Uint32 *adr = ((Uint32*)screen->pixels + y*pitch + x);
                
                
                for (I = intersect (generateRay4 (x, y), I.depth, I.depth >= Float(0.0)); 
                        (y!=(unsigned int)-1) & (I.depth > Float(0.0)).any(); 
                        --y
                ) {

                        const Ray ray = generateRay4 (x, y);
                        I = intersect (ray, I.depth, I.depth >= Float(0.0));
                        
                        const Vector3d ip = ray (I.depth); 
                        
                        // NOTE i get 15% more overall performance when we call skylight() only once instead of four times
                        // 16.47 fps
                        // l2 read misses 47.78 vs. 22.7
                        
                        #if 1
                        // Ambient lighting.                        
                        const float* skycol0 = skylight (I.vertex.normal.x[0], I.vertex.normal.y[0], I.vertex.normal.z[0]);                        
                        const Float 
                                sky_r = skycol0[0]*0.7,
                                sky_g = skycol0[1]*0.7,
                                sky_b = skycol0[2]*0.7
                        ;
                        
                        
                        // Aerial perspective.
                        const float* aerial0 = skylight (ray.direction.x[0], ray.direction.y[0], ray.direction.z[0]);                        
                        const Float 
                                aerial_r = aerial0[0] * 1.1, // TODO get rid of this multiplication?
                                aerial_g = aerial0[1] * 1.1,
                                aerial_b = aerial0[2] * 1.1
                        ;
                        
                        const Float
                                dp = dot (I.vertex.normal, sunDirection),
                                f = min (I.depth*I.depth*(0.02*0.02), 1.0),
                                r_ = lerp (f, (sky_r + sunColor[0] * dp), aerial_r),
                                g_ = lerp (f, (sky_g + sunColor[1] * dp), aerial_g),
                                b_ = lerp (f, (sky_b + sunColor[2] * dp), aerial_b),
                                r = saturate (r_ * _255, _0, _255),
                                g = saturate (g_ * _255, _0, _255),
                                b = saturate (b_ * _255, _0, _255)
                        ;
                        
                        #else
                        
                        // Ambient lighting.                        
                        const float* skycol0 = skylight (I.vertex.normal.x[0], I.vertex.normal.y[0], I.vertex.normal.z[0]);
                        const float* skycol1 = skylight (I.vertex.normal.x[1], I.vertex.normal.y[1], I.vertex.normal.z[1]);
                        const float* skycol2 = skylight (I.vertex.normal.x[2], I.vertex.normal.y[2], I.vertex.normal.z[2]);
                        const float* skycol3 = skylight (I.vertex.normal.x[3], I.vertex.normal.y[3], I.vertex.normal.z[3]);
                        const Float 
                                sky_r = Float(skycol3[0], skycol2[0],skycol1[0],skycol0[0]) * 0.3,
                                sky_g = Float(skycol3[1], skycol2[1],skycol1[1],skycol0[1]) * 0.3,
                                sky_b = Float(skycol3[2], skycol2[2],skycol1[2],skycol0[2]) * 0.3
                        ;
                        
                        
                        // Aerial perspective.
                        const float* aerial0 = skylight (ray.direction.x[0], ray.direction.y[0], ray.direction.z[0]);
                        const float* aerial1 = skylight (ray.direction.x[1], ray.direction.y[1], ray.direction.z[1]);
                        const float* aerial2 = skylight (ray.direction.x[2], ray.direction.y[2], ray.direction.z[2]);
                        const float* aerial3 = skylight (ray.direction.x[3], ray.direction.y[3], ray.direction.z[3]);
                        const Float 
                                aerial_r = Float(aerial3[0], aerial2[0],aerial1[0],aerial0[0]) * 1.1,
                                aerial_g = Float(aerial3[1], aerial2[1],aerial1[1],aerial0[1]) * 1.1,
                                aerial_b = Float(aerial3[2], aerial2[2],aerial1[2],aerial0[2]) * 1.1
                        ;
                        
                        const Float
                                dp = dot (I.vertex.normal, sunDirection),
                                f = min (I.depth*I.depth*(0.02*0.02), 1.0),
                                r_ = lerp (f, (sky_r + sunColor[0] * dp), aerial_r),
                                g_ = lerp (f, (sky_g + sunColor[1] * dp), aerial_g),
                                b_ = lerp (f, (sky_b + sunColor[2] * dp), aerial_b),
                                r = saturate (r_ * _255, _0, _255),
                                g = saturate (g_ * _255, _0, _255),
                                b = saturate (b_ * _255, _0, _255)
                        ;
                        #endif
                        
                        
                        for (int i=0; i<4; ++i) {
                                if (I.depth[i]>0.0)
                                        adr[i]  = ((unsigned char)(r[i])<<16)
                                                | ((unsigned char)(g[i])<<8)
                                                | ((unsigned char)(b[i])<<0);                                         
                                // portable way: SDL_MapRGB(screen->format,r,g,b)                                           
                        }                        
                        adr -= pitch;
                }
                
                if (1) {                        
                        for ( ; y!=(unsigned int)-1; --y) {
                                //Ray ray_ = generateRay4 (x, y);
                                float dir [3];
                                generateRayDirection (dir, x, y);
                                const float* skycol = skylight (dir[0], dir[1], dir[2]);
                                const Uint32 col = ((unsigned char)(255.0f*skycol[0])<<16)
                                                 | ((unsigned char)(255.0f*skycol[1])<<8)
                                                 | ((unsigned char)(255.0f*skycol[2])<<0);
                                for (int i=0; i<4; ++i)
                                        adr[i] = col;
                                adr -= pitch;
                        }
                }
        }

        if (SDL_MUSTLOCK(screen))
                SDL_UnlockSurface(screen);
        SDL_Flip (screen);
}

#else


void Renderer::render () {
        if (SDL_MUSTLOCK(screen) && SDL_LockSurface(screen)<0)
                return;
        
        const float clock = 1.0 * (float)::clock () / (float)CLOCKS_PER_SEC;
        currentPosition [0] = 100.0f * sin (clock * 0.01);
        currentPosition [1] = 2;//5+5.0f * sin (clock * 0.1);
        currentPosition [2] = 100.0f * sin (clock * 0.005);
        currentPosition [3] = 0.0;

        
        Uint32 const pitch = screen->pitch>>2;

        for (unsigned int x=0; x!=width; ++x) {

                float depth = 0.0;
                unsigned int y = height-1;                
                
                Uint32 *adr = ((Uint32*)screen->pixels + y*pitch + x);                
                
                for (; y!=(unsigned int)-1; --y) {

                        const ray_t ray = generateRay (x, y);                        
                        depth = intersect (ray, depth);
                        
                        if (depth < 0.0) {
                                break;
                        } else {
                                const float
                                        r_ = depth * 0.25f * 0.05f,
                                        g_ = depth * 0.125f* 0.05f,
                                        b_ = depth * 0.05125f* 0.05f,
                                        r = r_ < 0.0f ? 0.0f : r_ > 1.0f ? 1.0f : r_,
                                        g = g_ < 0.0f ? 0.0f : g_ > 1.0f ? 1.0f : g_,
                                        b = b_ < 0.0f ? 0.0f : b_ > 1.0f ? 1.0f : b_
                                ;

                                *adr =
                                        SDL_MapRGB(screen->format,
                                           (unsigned char)(255.0*r),
                                           (unsigned char)(255.0*g),
                                           (unsigned char)(255.0*b)
                                          );
                        }
                        adr -= pitch;
                }                
                
                const Uint32 col = 
                                SDL_MapRGB(screen->format,
                                   (unsigned char)(255.0*0.75),
                                   (unsigned char)(255.0*0.52),
                                   (unsigned char)(255.0*0.33)
                                  );

                for ( ; y!=(unsigned int)-1; --y) {
                        *adr = col;
                        adr -= pitch;
                }
                
        }

        if (SDL_MUSTLOCK(screen))
                SDL_UnlockSurface(screen);
        SDL_Flip (screen);
}
#endif




/*

inline Renderer::ray_t Renderer::generateRay (int x, int y) const {
        const float
                position [4] = {
                        currentPosition[0],
                        currentPosition[1],
                        currentPosition[2],
                        currentPosition[3]
                },
                d_x = (float)x / (float)width - 0.5f,
                d_y = 0.5 - (float)y / (float)height,
                d_z = 0.3f,
                d_pad = 0.0f,
                d_lensq = d_x*d_x + d_y*d_y + d_z*d_z,
                d_ilen = 1.0f / sqrtf (d_lensq),
                direction [4] = {d_x * d_ilen, d_y * d_ilen, d_z * d_ilen, d_pad * d_ilen}
        ;
                
        return ray_t(position, direction);
}
*/


/*
inline float Renderer::intersect (Renderer::ray_t ray, float depth) const {
        
        float step = 0.05;
        const float epsilon = 0.125;
        for ( ; depth < 24.0; depth += (depth*2.0*0.0125f)+epsilon) {
                const float 
                        position [4] = {
                                ray.position[0]+ray.direction[0]*depth,
                                ray.position[1]+ray.direction[1]*depth,
                                ray.position[2]+ray.direction[2]*depth,
                                ray.position[3]+ray.direction[3]*depth,
                        }
                ;
                        
                const bool does = position [1] < fun (position [0], position [2]).height;

                if (does)
                        return depth;
        }
        return -1.0f;
}
*/