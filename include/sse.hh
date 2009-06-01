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
#ifndef SSE_HH_INCLUDED_20090507
#define SSE_HH_INCLUDED_20090507

#include <iostream>
#include <xmmintrin.h>
#include <stdint.h>

namespace grind {

        namespace sse {
                //typedef float v4sf __attribute__((vector_size (16)));
                //typedef int v2si __attribute__((vector_size (8)));
                typedef __m128 v4sf;
                typedef __m64 v2si;
                /*struct float4_ {
                        v4sf v;
                        // i am a messy
                        float4_ (v4sf const &v_) : v(v_) {}
                        float4_ (float4_ const &v_) : v(v_.v) {}
                        float4_ () {}
                };
                typedef float4_ float4 __attribute__((aligned(16)));*/



                //typedef int v4si __attribute__((vector_size(16)));
                /*struct uint4_ {
                        v4si v;

                        // ough what a mess
                        uint4_ (v4si const &v_) : v(v_) {}
                        uint4_ (v4sf const &v_) { memcpy (&v, &v_, sizeof (v4sf)); }
                        uint4_ (uint4_ const &v_) : v(v_.v) {}
                        uint4_ (uint32_t _3, uint32_t _2, uint32_t _1, uint32_t _0) {
                                const uint32_t u[4] = { _3, _2, _1, _0 };
                                memcpy (&v, u, sizeof (v));
                        }
                        uint4_ (uint32_t all) {
                                const uint32_t u[4] = { all, all, all, all };
                                memcpy (&v, u, sizeof (v));
                        }
                        uint4_ () {}
                };
                typedef uint4_ uint4 __attribute__((aligned(16)));*/



                class Float {
                public:
                        enum { width = 4 };

                        Float () {}
                        Float (float _3, float _2, float _1, float _0) : value (_mm_set_ps(_3,_2,_1,_0)) {}
                        Float (float all) : value (_mm_set1_ps(all)) {}

                        Float (v4sf const &val) : value (val) {}
                        Float (Float const &val) : value (val.value) {}
                                
                        Float (int16_t *s_) {
                                union { int16_t s[4]; __m64 m64; } const val = {{s_[0], s_[1], s_[2], s_[3]}};
                                value = _mm_cvtpi16_ps (val.m64);
                        }
                        
                        Float (int32_t *s_) {
                                union { int32_t s[4]; __m64 m64[2]; } const val = {{s_[0], s_[1], s_[2], s_[3]}};
                                value = _mm_cvtpi32x2_ps (val.m64[0], val.m64[1]);
                        }

                        Float& operator=(v4sf const &val) { value = val;        return *this; }
                        Float& operator=(Float const &val){ value = val.value;  return *this; }
                        Float& operator=(float val)       { *this = Float(val); return *this; }

                        /*operator sse::v4sf () const {
                                return value.v;
                        }*/
                        const v4sf &get () const { return value; }

                        float operator [] (size_t i) const {
                                /*float ret[4];
                                memcpy (ret, &value, 16);
                                return ret[i];*/
                                union {
                                        v4sf v;
                                        float f[4];
                                };
                                v = value;
                                return f[i];
                                //return ((float*)(&value))[i];
                        }

                        void store_aligned (float *target) const {
                                _mm_store_ps (target, value);
                        }

                        void store_unaligned (float *target) const {
                                _mm_storeu_ps (target, value);
                        }

                        bool any () const {
                                union { const v4sf v; const uint32_t u[4]; } v = { value };
                                return (bool)v.u[0]
                                        | (bool)v.u[1]
                                        | (bool)v.u[2]
                                        | (bool)v.u[3]
                                ;
                        }

                        bool all () const {
                                union { const v4sf v; const uint32_t u[4]; } v = { value };
                                return (bool)v.u[0]
                                        & (bool)v.u[1]
                                        & (bool)v.u[2]
                                        & (bool)v.u[3]
                                ;
                        }

                        bool none () const {
                                return !any();
                        }

                private:
                        v4sf value;
                };



                class Vector3d {
                public:
                        Float x, y, z;

                        Vector3d() {}
                        Vector3d(Float const &x_, Float const &y_, Float const &z_)
                        : x(x_), y(y_), z(z_) {}
                };



                class Ray {
                public:
                        Vector3d position;
                        Vector3d direction;

                        Ray () {}

                        Ray (Vector3d const &position_, Vector3d const &direction_)
                        : position(position_), direction(direction_)
                        {}

                        Vector3d operator () (Float const &length) const ;
                };


                /*
                class Mask {
                public:
                        Mask () {}
                        Mask (v4sf const &val) {
                                union { const v4sf f; const v4si i; } v = { val };
                                value = v.i;
                        }
                        Mask (v4si const &val) : value (val) {}
                        Mask (Mask const &val) : value (val.value) {}

                        //Mask& operator=(sse::v4si const &val) { value = val; return *this; }
                        Mask& operator=(Mask const &val) { value = val.value; return *this; }

                        operator sse::v4si () const {
                                return value;
                        }


                        operator sse::v4sf () const {
                                // slowest
                                //v4sf ret;
                                //memcpy (&ret, &value.v, sizeof (ret));
                                //return ret;


                                // intermediate
                                union { const v4si i; const v4sf f; } const ret = {value};
                                return ret.f;

                                //return (v4sf)value;

                                // pointer-cast is fastest, but buggy
                                //return *(v4sf*)&value;
                        }

                        bool any () const {
                                const uint32_t * const u = (uint32_t*)&value;
                                return (bool)u[0]
                                        | (bool)u[1]
                                        | (bool)u[2]
                                        | (bool)u[3]
                                ;
                        }

                        bool all () const {
                                const uint32_t * const u = (uint32_t*)&value;
                                return (bool)u[0]
                                        & (bool)u[1]
                                        & (bool)u[2]
                                        & (bool)u[3]
                                ;
                        }

                        bool none () const {
                                return !any();
                        }
                private:
                        v4si value;
                };*/

                typedef Float Mask;


                // Comparison.
                inline Mask operator < (Float const &lhs, Float const &rhs) {
                        return _mm_cmplt_ps (lhs.get(), rhs.get());
                }

                inline Mask operator <= (Float const &lhs, Float const &rhs) {
                        return _mm_cmple_ps (lhs.get(), rhs.get());
                }

                inline Mask operator > (Float const &lhs, Float const &rhs) {
                        return _mm_cmpgt_ps (lhs.get(), rhs.get());
                }

                inline Mask operator >= (Float const &lhs, Float const &rhs) {
                        return _mm_cmpge_ps (lhs.get(), rhs.get());
                }


                // result = condition ? lhs : rhs
                template <typename T>
                inline Float ternary (Mask const &condition, T const &lhs, T const &rhs) {
                        return _mm_or_ps (
                                _mm_and_ps (condition.get(), lhs.get()),
                                _mm_andnot_ps (condition.get(), rhs.get())
                        );
                }
                
                // Arithmetic.
                inline Float operator + (Float const &lhs, Float const &rhs) {
                        return _mm_add_ps (lhs.get(), rhs.get());
                }

                inline Float operator += (Float &lhs, Float const &rhs) {
                        return lhs = lhs + rhs;
                }

                inline Float operator - (Float const &lhs, Float const &rhs) {
                        return _mm_sub_ps (lhs.get(), rhs.get());
                }

                inline Float operator - (Float const &var) {
                        return _mm_sub_ps (_mm_set1_ps(0.0), var.get());
                }

                inline Float operator * (Float const &lhs, Float const &rhs) {
                        return _mm_mul_ps (lhs.get(), rhs.get());
                }

                inline Float operator / (Float const &lhs, Float const &rhs) {
                        return _mm_div_ps (lhs.get(), rhs.get());
                }

                inline Float sqrt (Float const &val) {
                        return _mm_sqrt_ps (val.get());
                }

                inline Float rsqrt (Float const &val) {
                        return _mm_rsqrt_ps (val.get());
                }

                inline Float max (Float const &lhs, Float const &rhs) {
                        return _mm_max_ps (lhs.get(), rhs.get());
                }

                inline Float min (Float const &lhs, Float const &rhs) {
                        return _mm_min_ps (lhs.get(), rhs.get());
                }

                inline Float saturate (Float const &val, Float const &min_, Float const &max_) {
                        return min (max_.get(), max (val.get(), min_.get()));
                }

                inline Float trunc (Float const &val) {
                        //return Float ((float)(int)val[3], (float)(int)val[2], (float)(int)val[1], (float)(int)val[0]);
                        return _mm_cvtpi32x2_ps (
                                 _mm_cvttps_pi32 (val.get()),
                                 _mm_cvttps_pi32 (
                                   _mm_movehl_ps (_mm_setzero_ps(), val.get())));
                }
                
                inline void round_to_short_array (Float const &val, short *vals) {
                        // __m64 _mm_cvtps_pi16( __m128 a );
                        union {
                                __m64 m64;
                                short s[4];
                        };
                        m64 = _mm_cvtps_pi16 (val.get());
                        for (int i=0; i<4; ++i) {
                                vals [i] = s[i];
                        }
                }
                
                inline void trunc_to_short_array (Float const &val, short *vals) {
                        // __m64 _mm_cvtps_pi16( __m128 a );
                        union {
                                __m64 m64;
                                short s[4];
                        };
                        m64 = _mm_cvtps_pi16 (ternary (val>0.5,val-0.5,val).get());
                        for (int i=0; i<4; ++i) {
                                vals [i] = s[i];
                        }
                }

                inline void trunc_to_int_array (Float const &val, int *vals) {
                        // __m64 _mm_cvtps_pi16( __m128 a );
                        union {
                                __m64 m64[2];
                                int i[4];
                        };
                        m64[0] = _mm_cvttps_pi32 (val.get());
                        m64[1] = _mm_cvttps_pi32 (
                                   _mm_movehl_ps (_mm_setzero_ps(), val.get()));
                        for (int u=0; u<4; ++u) {
                                vals [u] = i[u];
                        }
                }

                


                // Bitty.
                inline Mask operator & (Mask const &lhs, Mask const &rhs) {
                        return _mm_and_ps (lhs.get(), rhs.get());
                }

                inline Mask operator &= (Mask &lhs, Mask const &rhs) {
                        return lhs = lhs & rhs;
                }

                inline Mask operator | (Mask const &lhs, Mask const &rhs) {
                        return _mm_or_ps (lhs.get(), rhs.get());
                }

                inline Mask operator |= (Mask &lhs, Mask const &rhs) {
                        return lhs = lhs | rhs;
                }

                inline Mask andnot (Mask const &lhs, Mask const &rhs) {
                        return _mm_andnot_ps (lhs.get(), rhs.get());
                }


                // Sugar.                

                // result = lhs + factor * (rhs-lhs)
                //        = (1-factor) * lhs + factor * lhs
                // TODO: Note sure which is really faster, profile.
                inline Float lerp (Float const &factor, Float const &lhs, Float const &rhs) {
                        return lhs + factor * (rhs - lhs);
                }



                // Vector3d
                inline Vector3d operator+ (Vector3d const &lhs, Vector3d  const &rhs) {
                        return Vector3d (lhs.x+rhs.x, lhs.y+rhs.y, lhs.z+rhs.z);
                }

                inline Vector3d operator* (Vector3d const &vec, Float const &factor) {
                        return Vector3d (vec.x*factor, vec.y*factor, vec.z*factor);
                }

                inline Vector3d operator* (Float const &factor, Vector3d const &vec) {
                        return Vector3d (vec.x*factor, vec.y*factor, vec.z*factor);
                }

                inline Float dot (Vector3d const &lhs, Vector3d const &rhs) {
                        return lhs.x*rhs.x + lhs.y*rhs.y + lhs.z*rhs.z;
                }

                inline Vector3d normalise (Vector3d const &vec) {
                        return vec * rsqrt (vec.x*vec.x + vec.y*vec.y + vec.z*vec.z);
                }

                inline Vector3d operator - (Vector3d const &var) {
                        return Vector3d (-var.x, -var.y, -var.z);
                }


                // Ray
                inline Vector3d Ray::operator () (Float const &length) const {
                        return position + length * direction;
                }


                // Stream IO.
                inline std::ostream& operator<< (std::ostream& o, const Float& c) {
                        _mm_empty();
                        o << c[3] << ", " << c[2] << ", " << c[1] << ", " << c[0];
                        return o;
                }

                /*inline std::ostream& operator<< (std::ostream& o, const Mask& c) {
                        o << (bool)c[3] << ", " << (bool)c[2] << ", " << (bool)c[1] << ", " << (bool)c[0];
                        return o;
                }*/
        }
}
#endif // #define SSE_HH_INCLUDED_20090507
