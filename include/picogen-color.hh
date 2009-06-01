/***************************************************************************
 *            images.h
 *
 *  Fri Oct 12 18:53:22 2007
 *  Copyright  2007  Sebastian Mach
 *  seb@greenhybrid.net
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; version 3 of the License, or (at your
 *  option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */


#ifndef _IMAGE_H
#define _IMAGE_H

#include "picogen-geometrics.hh"

namespace picogen {

class Color {
private:
        float r,g,b;
public:
        Color (
                float r_, 
                float g_,
                float b_
        )
        : r (r_), g (g_), b (b_) {                
        }
        
        Color (float x )
        : r (x), g (x), b (x) {                
        }
        
        Color (const Color &col) : r (col.r), g (col.g), b (col.b) {                
        }
        Color() : r (0), g (0), b (0) {                
        }
        inline void from_rgb (float r_, float g_, float b_) {
                r = r_;
                g = g_;
                b = b_;
        }
        inline void from_Yxy (float Y, float x, float y) {
                float CIE_X = x * ( Y / y );
                float CIE_Y = Y;
                float CIE_Z = (1-x-y) * ( Y / y );

                /*CIE_X /= 100.0;
                CIE_Y /= 100.0;
                CIE_Z /= 100.0;*/

                r = 2.565*CIE_X -1.167*CIE_Y -0.398*CIE_Z;
                g = -1.022*CIE_X +1.978*CIE_Y +0.044*CIE_Z;
                b =  0.075*CIE_X -0.252*CIE_Y +1.177*CIE_Z;

                if ( r > 0.0031308 ) r = 1.055 * ( pow(r,(1/2.4) ) ) - 0.055;
                else                 r = 12.92 * r;
                if ( g > 0.0031308 ) g = 1.055 * ( pow(g,(1/2.4) ) ) - 0.055;
                else                 g = 12.92 * g;
                if ( b > 0.0031308 ) b = 1.055 * ( pow(b,(1/2.4) ) ) - 0.055;
                else                 b = 12.92 * b;
        }
        inline void from_rgb (const float * const rgb) {
                r = rgb[0];
                g = rgb[1];
                b = rgb[2];
        }
        inline void to_rgb (float &r_, float &g_, float &b_) const {
                r_ = r;
                g_ = g;
                b_ = b;
        }
        inline void to_rgb (float * const rgb) const {
                rgb[0] = r;
                rgb[1] = g;
                rgb[2] = b;
        }
        inline Color saturate (const Color &min, const Color &max) {
                return Color (
                               r < min.r ? min.r : r > max.r ? max.r : r,
                               g < min.g ? min.g : g > max.g ? max.g : g,
                               b < min.b ? min.b : b > max.b ? max.b : b
                       );
        }
        inline Color operator + (const Color &col) const {
                return Color (r+col.r, g+col.g, b+col.b);
        }
        inline Color operator += (const Color &col) {
                return Color (r+=col.r, g+=col.g, b+=col.b);
        }
        inline Color operator * (const Color &col) const {
                return Color (r*col.r, g*col.g, b*col.b);
        }
        inline Color operator *= (const Color &col) {
                return Color (r*=col.r, g*=col.g, b*=col.b);
        }
        inline Color operator * (float f) const {
                return Color (r*f, g*f, b*f);
        }
        inline Color operator *= (float f) {
                return Color (r*=f, g*=f, b*=f);
        }
};
class AverageColor
{
private:
        Color m_color;
        float  m_alpha;
private:
        explicit AverageColor (const Color &col, float alpha) : m_color (col), m_alpha (alpha) {}
public:
        explicit AverageColor() : m_color (0,0,0), m_alpha (0) {}
        explicit AverageColor (float r, float g, float b) : m_color (r,g,b), m_alpha (1) {}
        explicit AverageColor (const Color &col) : m_color (col), m_alpha (1) {}
        inline int getAlpha() const {
                return (int) m_alpha;
        }
        inline AverageColor operator + (const Color &A) const {
                return AverageColor (m_color+A, m_alpha+1);
        }
        inline AverageColor operator += (const Color &A) {
                return AverageColor (m_color+=A, m_alpha+=1);
        }
        inline AverageColor operator = (const Color &A) {
                return AverageColor (m_color=A, m_alpha=1);
        }
        inline operator Color () const {
                return avg();
        }
        inline Color avg() const {
                return m_alpha!=0.0 ?
                       m_color * (1.0/m_alpha) :
                       Color (0,0,0);
        }
};

}

#endif /* _IMAGE_H */
