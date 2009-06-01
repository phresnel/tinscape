/***************************************************************************
 *            preetham.h
 *
 *  Mon Oct 15 18:15:23 2007
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


#ifndef _PREETHAM_H
#define _PREETHAM_H

#include "picogen-geometrics.hh"
#include "picogen-color.hh"
#include <cmath>

namespace picogen {

class Preetham {

private:
        Vector3d m_sunDirection;
        float m_T;
        float m_sunTheta, m_sunPhi;
        float m_zenith_x, m_zenith_y, m_zenith_Y;
        float m_perez_x[5], m_perez_y[5], m_perez_Y[5];
        float m_sunSolidAngle, m_sunSolidAngleFactor;
        Color m_sunColor, m_colorFilter;
        float m_beta;

        bool m_enableSunFalloffHack;
        //float m_sunFalloffHackMaxSolidAngleFactor;
        //float m_falloffHackExponent;
        float m_falloffHackFactor [3];


        bool m_enableFogHack;
        float m_fogHackFactor, m_fogHackSatDist;

private:
        static float perez (float Theta, float gamma, const float ABCDE[]);

public:

        // -- misc ------------------------------------------------
        Preetham();
        void invalidate();
        // --------------------------------------------------------



        // -- shading / sampling ----------------------------------
        void shade (Color &color, Ray const &ray) const;
        void sunShade (Color &color, Ray const &ray) const;
        void sunSample (Color &color, Ray &ray, float &p, Vector3d const &position) const;
        void atmosphereShade (Color &color, Color const &src_color, Ray const &ray, float distance) const;
        // --------------------------------------------------------



        // -- sun direction / solid angle -------------------------
        void setSunDirection (Vector3d const &ray);
        void setSunDirection (float lat, float longi, int sm, int jd, float tOfDay);
        void setSunSolidAngleFactor (float f);
        //Vector3d getSunDirection() const;
        void getSunDirection (Vector3d &direction) const;
        float getSunArealFactor () const;
        // --------------------------------------------------------



        // -- falloff hack ----------------------------------------
        //void setSunFalloffMaxSolidAngleFactor (float f);
        //void setSunFalloffExponent (float e);
        void setSunFalloffHackParameters (float a, float b, float c);
        void enableSunFalloffHack (bool enable);
        // --------------------------------------------------------



        // -- sun color -------------------------------------------
        void setSunColor (Color col);
        void setColorFilter (Color col);
        //Color getSunColor() const;
        void getSunColor (Color &color) const;
        // --------------------------------------------------------



        // -- atmosphere ------------------------------------------
        void setTurbidity (float t);
        void enableFogHack (bool enable, float f, float satDist);
        // --------------------------------------------------------

};

}

#endif /* _PREETHAM_H */
