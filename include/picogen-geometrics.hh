/***************************************************************************
 *            geometrics.h
 *
 *  Thu Oct 11 19:34:19 2007
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


#ifndef _GEOMETRICS_H
#define _GEOMETRICS_H

#include <cmath>
#include <limits>

namespace picogen {

class Vector3d {
private:        
        float m[3];
public:

        template <class T> inline float &operator [] (T u) {
                return m[u]; // do not make any type conversion, so compiler-warnings are produced
        }

        template <class T> inline float operator [] (T u) const {
                return m[u]; // do not make any type conversion, so compiler-warnings are produced
        }

        Vector3d (float X, float Y, float Z) {
                m[0] = X;
                m[1] = Y;
                m[2] = Z;
        }

        Vector3d() {
                m[0] = m[1] = m[2] = 0;
        }

        inline const Vector3d operator = (const Vector3d &v) {
                return Vector3d (m[0]=v.m[0], m[1]=v.m[1], m[2]=v.m[2]);
        }

        inline const Vector3d operator + (const Vector3d &v) const {
                return Vector3d (m[0]+v.m[0], m[1]+v.m[1], m[2]+v.m[2]);
        }

        inline const Vector3d operator - (const Vector3d &v) const {
                return Vector3d (m[0]-v.m[0], m[1]-v.m[1], m[2]-v.m[2]);
        }

        inline float operator * (const Vector3d &v) const {
                return (m[0]*v.m[0] + m[1]*v.m[1] + m[2]*v.m[2]);
        }

        inline const Vector3d operator * (const float f) const {
                return Vector3d (m[0]*f, m[1]*f, m[2]*f);
        }

        inline const Vector3d operator - () const {
                return Vector3d (-m[0], -m[1], -m[2]);
        }

        inline const Vector3d operator + () const {
                return *this;
        }

        inline const Vector3d computeCross (const Vector3d &v) const {
                return Vector3d (
                               m[1]*v.m[2] - m[2]*v.m[1],
                               m[2]*v.m[0] - m[0]*v.m[2],
                               m[0]*v.m[1] - m[1]*v.m[0]);
        }

        inline float computeLengthSq() const {
                return m[0]*m[0] + m[1]*m[1] + m[2]*m[2];
        }
        inline float computeLength() const {
                return sqrt (m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);
        }

        inline const Vector3d computeNormal() const {
                return *this * (1.0/length());
        }

        void computeCoordinateSystem (Vector3d &v1, Vector3d &v2, Vector3d &v3) const {
                v1 = computeNormal();
                // from PBRT
                if (fabs (v1.m[0]) > fabs (v1.m[1])) {
                        float ilen = 1.0 / sqrt (v1.m[0]*v1.m[0] + v1.m[2]*v1.m[2]);
                        v2 = Vector3d (-v1.m[2] * ilen, 0.0, v1.m[0] * ilen);
                } else {
                        float ilen = 1.0 / sqrt (v1.m[1]*v1.m[1] + v1.m[2]*v1.m[2]);
                        v2 = Vector3d (0.0, v1.m[2] * ilen, -v1.m[1] * ilen);
                }
                v3 = v1.computeCross (v2);
        }

        // --- below function are deprecated
        /// \todo get rid of below function
        inline Vector3d cross (const Vector3d &v) const {
                return Vector3d (
                               m[1]*v.m[2] - m[2]*v.m[1],
                               m[2]*v.m[0] - m[0]*v.m[2],
                               m[0]*v.m[1] - m[1]*v.m[0]);
        }

        /// \todo get rid of below function
        inline float lengthSq() const {
                return m[0]*m[0] + m[1]*m[1] + m[2]*m[2];
        }

        /// \todo get rid of below function
        inline float length() const {
                return sqrt (m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);
        }

        /// \todo get rid of below function
        inline Vector3d normal() const {
                return *this * (1.0/length());
        }
};



class Ray {
private:
private:
        Vector3d position, direction;
public:
        Ray() : position(), direction() {}

        Ray (const Vector3d &position_, const Vector3d &direction_)
                        : position (position_), direction (direction_) {}

        inline const Vector3d operator () (float f) const {
                return position + direction * f;
        }

        const inline Vector3d getPosition() const {
                return position;
        }

        const inline Vector3d getDirection() const {
                return direction;
        }

        void setPosition( const Vector3d &position ) {
                this->position = position;
        }

        void setDirection( const Vector3d &direction ) {
                this->direction = direction;
        }

        // --- below functions are deprecated

        /// \todo get rid of below function
        const inline Vector3d x() const {
                return position;
        }

        /// \todo get rid of below function
        const inline Vector3d w() const {
                return direction;
        }

        /// \todo get rid of below function
        inline Vector3d &x() {
                return position;
        }

        /// \todo get rid of below function
        inline Vector3d &w() {
                return direction;
        }
};


class BoundingBox {
private:
        Vector3d bbmin,bbmax;

public:
        BoundingBox() :
                bbmin (std::numeric_limits<float>::max(),std::numeric_limits<float>::max(),std::numeric_limits<float>::max()),
                bbmax (-std::numeric_limits<float>::max(),-std::numeric_limits<float>::max(),-std::numeric_limits<float>::max())
        {
        }

        BoundingBox (const Vector3d & bbmin, const Vector3d & bbmax) :
                bbmin (bbmin),
                bbmax (bbmax)
        {
        }

        /*BoundingBox (Vector3d bbmin, Vector3d bbmax) :
                bbmin (bbmin),
                bbmax (bbmax) {
        }*/

        void reset() {
                bbmin = Vector3d (std::numeric_limits<float>::max(),std::numeric_limits<float>::max(),std::numeric_limits<float>::max());
                bbmax = Vector3d (-std::numeric_limits<float>::max(),-std::numeric_limits<float>::max(),-std::numeric_limits<float>::max());
        }

        int getAxisOfMaxExtent() const {
                const Vector3d diff = bbmax - bbmin;
                int ret = 0;
                if (fabs (diff[1]) > fabs (diff[0]) && fabs (diff[1]) > fabs (diff[2])) ret = 1;
                if (fabs (diff[2]) > fabs (diff[0]) && fabs (diff[2]) > fabs (diff[1])) ret = 2;
                return ret;
        }

        Vector3d computeMedian() const {
                return (bbmin+bbmax) *0.5;
        }

        Vector3d getMin() const {
                return bbmin;
        }

        Vector3d getMax() const {
                return bbmax;
        }

        void setMin (const Vector3d &b) {
                bbmin = b;
        }

        void setMax (const Vector3d &b) {
                bbmax = b;
        }

        float computeWidth () const {
                return bbmax[0] - bbmin[0];
        }

        float computeHeight () const {
                return bbmax[1] - bbmin[1];
        }

        float computeDepth () const {
                return bbmax[2] - bbmin[2];
        }

        float computeSize (unsigned int axis) const {
                return bbmax[axis] - bbmin[axis];
        }

        Vector3d computeCenter () const {
                return (bbmin + bbmax) * 0.5;
        }

        void update (Vector3d const &x) {
                if (x[0] < bbmin[0]) bbmin[0] = x[0];
                if (x[1] < bbmin[1]) bbmin[1] = x[1];
                if (x[2] < bbmin[2]) bbmin[2] = x[2];
                if (x[0] > bbmax[0]) bbmax[0] = x[0];
                if (x[1] > bbmax[1]) bbmax[1] = x[1];
                if (x[2] > bbmax[2]) bbmax[2] = x[2];
        }

        bool intersects (Vector3d const &x) const {
                return x[0]>=bbmin[0]
                       && x[1]>=bbmin[1]
                       && x[2]>=bbmin[2]
                       && x[0]<=bbmax[0]
                       && x[1]<=bbmax[1]
                       && x[2]<=bbmax[2]
                       ;
        }

        bool intersect (float &t_min, float &t_max, Ray const &r) const {
                float tmin, tmax, tymin, tymax, tzmin, tzmax;

                float inx = 1.0 / r.w() [0];
                if (r.w() [0] >= 0) {
                        tmin = (bbmin[0] - r.x() [0]) * inx;
                        tmax = (bbmax[0] - r.x() [0]) * inx;
                } else {
                        tmin = (bbmax[0] - r.x() [0]) * inx;
                        tmax = (bbmin[0] - r.x() [0]) * inx;
                }

                float iny = 1.0 / r.w() [1];
                if (r.w() [1] >= 0) {
                        tymin = (bbmin[1] - r.x() [1]) * iny;
                        tymax = (bbmax[1] - r.x() [1]) * iny;
                } else {
                        tymin = (bbmax[1] - r.x() [1]) * iny;
                        tymax = (bbmin[1] - r.x() [1]) * iny;
                }

                if ( (tmin > tymax)  || (tymin > tmax))
                        return false;

                if (tymin > tmin)
                        tmin = tymin;
                if (tymax < tmax)
                        tmax = tymax;

                float inz = 1.0 / r.w() [2];
                if (r.w() [2] >= 0) {
                        tzmin = (bbmin[2] - r.x() [2]) * inz;
                        tzmax = (bbmax[2] - r.x() [2]) * inz;
                } else {
                        tzmin = (bbmax[2] - r.x() [2]) * inz;
                        tzmax = (bbmin[2] - r.x() [2]) * inz;
                }

                if ( (tmin > tzmax)  || (tzmin > tmax))
                        return false;
                if (tzmin > tmin)
                        tmin = tzmin;
                if (tzmax < tmax)
                        tmax = tzmax;

                t_min = tmin;
                t_max = tmax;

                return t_min>0.0001 || t_max > 0.0001;
        }
};

// #include "geometrics/Transformation.h"

}

#endif /* _GEOMETRICS_H */
