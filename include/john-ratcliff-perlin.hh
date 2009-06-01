// lifted and modified from John Ratcliffs perlin code on flipcode.

// made octave count a template param
// const corrected

#ifndef PERLIN_H_

#define PERLIN_H_

#include <stdlib.h>


#define SAMPLE_SIZE 1024

template <int OCTAVES> class Perlin
{
public:

  Perlin(float freq,float amp,int seed);


  float Get(float x,float y) const
  {
    float vec[2];
    vec[0] = x;
    vec[1] = y;
    return perlin_noise_2D(vec);
  };

private:
  void init_perlin(int n,float p) const;
  float perlin_noise_2D(const float vec[2]) const;

  float noise1(float arg) const ;
  float noise2(const float vec[2]) const ;
  float noise3(const float vec[3]) const ;
  void normalize2(float v[2]) const ;
  void normalize3(float v[3]) const ;

  float mFrequency;
  float mAmplitude;
  int   mSeed;

  int p[SAMPLE_SIZE + SAMPLE_SIZE + 2];
  float g3[SAMPLE_SIZE + SAMPLE_SIZE + 2][3];
  float g2[SAMPLE_SIZE + SAMPLE_SIZE + 2][2];
  float g1[SAMPLE_SIZE + SAMPLE_SIZE + 2];

};

#endif
