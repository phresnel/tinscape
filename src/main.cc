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

#include <cstdlib>
#include <iostream>
#include <ctime>
#include <SDL/SDL.h>

#include <omp.h>

#include "main.hh"

namespace {

        
inline double walltime(/* double *t0 */)
{
  /*double mic, time;
  double mega = 0.000001;
  struct timeval tp;
  struct timezone tzp;
  static long base_sec = 0;
  static long base_usec = 0;

  (void) gettimeofday(&tp,&tzp);
  if (base_sec == 0)
    {
      base_sec = tp.tv_sec;
      base_usec = tp.tv_usec;
    }

  time = (double) (tp.tv_sec - base_sec);
  mic = (double) (tp.tv_usec - base_usec);
  time = (time + mic * mega);// - *t0;
  return(time);*/
        return omp_get_wtime();
}


int interactive () {
        using namespace std;
        
        Renderer renderer (512, 320);
        //Renderer renderer (1680, 1050);

        
        struct kkkk {
                int up, left, forward;                
                kkkk () :up(0), left(0), forward(0) {}
        } keystate;
        float px=0.0, py=0.0, pz=0.0;
        double lastMove = walltime();
        // bah.
        const clock_t numUpdatesPerSec = 30;        
        //const clock_t updateStep = CLOCKS_PER_SEC / numUpdatesPerSec;
        const float updateStepf = 1.0 / (double)numUpdatesPerSec;//(float)updateStep / (float)CLOCKS_PER_SEC;


        //////////////////////////////////////////////////////////////////////
        // Program main loop
        //////////////////////////////////////////////////////////////////////
        {
                // Fps.
                unsigned int fpsNumFrames = 0;
                float fpsFPS = 0.0f;
                double fpsLastMeasure = walltime();

                bool done = false;
                SDL_WM_SetCaption ("Pipapo.", "Pipapo.");
                while (!done /*&& ++f<100*/ ) {

                        // Lame fps checking.
                        if (walltime() - fpsLastMeasure > 1.0) {
                                fpsLastMeasure = walltime();
                                fpsFPS = static_cast<float>(fpsNumFrames) * 0.5;
                                fpsNumFrames = 0;
                                //sprintf( caption, "fps:%.2f", fpsFPS );
                                std::cout << fpsFPS << "fps" << std::endl;
                                //SDL_WM_SetCaption( caption, "..." ); // Can set caption with this.
                        }

                        // message processing loop
                        SDL_Event event;
                        while (SDL_PollEvent(&event)) {
                                // check for messages
                                switch (event.type) {
                                        // exit if the window is closed
                                case SDL_QUIT:
                                        done = true;
                                        break;

                                        // check for keydowns
                                case SDL_KEYDOWN: {
                                        // Exit if ESCAPE is pressed.
                                        if (event.key.keysym.sym == SDLK_ESCAPE)
                                                done = true;
                                        if (event.key.keysym.sym == SDLK_UP)
                                                keystate.up = 1;
                                        if (event.key.keysym.sym == SDLK_DOWN)
                                                keystate.up = -1;
                                        if (event.key.keysym.sym == SDLK_a)
                                                keystate.left = -1;
                                        if (event.key.keysym.sym == SDLK_d)
                                                keystate.left = 1;
                                        if (event.key.keysym.sym == SDLK_w)
                                                keystate.forward = 1;
                                        if (event.key.keysym.sym == SDLK_s)
                                                keystate.forward = -1;
                                        break;
                                }
                                // check for keyups
                                case SDL_KEYUP: {
                                        // Exit if ESCAPE is pressed.
                                        if (event.key.keysym.sym == SDLK_ESCAPE)
                                                done = true;
                                        if (event.key.keysym.sym == SDLK_UP)
                                                keystate.up = 0;
                                        if (event.key.keysym.sym == SDLK_DOWN)
                                                keystate.up = 0;
                                        if (event.key.keysym.sym == SDLK_a)
                                                keystate.left = 0;
                                        if (event.key.keysym.sym == SDLK_d)
                                                keystate.left = 0;
                                        if (event.key.keysym.sym == SDLK_w)
                                                keystate.forward = 0;
                                        if (event.key.keysym.sym == SDLK_s)
                                                keystate.forward = 0;
                                        break;
                                }


                                //~ case SDL_MOUSEMOTION: {
                                    //~ Uint8 keystate = SDL_GetMouseState(NULL, NULL);

                                    //~ if( keystate & SDL_BUTTON(1) ){
                                        //~ lumiX += 0.01*static_cast<float>(event.motion.xrel);
                                        //~ lumiY -= 0.01*static_cast<float>(event.motion.yrel);
                                    //~ }
                                //~ } break;

                                }
                        }
                        
                        const double currTime = walltime();
                        for ( ; lastMove < currTime; lastMove+=updateStepf) {
                                px += updateStepf * 5.0*(keystate.left<0 ? -1.0 : keystate.left>0 ? 1.0 : 0.0);                                
                                py += updateStepf * 5.0*(keystate.up<0 ? -1.0 : keystate.up>0 ? 1.0 : 0.0);                                
                                pz += updateStepf * 5.0*(keystate.forward<0 ? -1.0 : keystate.forward>0 ? 1.0 : 0.0);                                
                        }
                        
                        renderer.render (px, py, pz);
                        //draw (screen, voxbox);
                        ++fpsNumFrames;
                }
        }

        // All is well ;)
        cout << "exited cleanly." << endl;
        return 0;
}



template <int NUM_FRAMES> int benchmark () {
        using namespace std;
        
        const int width = 512, height = 320;
        Renderer renderer (width, height);
        renderer.render (0.0f, 0.0f, 0.0f); // hack to awake heightmap and skymap creation before actual benchmark

        
        struct kkkk {
                int up, left, forward;                
                kkkk () :up(0), left(0), forward(0) {}
        } keystate;
        float px=0.0, py=0.0, pz=0.0;
        
        // bah.
        double lastMove = walltime();
        // bah.
        const clock_t numUpdatesPerSec = 30;        
        //const clock_t updateStep = CLOCKS_PER_SEC / numUpdatesPerSec;
        const float updateStepf = 1.0 / (double)numUpdatesPerSec;//(float)updateStep / (float)CLOCKS_PER_SEC;


        //////////////////////////////////////////////////////////////////////
        // Program main loop
        //////////////////////////////////////////////////////////////////////
        {
                // Fps.
                unsigned int numFrames = 0;               
                float ftime = 0.0f;

                bool done = false;
                SDL_WM_SetCaption ("voxy benchmark.", "voxy benchmark.");
                
                const double startTime = walltime();
                while (numFrames < NUM_FRAMES) {
                        px = 10 + 10.0f * sinf (ftime);
                        py =  5 +  3.0f * sinf (0.3f * ftime);
                        pz = 10 + 10.0f * sinf (0.7f * ftime);
                        ftime += 0.1f;
                        
                        renderer.render (px, py, pz);
                        ++numFrames;
                }
                const double endTime = walltime();
                const float
                        time = (float)(endTime-startTime),
                        fps = (float)numFrames/time
                ;
                cout << "completed " << numFrames << " frames in " 
                     << time << " seconds == " << fps << "fps @ " << width << "x" << height << "\n";
        }

        // All is well ;)
        cout << "exited cleanly." << endl;
        return 0;
}


}



///////////////////////////////////////////////////////////////////////////////
// Main
///////////////////////////////////////////////////////////////////////////////
int main (int argc, char** argv) {
        --argc; ++argv;
        for (int i=0; i<argc; ++i) {
                if (!strcmp (argv[i], "-b"))
                        benchmark<1000>();
                else if (!strcmp (argv[i], "-p"))
                        benchmark<10>();
                else if (!strcmp (argv[i], "-i"))
                        interactive();
        }
}
