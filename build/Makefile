#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Copyright (C) 2009  Sebastian Mach (*1983)
# * phresnel/at/gmail/dot/com
# * http:#phresnel.org
# * http:#picogen.org
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DIR_LIB := ./lib
DIR_OBJ := ./obj
DIR_BIN := ./bin

LIBS_SDL := `sdl-config --static-libs`
LIBS := $(LIBS_SDL) -lgomp
INCLUDEDIRS := -I../include


#ifneq (,$(findstring profile,$(MAKEFLAGS)))
        #echo "profile build"
#else
#endif

ifeq ($(debug),yes)
	CXXFLAGS += -g -O0 -msse -fno-inline
else
	CXXFLAGS += -O3 -march=native -funroll-loops
	CXXFLAGS += -ffast-math -fno-math-errno -fno-trapping-math -fsingle-precision-constant
	CXXFLAGS += -fno-signed-zeros -freciprocal-math -march=native -fassociative-math 
	CXXFLAGS += -mfpmath=sse
	CXXFLAGS += -ftree-vectorize
	CXXFLAGS += -mrecip -funsafe-math-optimizations -ffinite-math-only
	CXXFLAGS += -fno-rtti -fno-exceptions
	CXXFLAGS += -fomit-frame-pointer
	CXXFLAGS += -fstrict-aliasing -fargument-noalias-anything
	CXXFLAGS += -funsafe-loop-optimizations 
       
        CXXFLAGS += -fopenmp
	#CXXFLAGS += -g -ggdb -O0 -msse -fno-inline

	#CXXFLAGS += -fverbose-asm -S
	#CXXFLAGS += -fbranch-probabilities
	#CXXFLAGS += -ftest-coverage -fprofile-arcs -fvpt -fprofile-values
	#CXXFLAGS += -ggdb3 -fprofile-arcs -ftest-coverage
	# add --coverage to link and compiöe
        ifeq ($(profile),yes)
                CXXFLAGS += -g
        endif
endif

CXXFLAGS += $(INCLUDEDIRS)
CXXFLAGS += `sdl-config --cflags` 
CXXFLAGS += -Wall -Wextra -fdiagnostics-show-option -Wno-unused -pedantic
        

LINK := g++-4.4 $(LIBS)
AR   := ar -rs 
CXX  := g++-4.4 -c $(CXXFLAGS) $(LIBS)

AMALGAM := g++-4.4 $(CXXFLAGS) \
	-combine -fwhole-program -funit-at-a-time -DAMALGAM \
	 -fipa-cp -fipa-matrix-reorg  \
	$(LIBS)

OBJFILES :=
OBJFILES += $(DIR_OBJ)/main.o
OBJFILES += $(DIR_OBJ)/renderer.o
OBJFILES += $(DIR_OBJ)/preetham.o

.PHONY: all
all: voxbox amalgam

.PHONY: voxbox
voxbox: $(OBJFILES)
	$(LINK) -o $(DIR_BIN)/voxbox $^

amalgam: ../src/*.cc
	cat ../src/*.cc > amalgam.cc 
	$(AMALGAM) -o $(DIR_BIN)/voxbox-amalgam amalgam.cc

#.PHONY: profile
#profile: $(OBJFILES)
#	echo $(CXXFLAGS)
#	echo $@
#	echo $(target_name)
	#$(LINK) -o $(DIR_BIN)/voxbox-profile $^



$(DIR_OBJ)/main.o: ../src/main.cc
	$(CXX) -o $@ $<	

$(DIR_OBJ)/renderer.o: ../src/renderer.cc
	$(CXX) -o $@ $<	

$(DIR_OBJ)/preetham.o: ../src/preetham.cc
	$(CXX) -o $@ $<	


.PHONY: clean
clean:
	rm -f bin/*
	rm -f obj/*
	rm -f lib/*
