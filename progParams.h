// This is start of the header guard.  initial_H can be any unique name.  By convention, we use the name of the header file.
#ifndef __PROGPARAMS_H__
#define __PROGPARAMS_H__

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cstring>
#include <math.h>

#define DEBUG 0

// This is the content of the .h file, which is where the declarations go

class programParams
{
public:
    programParams();
    int
    initialize(char * filename);
    void
    announce();

    int matrixsize;
    int stepnumber;
    char * finfile;
    char * outfile;
    int initialconfig;
    int measure;
    char * inifile;
    double gD2;
    double gD22;
    double gD4;
    double wmoveA;
    int Type;
    int moveT;
};


// This is the end of the header guard
#endif // ifndef __PROGPARAMS_H__
