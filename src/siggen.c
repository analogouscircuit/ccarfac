#include <math.h>
#include <stdlib.h>
#include "siggen.h"

double * sinusoid(double dur, double fs, double f, double phi, double gain)
{
    int num_points = (int) dur*fs;
    double df = 2.0*M_PI*f/fs;
    double * vals = malloc(sizeof(double)*num_points);

    vals[0] = cos(phi);
    for(int t=1; t < num_points; t++) {
        phi += df;
        vals[t] = cos(phi);
    }

    return vals;
}

float * sinusoid_f(float dur, float fs, float f, float phi, float gain)
{
    int num_points = (int) dur*fs;
    float df = 2.0*M_PI*f/fs;
    float * vals = malloc(sizeof(float)*num_points);

    vals[0] = cos(phi);
    for(int t=1; t < num_points; t++) {
        phi += df;
        vals[t] = cos(phi);
    }

    return vals;
}

double * swept_sine(double dur, double fs, double f1, double f2, double phi, double gain, int method)
{
    int num_points = (int) dur*fs;
    double dt = 1.0/fs;
    double f = f1;
    double coef = 2*M_PI*dt;
    double * samples = malloc(sizeof(double)*num_points);
    //double samples[num_points];

    // Main Loop
    samples[0] = sin(phi);  // initialize first sample

    if(method==0) {
        // linear accumulation
        double f_delta = (f2-f1)/num_points;
        for(int k = 1; k < num_points; k++) {
            f += f_delta;
            phi += coef*f;
            samples[k] = sin(phi);
        }
    } else {
        // logarithmic accumulation
        double alpha = exp( (1.0/num_points) * log(f2/f1) );
        for(int k = 1; k < num_points; k++) {
            f *= alpha;
            phi += coef*f;
            samples[k] = sin(phi);
        }
    }
    return samples;
}


float * swept_sine_f(float dur, float fs, float f1, float f2, float phi, float gain, int method)
{
    int num_points = (int) dur*fs;
    float dt = 1.0/fs;
    float f = f1;
    float coef = 2*M_PI*dt;
    float * samples = malloc(sizeof(float)*num_points);
    //float samples[num_points];

    // Main Loop
    samples[0] = sin(phi);  // initialize first sample

    if(method==0) {
        // linear accumulation
        float f_delta = (f2-f1)/num_points;
        for(int k = 1; k < num_points; k++) {
            f += f_delta;
            phi += coef*f;
            samples[k] = sin(phi);
        }
    } else {
        // logarithmic accumulation
        float alpha = exp( (1.0/num_points) * log(f2/f1) );
        for(int k = 1; k < num_points; k++) {
            f *= alpha;
            phi += coef*f;
            samples[k] = sin(phi);
        }
    }
    return samples;
}
