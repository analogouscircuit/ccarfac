#define LINEAR 0
#define LOGARITHMIC 1

double * sinusoid(double dur, double fs, double f, double phi, double gain);
float * sinusoid_f(float dur, float fs, float f, float phi, float gain);
double * swept_sine(double dur, double fs, double f1, double f2, double phi, double gain, int method);
float * swept_sine_f(float dur, float fs, float f1, float f2, float phi, float gain, int method);
