#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265358979323846

/*  butterworth filter coeffients */
//generated using MATLAB
double b0 = 0.000098259;
double b1 = 0.0000019652;
double b2 = 0.00000098259;
double a0 = 1.0000;
double a1 = -1.9911;
double a2 = 0.9912;

// Butterworth filter coeffients generate function 
void butterworth(int n, double fc, int fs, double *b, double *a) {
    double c = 1 / tan(PI * fc / fs);
    double a0 = 1 + sqrt(2) * c + c * c;
    b[0] = 1 / a0;
    b[1] = 2 * b[0];
    b[2] = b[0];
    a[0] = 1;
    a[1] = 2 * (c * c - 1) / a0;
    a[2] = (1 - sqrt(2) * c + c * c) / a0;
}

// Filter function
void filter(double *x, double *y, int len, double *b, double *a) {
    y[0] = b0 * x[0];
    y[1] = b0 * x[1] + b1 * x[0] - a1 * y[0];
    for (int i = 2; i < len; i++) {
        y[i] = b0 * x[i] + b1 * x[i - 1] + b2 * x[i - 2] - a1 * y[i - 1] - a2 * y[i - 2];
    }
}

// amplitude (i_phase ^2 + q_phase ^2)
void amplitude(double* m, double* n, int len){
    double amp[len];
    for(int i;i<len;i++){
        // sqrt of i_phase ^2 + q_phase ^2
        amp[i]=sqrt(pow(m[i],2)+pow(n[i],2));
        printf("amplitude = %lf\n",amp[i]);
}
}

int main() {

    // Define time and frequency parameters
    int fs = 100000; // Sampling rate
    double t[fs];
    double f1 = 6000; // Frequency of 6kHz wave

    // Create time vector
    for (int i = 0; i < fs; i++) {
        t[i] = i / (double)fs;
    }

    // Create signals
    double x[fs];
    double s1[fs];
    double s2[fs];

    for (int i = 0; i < fs; i++) {
        x[i]  = 3  * sin(2 * PI * f1 * t[i]); // 6kHz wave
        s1[i] = 10 * cos(2 * PI * f1 * t[i]); // Sine wave
        s2[i] = 10 * sin(2 * PI * f1 * t[i]); // Sine wave
    }

    // Multiply signals
    double y1[fs];
    double y2[fs];


    for (int i = 0; i < fs; i++) {
        y1[i] = x[i] * s2[i]; // Sine multiplication
    }

    for (int i = 0; i < fs; i++) {
        y2[i] = x[i] * s1[i]; // Sine multiplication
    }

    // Define filter parameters
    double f_cutoff = 100; // Cutoff frequency
    double b[3], a[3];
    
   // butterworth(2, f_cutoff, fs, b, a); // Call butterworth function

    // Filter the signal in phase
    double yfs[fs];
    double yfc[fs];
    filter(y1, yfs, fs, b, a); // Call filter function
    filter(y2, yfc, fs, b, a);
    amplitude(yfs,yfc,fs);
    return 0;
}