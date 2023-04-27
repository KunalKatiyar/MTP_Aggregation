#include <iostream>
#include <cmath>

using namespace std;

const int N = 100; // number of bins
const double L = 10.0; // total volume
const double dt = 0.01; // time step
const double alpha = 0.3; // aggregation rate constant
const double beta = 0.2; // breakup rate constant
const double k = 1.0; // constant

// calculate the derivative of the population in each bin
void deriv(double t, double y[], double dydt[]) {
    double dp = L / N; // bin size
    for (int i = 0; i < N; i++) {
        double p = (i + 0.5) * dp; // bin center
        dydt[i] = k * (y[i-1] - y[i]) / dp + alpha * (y[i-1] * y[N-i-1] - y[i] * y[N-i-2]) + beta * (y[i+1] + y[N-i-3] - 2.0 * y[i]);
    }
}

int main() {
    double t = 0.0; // initial time
    double y[N]; // population in each bin
    double dydt[N]; // derivative of population in each bin

    // set initial conditions
    for (int i = 0; i < N; i++) {
        double p = (i + 0.5) * L / N; // bin center
        y[i] = exp(-pow(p, 2)); // Gaussian distribution
    }

    // solve the population balance equation using the fourth-order Runge-Kutta method
    while (t < 10.0) {
        deriv(t, y, dydt);
        double k1[N], k2[N], k3[N], k4[N];
        for (int i = 0; i < N; i++) {
            k1[i] = dt * dydt[i];
            k2[i] = dt * (dydt[i] + 0.5 * k1[i]);
            k3[i] = dt * (dydt[i] + 0.5 * k2[i]);
            k4[i] = dt * (dydt[i] + k3[i]);
            y[i] += (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
        }
        t += dt;
    }

    // print the final population distribution
    for (int i = 0; i < N; i++) {
        double p = (i + 0.5) * L / N; // bin center
        cout << "p = " << p << ", y = " << y[i] << endl;
    }

    return 0;
}