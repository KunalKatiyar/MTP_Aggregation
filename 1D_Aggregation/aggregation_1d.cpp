#include "aggregation_1d_solver.h"
#include <chrono>

using namespace std;
using namespace std::chrono;

int main() {

    auto start = high_resolution_clock::now();

    double r=2; //Geometric Ratio
    int N0=1; //Initial Total Number
    double v0 = 0.01, v_max = 10000.0, v_min = 0.000005; // Initial Average Size, Maximum size/volume, Minimum size/volume
    double t_start = 0; //Time span start (for integration)
    double t_end = 10000; // Time span end (for integration)
    double a0 =1.0; // Constant Kernel
      
    // Initialize Object
    Aggregation_1D_Solver obj(r, N0, t_start, t_end, a0, v0, v_max, v_min);

    //Solve using the function
    obj.initializeSolver();

    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(stop - start);
 
    cout << "Time taken by function: " << duration.count() << " microseconds" << endl;
 
    return 0;
}