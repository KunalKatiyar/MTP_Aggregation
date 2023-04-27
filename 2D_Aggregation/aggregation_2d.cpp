#include "aggregation_1d_solver.h"

int main() {

    int r=2, N0=1;
    double v0 = 0.05, v_max = 10000.0, v_min = 0.00001;
    double m0 = 0.01, m_max = 10000.0, m_min = 0.00001;
    double t_start = 0;
    double t_end = 10000;
    double a0 =1.0;
      
    // Initialize Object
    Aggregation_Solver_2d obj1(r, N0, t_start, t_end, a0, v0, v_max, v_min, m_max, m_min);

    //Solve using the function
    obj1.initializeSolver();

    return 0;
}