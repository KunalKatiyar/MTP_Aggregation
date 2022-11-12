#ifndef AGGREGATION_1D_SOLVER_H
#define AGGREGATION_1D_SOLVER_H  


#include "nr3.h"
#include "stepper.h"
#include "stepperdopr5.h"
#include "odeint.h"
#include <stdio.h>
#include <iostream>
#include <string>

class Aggregation_1D_Solver
{
    private:
        double r;
        int N0;
        double t_start, t_end;
        double v0, v_max, v_min, a0;
        vector<double> tspan;

        int m,mv;
        vector<double> vi,wi,xi,Ni;

        vector<vector<double>> diag; 
        vector<vector<vector<double>>> eta;

        

        vector<vector<double>> finaly; 
        vector<double> time;
        vector<double> Ntol;



        int nt = 45;
        Doub h1 = 2;
        Doub hmin = 0.045500;
        
        Doub atol = 1.0e-6; // absolute tolerance
        Doub rtol = 1.0e-3; // relative tolerance

        FILE *fp = NULL;
        FILE *fpnew = NULL;
        FILE *fp2 = NULL;
        FILE *fp3 = NULL;
        FILE *fp4 = NULL;
        FILE *fp5 = NULL;

    public:
        Aggregation_1D_Solver(double r, int N0, double t_start, double t_end, double a0, double v0, double v_max, double v_min);

        double round(double var);

        int factorial(int n);

        struct rhs_van;

        void geometricGrid();

        void initialConditions();

        void generatingETA();

        void solvingODE ();

        void storingData ();

        void initializeSolver();
};

#endif