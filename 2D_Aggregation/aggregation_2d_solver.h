#ifndef AGGREGATION_2D_SOLVER_H
#define AGGREGATION_2D_SOLVER_H  


class Aggregation_2D_Solver
{
    private:
        int r, N0;
        double t_start, t_end;
        double v0, v_max, v_min, a0, m0, m_max, m_min;
        vector<double> tspan;

        int m,mv,n,nv;
        vector<double> vi,wi,xi,yi,wyi;
        vector<vector<double>> Nij;

        vector<vector<double>> diag; 
        vector<vector<vector<double>>> eta;

        int nt = 45;

        vector<vector<double>> finaly; 
        vector<double> time;
        vector<double> Ntol;

        Doub h1 = 2;
        Doub hmin = 0.045500;
        

        Doub atol = 1.0e-6; // absolute tolerance
        Doub rtol = 1.0e-3; // relative tolerance

        FILE *fp = NULL;
        FILE *fp2 = NULL;
        FILE *fp3 = NULL;
        FILE *fp4 = NULL;
        FILE *fp5 = NULL;
    public:

        Aggregation_2D_Solver(int r, int N0, double t_start, double t_end, double a0, double v0, double v_max, double v_min, double m0, double m_max, double m_min);

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