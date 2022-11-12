#include "aggregation_1d_solver.h"
#include <stdio.h>
#include <iostream>
#include <string>

using namespace std;


double Aggregation_1D_Solver::round(double var)
{
    // double value = (int)(var * 10000 + .5);
    // return (double)value / 10000;
    return var;
}

int Aggregation_1D_Solver::factorial(int n)    
{    
    if(n<0)    
        return -1;   
    if(n==0)    
        return 1 ;  
    else {    
        return n*factorial(n-1);        
    }    
}  


struct Aggregation_1D_Solver::rhs_van
{
    int len;

    vector<vector<vector<double>>> eta;
    vector<vector<double>> diag;

    rhs_van(int len, vector<vector<vector<double>>> eta, vector<vector<double>> diag) : len(len),eta(eta), diag(diag) {}

    void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx)
    {
        for(int i=0;i<len;i++){
            double birth=0.0;
            for(int j=0;j<len;j++){
                for(int k=j;k<len;k++){
                    birth+=((1-0.5*diag[j][k])*eta[i][j][k]*y[j]*y[k]);
                }
            }
            double sum=0.0;
            for(int p=0;p<len;p++){
                sum+=y[p];
            }
            dydx[i]=(-1.0)*y[i]*sum + birth;
        }
    }
};

Aggregation_1D_Solver::Aggregation_1D_Solver(double r, int N0, double t_start, double t_end, double a0, double v0, double v_max, double v_min):
    r(r), N0(N0), t_start(t_start), t_end(t_end), a0(a0), v0(v0), v_max(v_max), v_min(v_min) {}

void Aggregation_1D_Solver::geometricGrid(){
    tspan.push_back(t_start);
    tspan.push_back(t_end);
    m = ceil(log(v_max/v_min)/log(r));
    mv = m+1;
    for(int i=0;i<mv;i++){
        vi.push_back(0);
    }
    vi[0] = v_min;
    for(int i=1;i<mv;i++){
        vi[i] = vi[i-1]*r;
    }
    for(int i=0;i<m;i++){
        xi.push_back(0);
    }
    for(int i=0;i<m;i++){
        wi.push_back(0);
    }
    for(int i=0;i<m;i++){
        xi[i] = (vi[i]+vi[i+1])/2;
        wi[i] = vi[i+1] - vi[i];
    }
}

void Aggregation_1D_Solver::initialConditions() {
    for(int i=0;i<m;i++){
        Ni.push_back(0);
    }
    // for(int i=0;i<m;i++){
    //     Ni[i]= (N0)*(exp(-vi[i]/v0)-exp(-vi[i+1]/v0));
    //     Ni[i] = round(Ni[i]);
    // }
    for(int i=0;i<m;i++){
        Ni[i]= N0*(exp(-vi[i]/v0)-exp(-vi[i+1]/v0));
        Ni[i] = round(Ni[i]);
    }
    for(int i=0;i<m;i++){
        double val = (N0/v0)*(vi[i]/v0)*exp(-vi[i]/v0);
        fprintf(fpnew,"%lf %e\n", xi[i], val);
    }
}

void Aggregation_1D_Solver::generatingETA() {
    diag.resize(m,vector<double>(m));

    for(int i = 0; i < m; i++)
    {
        for(int j = 0; j < m; j++)
        {
            if(i==j){
                diag[i][j] =1.0;
            }
        }
    }

    eta.resize(m,vector<vector<double> >(m,vector<double>(m)));
    
    for(int j=0;j<m;j++){
        for(int k=j;k<m;k++){
            double v = xi[j] + xi[k];
            for(int i=0;i<m-1;i++){
                if(v>xi[i] && v<=xi[i+1]){
                    eta[i+1][j][k] = (v-xi[i])/(xi[i+1]-xi[i]);
                    eta[i][j][k] = (xi[i+1]-v)/(xi[i+1]-xi[i]);
                    break;
                }
            }
        }
    }
}

void Aggregation_1D_Solver::solvingODE () {

    VecDoub ystart(m);
    for(int i=0;i<m;i++){
        ystart[i]=Ni[i];
    }

    Output out(nt-1);

    rhs_van d(m, eta, diag);

    Odeint<StepperDopr5<rhs_van> > ode(ystart, tspan[0], tspan[1], atol, rtol, h1, hmin, out, d);

    ode.integrate();

    finaly.resize(nt,vector<double>(m));

    for(int i=0;i<nt;i++){
        out.xsave[i]=round(out.xsave[i]);
        time.push_back(out.xsave[i]);
        for(int j=0;j<m;j++){
            out.ysave[j][i]=round(out.ysave[j][i]);
            finaly[i][j] = out.ysave[j][i];
        }
    }

    
    for(int i=0;i<nt;i++){
        double val=0.0;
        for(int j=0;j<m;j++){
            val+=finaly[i][j];
        }
        Ntol.push_back(val);
    }
}

void Aggregation_1D_Solver::storingData () {

    // Figure 1

    for(int i=0;i<nt;i++){
        fprintf(fp,"%e %e\n", time[i], Ntol[i]);
    }

    vector<double> Ntol_ana;
    for(int i=0;i<nt;i++){
        Ntol_ana.push_back(1/(1/N0+a0*time[i]/2));
        Ntol_ana[i]=round(Ntol_ana[i]);
    }

    for(int i=0;i<nt;i++){
        fprintf(fp2,"%llf %llf\n", time[i], Ntol_ana[i]);
    }


    // Figure 2

    vector<vector<double>> numden( nt , vector<double> (m, 0)); 
    for(int i=0;i<nt;i++){
        for(int j=0;j<m;j++){
            numden[i][j] = finaly[i][j]/wi[j];
        }
    }


    for(int i=0;i<m;i++){
        cout << xi[i] << " " << numden[nt-1][i] << "\n";
        cout << time[11] << "\n";
        fprintf(fp3,"%lf %e\n", xi[i], numden[11][i]);
    }

    vector<vector<double>> anaden( nt , vector<double> (m, 0)); 
    for(int i=0;i<nt;i++){
        double tau = a0*N0*time[i];
        for(int j=0;j<m;j++){
            double v = xi[j];
            double sumn = 0;
            for(int k=1;k<11;k++){
                sumn+=(pow((tau/(tau+2)),(k+1))*pow(4,k)*pow((v/v0),(2*k-2))/factorial(2*k+1));
            }
            anaden[i][j]=(4*N0/(pow((tau+2),2)))*(4/v0)*(exp(-2*v/v0))*sumn;
        }
    }

    for(int i=0;i<m;i++){
        fprintf(fp4,"%lf %lf\n", xi[i], anaden[0][i]);
    }

    for(int i=0;i<m;i++){
        fprintf(fp5,"%lf %lf\n", xi[i], anaden[nt-1][i]);
    }
}


//Initializes Solver and uses different functions to generate result
void Aggregation_1D_Solver::initializeSolver(){
    
    fp = fopen("Ntol.tmp", "w");
    fpnew = fopen("Nneww.tmp", "w");
    fp2 = fopen("Ntol_ana.tmp", "w"); 
    fp3 = fopen("numden2_new.tmp", "w");
    fp4 = fopen("anaden_ini.tmp", "w");
    fp5 = fopen("anaden_final.tmp", "w"); 


    cout << "Creating the geometric grid" << "\n";
    geometricGrid();

    cout << "Setting up Intitial Conditions" << "\n";
    initialConditions();

    cout << "Generating eta" << "\n";
    generatingETA();

    cout << "Solving the ODE" << "\n";
    solvingODE();

    cout << "Storing Data" << "\n";
    storingData();

    return;
}