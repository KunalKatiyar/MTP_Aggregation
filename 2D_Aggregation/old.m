function [] = aggregation_code()
clear all;
close all;
clc;
% ======== User handle starts ==================
%% Initializing
v_max = 100.0; % Maximum size/volume
v_min =1e-5; % Minimum size/volume
r = 2; % Geometric ratio
v0 = 0.01; % Initial average size
N0 = 1; % Initial total number
tspan = [0 2]; % time span for integration
%% Generating discretized grid
disp('Creating the geometric grid..')
m = ceil(log(v_max/v_min)/log(r)); % No pivots
mv = m + 1; % Number of volume boundaries
vi = zeros(mv,1);
vi(1)=v_min;
for i=2:mv
vi(i)=vi(i-1)*r;
end
x=zeros(m,1); % Creating pivots
w=zeros(m,1); % Bin width
for i = 1:m
x(i) = (vi(i)+vi(i+1))/2;
w(i) = vi(i+1) - vi(i);
end

disp('Setting-up initial conditions...')
Ni0 = zeros(m,1);
for i = 1:m
Ni0(i) = N0*(exp(-vi(i)/v0)-exp(-vi(i+1)/v0));
end
%%
%% Aggregation Kernel
% Constant Kernel
a0 = 1.0;
%% ======== User handle Ends ==================
%% ======== \eta_ijk calculation ==============
% n(i,j,k) is contribution to ith pivot when
% particles on jth and kth pivot aggregate.
% For details see Kumar and Ramkrishna 1996
disp('Generating \eta_ijk..')
eta = zeros(m,m,m);
for j = 1:m
for k = j:m
v = x(j) + x(k);
for i = 1:m-1
if(v>x(i) && v<= x(i+1))
eta(i+1,j,k) = (v-x(i))/(x(i+1)-x(i));
eta(i,j,k)= (x(i+1)-v)/(x(i+1)-x(i));
break;
end
end
end
end
%% =============================================
%% solver
disp('Solving the ODE..')
tic
[t,N] = ode45(@pbe,tspan,Ni0);
toc
%%
%% ========= Discretized equations =============
function dN_dt = pbe(t,N)
dN_dt = zeros(m,1);
dgnl = eye(m); %This works as Kronecker delta
for i = 1:m
birth = 0;
for j = 1:m
for k = j:m
birth = birth + ...
(1-0.5*dgnl(j,k))*eta(i,j,k)*N(j)*N(k);
end
end
dN_dt(i) = - N(i)*sum(N)+ birth;
end
end
%% =============================================
%% Presentation of data
% Plotting and comparing zeroth moment. Closed
% form analytical expression is available for
% zeroth moment for constant kernel.
nt = size(t,1);
figure(1)
disp('Comparing zeroth moment with analytical.')
Ntol = sum(N,2);
loglog(t,Ntol,'s')
xlabel('time')
ylabel('N(t)/N(0)')
hold on
Ntot_ana = zeros(nt,1);
for i=1:nt
Ntot_ana(i) = 1/(1/N0+a0*t(i)/2);
end
loglog(t,Ntot_ana)
legend('Numerical','Analytical')
title('Comparison with analytical zeroth moment')
%
%% Plotting and comparing number density
figure(2)
% Convirt the Ni to density
numden=zeros(nt,m); % Numerical numberdensity
for i=1:nt
for j=1:m
numden(i,j) = N(i,j)/w(j);
end
end
% Plot numerical number density
loglog(x,numden(nt,:),'*')
axis([1e-4 1e1 1e-1 1e3])
hold on
% Calculate the analytical density.
% See Chakraborty and Kumar, Chemical
% Engineering Science 2007 for the
% expression and reference for the same.
anaden=zeros(nt,m);
for i =1:nt
tau = a0*N0*t(i);
for j=1:m
v=x(j);
sumn = 0;
for k=1:12
sumn = sumn + ...
((tau/(tau+2))^(k-1))*((v/v0)^(k-1))/factorial(k-1);
end
anaden(i,j)=(4*N0/((tau+2)^2))*(1/v0)*exp(-v/v0)*sumn;
end
end
% Plot the initial condition.
loglog(x,anaden(1,:))
% Plot the analytical number density
loglog(x,anaden(nt,:))
legend('Numerical','Initial','Analytical')
title('Comparison with analytical number density')
end