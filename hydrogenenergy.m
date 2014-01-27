%% parameter assignment (atomic units)
M = 1837.2;
DE = 0.17639;
BETA = 1.02423;
RE = 1.40104;
RMAX = 20.0;
N = 1000;
%% evenly spaced grid
r = RMAX/N:RMAX/N:RMAX;
% V(r)
pot = DE*(1-exp(-BETA*((r-RE).^2)))-DE;
plot(r, pot)
axis([0 15 -0.2 0.1]);
%% Kinetic Energy Matrix
% Main diagonal
KE = diag(-2*ones(1, N));
% -1th diagonal
KE = KE + diag(ones(1, N-1), -1);
% 1th diagonal
KE = KE + diag(ones(1, N-1), 1);
% constants
t = KE*(-1/(M*(RMAX/N)^2));
%% Potential Energy Matrix
v = diag(pot);
h = t+v;
