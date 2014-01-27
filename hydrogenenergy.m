clear all
close all
%% parameter assignment (atomic units)
M = 1837.2;
DE = 0.17639;
BETA = 1.02423;
RE = 1.40104;
RMAX = 20.0;
N = 1000;
DELTAR = RMAX/N;
%% evenly spaced grid
r = RMAX/N:DELTAR:RMAX;
% V(r)
pot = DE*(1-exp(-BETA*(r-RE))).^2-DE;
subplot(2, 2, 1);
plot(r, pot)
axis([0 15 -0.2 0.1]);
xlabel('Radial coordinate (r) / a.u')
ylabel('Morse Potential (V(r)) / a.u.')
title('Morse Potential')
%% Kinetic Energy Matrix
% Main diagonal
KE = diag(-2*ones(1, N));
% -1th diagonal
KE = KE + diag(ones(1, N-1), -1);
% 1th diagonal
KE = KE + diag(ones(1, N-1), 1);
% constants
prefactor = (-1/M)*(1/(RMAX/N)^2);
t = prefactor*KE; %%
%% Potential Energy Matrix
v = diag(pot);
h = t+v;
%% eigenfunctions
% vs vibrational quantum numbers
[s, e] = eig(h);
vs = 0:length(e)-1;
subplot(2, 2, 2);
plot(vs, diag(e), '.-')
axis([0 15 -0.2 0]);
xlabel('Vibrational Quantum Number')
ylabel('Energy / a.u.')
title('Quantum Number vs Energy')
% Normalisation choose negative sqrt? otherwise graph is upside down in
% comparison to assignment
s = -s*(1/sqrt(RMAX/N));
% Plot first 3 W.Fs fig.3 on assignment
subplot(2, 2, [3 4]);
plot(r, s(:, 1));
hold on 
plot(r, s(:, 2), '-.r');
hold on
plot(r, s(:, 3), '--g');
xlabel('Radial coordinate (r) / a.u')
ylabel('Wavefunction / a.u.')
title('Ground, first and second excited wavefunctions for Hydrogen atom')
legend('ground','first', 'second')
axis([0 3 -1.5 1.6])