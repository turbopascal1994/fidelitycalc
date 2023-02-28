% Just calculating the sequence fidelity
clear;
tstep = 1e-15; % Schrodinger equation step

% Qubit frequencies
w01 = 4.73047*2*pi*1e9;
w12 = w01 - 0.25*2*pi*1e9; 
wt = 25*2*pi*1e9; % Generator frequency
w = 4*1e-12; % Pulse length
T = 2*pi/wt; % Period
Theta = 0.032; % Angle
NumberofCycles = 8;
line = '1110111100011000100001000110001111011';
len = length(line);
% Convert string into sequence array
Sequence = EnterSeq(line,NumberofCycles);
L = Leakage(w01, w12, w, T, tstep, Sequence, Theta);
F = Fidelity(w01, w12, w, T, tstep, Sequence, Theta);
disp(['W2 = ', num2str(L)]);
disp(['F = ', num2str(1 - F)]);

