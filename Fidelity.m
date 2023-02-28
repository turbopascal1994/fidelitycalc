function F = ...
    Fidelity(w01, w12, w, T, tstep, SignalString, Theta)

CellsNumber = length(SignalString);
h = 1.054e-34; % Planck constant
C1 = 1e-12; % Qubit capacitance
F0 = 2.06e-15; % Magnetic flux quantum
Id = [1 0 0; 0 1 0; 0 0 1]; % Identity matrix
V = F0/w; % Voltage
Cc = Theta/(F0*sqrt(2*w01/(h*C1))); % Connection capacitance
Amp = Cc*V*sqrt(h*w01/(2*C1)); % Pulse amplitude
    
% Initial conditions
% |z+> = [0 1 0]
% |z-> = [1 0 0]
% |x+> = [1/sqrt(2)  1/sqrt(2)  0]
% |x-> = [1/sqrt(2) -1/sqrt(2)  0]
% |y+> = [1/sqrt(2)  1i/sqrt(2) 0]
% |y-> = [1/sqrt(2) -1i/sqrt(2) 0]
InitStates = ...
   [1 0 1/sqrt(2) 1/sqrt(2) 1/sqrt(2) 1/sqrt(2);...
    0 1 1/sqrt(2) -1/sqrt(2) 1i/sqrt(2) -1i/sqrt(2);...
    0 0 0 0 0 0];

H0 = [0 0 0; 0 h*w01 0; 0 0 h*w01+h*w12]; % Undisturbed hamiltonian
[EigVec, EigVal] = eig(H0); % Finding eigenvectors and eigenvalues
[~, ind] = sort(diag(EigVal)); % Sorting eigenvalues
EigVals = EigVal(ind,ind);
EigVecs = EigVec(:,ind); % Sorting eigenvectors

WF0 = EigVec(:,1); % State |0>
WF1 = EigVec(:,2); % State |1>
WF2 = EigVec(:,3); % State |2>

% Ideal gate matrix
Y = 1/sqrt(2)*(WF1*ctranspose(WF1) + WF0*ctranspose(WF0) + ...
    WF1*ctranspose(WF0) - WF0*ctranspose(WF1));

% Hamiltonians
Hmatrix = [0 -1 0; 1 0 -sqrt(2); 0 sqrt(2) 0]; % Disturbed hamiltonian
HrPlus = 1i*Amp*Hmatrix + H0;
HrMinus = -1i*Amp*Hmatrix + H0;
HrZero = H0;

% Evolution operators (EO)
% Umartix function turns Hamiltonian into the evolution operator
UPlusStep = UMatrix(HrPlus,tstep); % EO for the grid step
UMinusStep = UMatrix(HrMinus,tstep);
UZeroStep = UMatrix(HrZero,tstep);
UPlus = Id;
UMinus = Id;
UZero = Id;
UT = Id;

% Constructing EO for the whole pulse
for i = 1:1:(w/tstep)
    UPlus = UPlusStep*UPlus;
    UMinus = UMinusStep*UMinus;
    UZero = UZeroStep*UZero;
end
UTStep = UMatrix(HrZero,tstep);

% EO for the space between pulses
for i = 1:1:((T - w)/tstep)
    UT = UTStep*UT;
end

% Combining EOs
UTPlus = UT*UPlus; % 1 pulse + free precession
UTMinus = UT*UMinus; % -1 pulse + free precession
UTZero = UT*UZero; % 0 pulse + free precession
Umatrix = Id;

for j = 1:1:CellsNumber
	switch SignalString(j)
        case 1
            U = UTPlus;
        case -1
            U = UTMinus;
        case 0
            U = UTZero;
	end
	Umatrix = U*Umatrix; % Applying evolution operator
end

% Fidelity, method |<PSI_ideal|PSI>|^2
F = 0;
for IS = 1:1:6
    alpha = InitStates(:,IS);
    ket = Umatrix*alpha; % WF after the calculated gate implementation
    r_ket = Y*alpha; % WF after the ideal gate implementation
    F = F + (abs(ctranspose(r_ket)*ket))^2;
end
F = F/6;
end
