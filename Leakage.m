function L = ...
    Leakage(w01, w12, w, T, tstep, SignalString, Theta)

N = length(SignalString);
h = 1.054e-34; % Planck constant
C1 = 1e-12; % Qubit capacitance
F0 = 2.06e-15; % Magnetic flux quantum
Id = [1 0 0; 0 1 0; 0 0 1]; % Identity matrix
V = F0/w; % Voltage
Cc = Theta/(F0*sqrt(2*w01/(h*C1))); % Connection capacitance
Amp = Cc*V*sqrt(h*w01/(2*C1)); % Pulse amplitude
    
% Applying signal to a qubit
% Undisturbed hamiltonian
H0 = [0 0 0; 0 h*w01 0; 0 0 h*w01 + h*w12];
[EigVec, ~] = eig(H0); % Finding eigenvectors and eigenvalues
WF0 = EigVec(:,1); % State |0>
WF1 = EigVec(:,2); % State |1>
WF2 = EigVec(:,3); % State |2>

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

% Different i.c. leakage arrays
Leak1 = zeros(1,N); Leak2 = zeros(1,N); Leak3 = zeros(1,N); 
Leak4 = zeros(1,N); Leak5 = zeros(1,N); Leak6 = zeros(1,N);

% Hamiltonians
% Disturbed hamiltonian
Hmatrix = [0 -1 0; 1 0 -sqrt(2); 0 sqrt(2) 0]; 
HrPlus = 1i*Amp*Hmatrix + H0;
HrMinus = -1i*Amp*Hmatrix + H0;
HrPlus2 = 1i*2*Amp*Hmatrix + H0;
HrMinus2 = -1i*2*Amp*Hmatrix + H0;
HrZero = H0;

% Evolution operators
UPlusStep = UMatrix(HrPlus,tstep);
UMinusStep = UMatrix(HrMinus,tstep);
UPlus2Step = UMatrix(HrPlus2,tstep);
UMinus2Step = UMatrix(HrMinus2,tstep);
UZeroStep = UMatrix(HrZero,tstep);
UPlus = Id; UMinus = Id; UZero = Id; UPlus2 = Id; UMinus2 = Id;

for i = 1:1:(w/tstep)
    UPlus = UPlusStep*UPlus;
    UMinus = UMinusStep*UMinus;
    UPlus2 = UPlus2Step*UPlus2;
    UMinus2 = UMinus2Step*UMinus2;
    UZero = UZeroStep*UZero;
end

% Evolution operator for the whole period
UTStep = UMatrix(HrZero,tstep);
UT = Id;

for i = 1:1:((T - w)/tstep)
    UT = UTStep*UT;
end
UTPlus = UT*UPlus;
UTMinus = UT*UMinus;
UTPlus2 = UT*UPlus2;
UTMinus2 = UT*UMinus2;
UTZero = UT*UZero;

for IS = 1:1:6
    Prob1 = zeros(1,N); Prob2 = zeros(1,N); Prob3 = zeros(1,N);
    WF = InitStates(:,IS); % WF i.c.
    for j = 1:1:N
        switch SignalString(j)
            case 2
                U = UTPlus2;
            case -2
                U = UTMinus2;
            case 1
                U = UTPlus;
            case -1
                U = UTMinus;
            case 0
                U = UTZero;
        end
        WF = U*WF; % Applying evolution operator
        % Calculating state probabilities
        Prob1(j) = abs(ctranspose(WF0)*WF)^2;
        Prob2(j) = abs(ctranspose(WF1)*WF)^2;
        Prob3(j) = abs(ctranspose(WF2)*WF)^2;
    end   
    switch IS
        case 1
            Leak1 = Prob3;
        case 2
            Leak2 = Prob3;
        case 3
            Leak3 = Prob3;
        case 4
            Leak4 = Prob3;
        case 5
            Leak5 = Prob3;
        case 6
            Leak6 = Prob3;
    end
end
L = 1/6*(Leak1(N) + Leak2(N) + Leak3(N) + Leak4(N) + Leak5(N) + Leak6(N));
end
