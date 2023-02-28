function U = UMatrix(H,t)
Id = [1 0 0; 0 1 0; 0 0 1];
h = 1.054e-34;
U = (Id - 1i*H*t/(2*h))/(Id + 1i*H*t/(2*h));