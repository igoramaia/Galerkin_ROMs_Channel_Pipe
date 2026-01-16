function [U_tilde, S, V_tilde] = resolvent(U, r, Re, alfa, omega, m, D, D2, W_half)
% Resolvent analysis of pipe flow  

% Identity and null matrices
I           = eye(size(D));
Z           = zeros(size(D));

% Auxiliary delta matrix
Delta       = -alfa^2*I - (m^2 + 1)*r^(-2) + r^(-1)*D + D2;

% Input-output matrices
% M*(dq/dt) = Aq + Bf / q --> saída,      f --> forçagem
%     z     = Cq + n  / z --> observação, n --> ruído das medidas

B           = [I Z Z
               Z I Z
               Z Z I
               Z Z Z];

C           = [I Z Z Z
               Z I Z Z
               Z Z I Z];

M           = [I Z Z Z
               Z I Z Z
               Z Z I Z
               Z Z Z Z];

% ------------------------------------------------------------------------------------------------------------- %
A = [(-1i*alfa*U + 1/Re*(Delta + r^(-2)))            -D*U                      Z              -1i*alfa*I;       %
    Z                       (-1i*alfa*U + Delta/Re)     -2i*m*r^(-2)/Re           -D;           %
    Z                           2i*m*r^(-2)/Re      (-1i*alfa*U + Delta/Re)  -1i*m*r^(-1);      %
    1i*alfa*I                        (D + r^(-1))             1i*m*r^(-1)              Z        ];  %
% ------------------------------------------------------------------------------------------------------------- %

% Linearised Navier-Stokes operator in cylindrical coordinates

L           = -(A + 1i*omega*M);

% Boundary conditions
N           = size(r, 1);

% ux(r) = 0
L(1, :)     = 0; L(1, 1) = 1;
B(1, :)     = 0;

% ur(r) = 0
L(N+1, :)   = 0; L(N+1, N+1) = 1;
B(N+1, :)   = 0;

% um(r) = 0
L(2*N+1, :) = 0; L(2*N+1, 2*N+1) = 1;
B(2*N+1, :) = 0;

% Additional boundary condition for mean flow 
if alfa == 0 && m == 0
    % p(r) = 0
    L(3*N+1, :) = 0; L(3*N+1, 3*N+1) = 1;
    B(3*N+1, :) = 0;
end

% Resolvent
R           = C/L*B;
R_tilde     = W_half*R/W_half;
% Resolvent modes
try
    [U_tilde, S, V_tilde] = svd(R_tilde);
catch err
    % To avoid conflict of versions
    if err.identifier == "MATLAB:svd:svdNoConvergence"
        [U_tilde, S, V_tilde]   = svd(R_tilde + 1e-14*norm(R_tilde)*randn(RandStream("mt19937ar"), size(R_tilde)));
    else
        rethrow(err)
    end
end
end