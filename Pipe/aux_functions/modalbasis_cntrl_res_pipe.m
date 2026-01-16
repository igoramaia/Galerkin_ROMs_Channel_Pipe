function U_control = modalbasis_cntrl_res_pipe(U, r, x, theta, Re, alfa, omega, m, D, D2, W_half, Wg, Nr, Nx, Nm, modes_per_omega, Nmodos)
% Compute controllability modes for pipe flow

% Initialising matrices
nomega             = length(omega);
X                  = zeros(3*Nr*Nx*Nm, 2*modes_per_omega*(nomega - 1));
q_field            = zeros(3*Nr, Nx, Nm, modes_per_omega); 
col                = 1; % counter

% Compute comtrollability Gramian
for i = 1:nomega
    % Análise resolvente
    [U_tilde, S, ~]     = resolvent(U, r, Re, alfa, omega(i), m, D, D2, W_half);
    gains               = diag(S);
    
    % Most amplified response mdoes
    q                   = W_half\U_tilde(:, 1:modes_per_omega);
   
    % Building global modes
    for j = 1:Nx
        for k = 1:Nm
            % As 3 primeiras dimensões dessa matriz correspondem às coordenadas (r, x, θ) do campo de velocidades, e a última corresponde ao modo
            q_field(:, j, k, :)     = q*exp(1i*(alfa*x(j) + m*theta(k)));
        end
    end

    q_field_vec                     = reshape(q_field, [3*Nr*Nx*Nm, modes_per_omega]);
    
    % Build snapshot matrix
    for j = 1:modes_per_omega
        % Parte real do j-ésimo modo mais amplificado multiplicado pelo seu respectivo ganho
        X(:, col)                   = gains(j)*real(q_field_vec(:, j)); col = col + 1;
        % Se omega = 0, a parte imaginária é nula
        if omega ~= 0
            % Parte imaginária do j-ésimo modo mais amplificado multiplicado pelo seu respectivo ganho
            X(:, col)               = gains(j)*imag(q_field_vec(:, j)); col = col + 1;
        end
    end
end

% SVD of the Gramiano
G               = X'*Wg*X; % Gramiano
[~, Sigma, Psi] = svd(G);
% Modos de controlabilidade
U_control       = X*Psi*Sigma^(-1/2);
U_control       = U_control(:, 1:Nmodos);
end