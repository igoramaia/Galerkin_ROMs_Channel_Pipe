function M = orthogonality_test(Phi, W)
Nmodes = size(Phi, 5);
M = zeros(Nmodes, Nmodes);
for i = 1:Nmodes
    for j = 1:Nmodes
        mode_i = squeeze(Phi(:, :, :, :, i)); mode_i = mode_i(:);
        mode_j = squeeze(Phi(:, :, :, :, j)); mode_j = mode_j(:);
        M(i, j) = mode_i'*W*mode_j;
    end
end
end