function jac = jacob(~,X,L,Q,Re);
J = L/Re;
% evaluate Jacobian
for j=1:size(L,1)
    J = J + squeeze(Q(:,j,:))*X(j)+squeeze(Q(:,:,j))*X(j);
end

jac = sparse(J);

