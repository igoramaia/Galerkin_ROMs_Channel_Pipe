function jpat = jacob_pat(~,X,L,Q,Re);
J = L/Re;

% evaluate Jacobian
for j=1:size(L,1)
    J = J + squeeze(Q(:,j,:))*X(j)+squeeze(Q(:,:,j))*X(j);
end
J= sparse(J);
[jr, jc] = size(J);
jpat = sparse(jr, jc);
jpat(find(J~=0)) = 1;
