function inp = inprod(f1,f2,dx,dz,W)
[Ny,Nx,Nz,~]=size(f1);
inty = zeros(Nx,Nz);
for i=1:Nz
%    for j=1:Nx
        inty(:,i) = diag(squeeze(f1(:,:,i,1))'*W*squeeze(f2(:,:,i,1)) + squeeze(f1(:,:,i,2))'*W*squeeze(f2(:,:,i,2)) + squeeze(f1(:,:,i,3))'*W*squeeze(f2(:,:,i,3)));
        %inty(j,i) = (squeeze(f1(:,j,i,1))'*W*squeeze(f2(:,j,i,1)) + squeeze(f1(:,j,i,2))'*W*squeeze(f2(:,j,i,2)) + squeeze(f1(:,j,i,3))'*W*squeeze(f2(:,j,i,3)));
    %end
end
inp = sum(sum((inty)))*dx*dz/(2*Nx*dx*Nz*dz);
