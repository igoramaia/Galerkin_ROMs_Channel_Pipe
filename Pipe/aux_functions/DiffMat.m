function [D_1, D_2, D_3] = DiffMat(N, x, M)
if M == N
    D_1     = zeros(N,N); D_2=zeros(N,N); D_3=zeros(N,N);
else
    D_1     = sparse(N,N); D_2=sparse(N,N); D_3=sparse(N,N);
end
%
for i=1:N
    if M==N
        xb  = x; i_s= i;
        [D_1(i,1:M), D_2(i,1:M), D_3(i,1:M)]    = Dbs(i_s,xb,M);
    else
        M_c=(M+1)/2;
        if i < M_c
            xb =x(1:M); i_s=i;
            [D_1(i,1:M), D_2(i,1:M), D_3(i,1:M)] = Dbs(i_s,xb,M);
        end
        if i >= M_c && i<= N-M_c
            xb = x((i-M_c+1):(i+M_c-1)); i_s = M_c;
            [D_1(i,(i-M_c+1):(i+M_c-1)), D_2(i,(i-M_c+1):(i+M_c-1)), ...
                D_3(i,(i-M_c+1):(i+M_c-1))] = Dbs(i_s,xb,M);
        end
        if i > N-M_c
            xb =x((N-M+1):N); i_s = i+M-N;
            [D_1(i,(N-M+1):N), D_2(i,(N-M+1):N), D_3(i,(N-M+1):N)] = Dbs(i_s,xb,M);
        end
    end
end
%
%
%
function[db_1, db_2, db_3]=Dbs(i,xb,M)
% Differentiation matrices:
Db_1=zeros(M,M); Db_2=zeros(M,M); Db_3=zeros(M,M);
%Compute factors a_i/a_m and store values for each m not equal to i
aiam=zeros(1,M);
F_i=xb(i)-xb; F_i(i)=1; n_i=length(find(F_i<0));
for m=1:M
    if abs(m-i)>0
        %    else
        F_m=xb(m)-xb; F_m(m)=1; n_m=length(find(F_m <0));
        aiam(m)=(-1)^(n_i-n_m)*prod(sort(abs(F_i))./sort(abs(F_m)));
    end
end
%
% s=1:
for m=1:M
    if abs(m-i)>0
        Db_1(i,m)=aiam(m)/(xb(i)-xb(m));
    end
end
Db_1(i,i)=-sum(Db_1(i,:));
% s=2
for m=1:M
    if abs(m-i)>0
        Db_2(i,m)=2/(xb(i)-xb(m))*(aiam(m)*Db_1(i,i)-Db_1(i,m));
    end
end
Db_2(i,i)=-sum(Db_2(i,:));
% s=3
for m=1:M
    if abs(m-i)>0
        Db_3(i,m)=3/(xb(i)-xb(m))*(aiam(m)*Db_2(i,i)-Db_2(i,m));
    end
end
Db_3(i,i)=-sum(Db_3(i,:));
%
% Output
db_1=Db_1(i,:); db_2=Db_2(i,:); db_3=Db_3(i,:);

