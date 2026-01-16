function [q,qH] = modalbasis_bal_res(alpha,beta,U,Re,nmodesy,x,z)
% Function that computes balanced modes resolvent analysis

% Grid parameters
N           = length(U);
Nx          = length(x);
Nz          = length(z);
[~,A]       = chebdif(N,2);
D           = squeeze(A(:,:,1));
D2          = squeeze(A(:,:,2));
np          = N;

ZZ          = zeros(size(D));
[~,W]       = clencurt(N-1);
W           = diag(W);
Wfull       = [W ZZ ZZ;
               ZZ W ZZ;
               ZZ ZZ W];

% Observation and forcing restriction matrices
ZZ          = zeros(size(D)); II = eye(size(D));
IIe         = (II+fliplr(II))/2; % consider even forcing
IIo         = (II-fliplr(II))/2; % consider odd forcing

nm          = 50; % Number of resolvent modes per frequency

% Initialize response and forcing modes
Ufield      = zeros(3*np,Nx,Nz,nm);
Vfield      = zeros(3*np,Nx,Nz,nm);

% Define frequency vector 
omega       = linspace(-2*abs(alpha)-1,2*abs(alpha)+1,21);
% "Snapshot" matrices
X           = zeros(3*np*Nx*Nz,2*nm*length(omega));
Y           = zeros(3*np*Nx*Nz,2*nm*length(omega));

% Resolvent analysis
for k=1:length(omega)

    C       = -1i*omega(k)*II + 1i*alpha*diag(U);
    V       = -1/Re*(-alpha^2*II + D2 -beta^2*II);

    L       = [1i*alpha*II D 1i*beta*II ZZ;
               C+V diag(D*U) ZZ 1i*alpha*II;
               ZZ C+V ZZ D;
               ZZ ZZ C+V 1i*beta*II];

    H       = [II ZZ ZZ ZZ;
               ZZ II ZZ ZZ;
               ZZ ZZ II ZZ];

    B       = [ZZ ZZ ZZ;
               IIe ZZ ZZ;
               ZZ IIo ZZ;
               ZZ ZZ IIe];

    % Boundary conditions
    % u = 0
    L(np+1,:)     = 0; L(np+1,1) = 1; B(np+1,:) = 0;
    L(2*np,:)     = 0; L(2*np,np) = 1; B(2*np,:) = 0;
    % v = 0
    L(2*np+1,:)   = 0; L(2*np+1,np+1) = 1; B(2*np+1,:) = 0;
    if alpha==0 && beta==0
        %p=0
        L(1*np,:) = 0; L(1*np,4*np) = 1; B(1*np,:) = 0;
    else
        %v=0
        L(3*np,:) = 0; L(3*np,2*np) = 1; B(3*np,:) = 0;
    end
    % w = 0
    L(3*np+1,:)   = 0; L(3*np+1,2*np+1) = 1; B(3*np+1,:) = 0;
    L(4*np,:)     = 0; L(4*np,3*np) = 1; B(4*np,:) = 0;

    % Cholesky factorization
    Wsq            = sqrtm(Wfull);

    % Modified resolvent operator
    R              = H*inv(L)*B;
    Rtil           = (Wsq)*R*inv(Wsq);
    % Singular-value decomposition
    [Util,S,Vtil]  = svdecon(Rtil);

    % Response and forcing modes
    Up             = inv(Wsq)*Util;
    Vp             = inv(Wsq)*Vtil;

    % Singular values
    lambda         = diag(S(1:nm,1:nm));


    % Building global modes
    for ix=1:length(x)
        for iz=1:length(z)
            Ufield(:,ix,iz,:)   = Up(:,1:nm).*exp(1i*alpha*x(ix))*exp(1i*beta*z(iz));
            Vfield(:,ix,iz,:)   = Vp(:,1:nm).*exp(1i*alpha*x(ix))*exp(1i*beta*z(iz));
        end
    end

    Ufield          = reshape(Ufield,[3*np*Nx*Nz,nm]);
    Ufield_re       = real(Ufield.*(lambda).');
    Ufield_imag     = imag(Ufield.*(lambda).');
    Vfield          = reshape(Vfield,[3*np*Nx*Nz,nm]);
    Vfield_re       = real(Vfield.*(lambda).');
    Vfield_imag     = imag(Vfield.*(lambda).');

    % Fill snapshot matrices
    id_re           = 2*(k-1)*nm+1:2:k*(2*nm)-1;
    id_imag         = 2*(k-1)*nm+2:2:k*2*nm;
    X(:,id_re)      = Ufield_re;
    X(:,id_imag)    = Ufield_imag;
    Y(:,id_re)      = Vfield_re;
    Y(:,id_imag)    = Vfield_imag;
    clear Ufield Vfield
end

% Integration weights
[~,W]               = clencurt(N-1);
Wcol                = repmat(W.',Nx*Nz,1);
Wf                  = 0.5/(Nx*Nz)*spdiags(Wcol,0,length(Wcol),length(Wcol));

ZZ                  = zeros(size(Wf));

Wff                 = [Wf, ZZ, ZZ;
                       ZZ, Wf, ZZ;
                       ZZ, ZZ, Wf];
    
% Balanced modes
[Ubal,Sbal,Vbal]    = svdecon(Y'*Wff*X);
[~,ord]             = sort(diag(Sbal),'descend');
Vcnt                = Vbal(:,ord);Ubal = Ubal(:,ord);Sbal = Sbal(ord,ord);
iSbal               = diag(1./(sqrt(diag(Sbal))));
phi_dir             = X*Vbal*iSbal;
phi_adj             = Y*Ubal*iSbal;

q                   = phi_dir(:,1:nmodesy);
qH                  = phi_adj(:,1:nmodesy);
end
