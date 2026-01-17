%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Galerkin system of ODE's and reconstruct velocity field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%
%
%
clear
close all;
clc

addpath(genpath('aux_functions'));

set(groot,'defaulttextInterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaultColorbarTickLabelInterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex');
set(0,'defaultfigureunits','inches');
set(0,'defaultfigureposition',[7 7 6 4.5])
set(0,'DefaultLineLinewidth',2);
set(0,'DefaultAxesFontSize',18);
set(0,'DefaultFigureColor','remove')
set(0,'DefaultFigureColor',[1 1 1])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Loading model');
tic;
load('galerkin_model_cntrl_2_48_3_Lx_1_Lz_0.5_Re_100.mat');
toc

nmodes               = size(phi,5);
nalphas              = 2; ngammas = 2;
[x,Dx2]              = fourdif(Nx,2);
[x,Dx]               = fourdif(Nx,1);
[z,Dz2]              = fourdif(Nz,2);
[z,Dz]               = fourdif(Nz,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set operators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Retau               = 180;                  % Friction Reynolds number for ROM simulation
L0                  = L;
L                   = L0 + Re*L2;           % Complete linear term (diffusion + energy extraction from base flow)


dx                  = x(2)-x(1);
dz                  = z(2)-z(1);

% Initialise pressure gradient term
F                   = zeros(nmodes,1);
bf                  = zeros(Ny,Nx,Nz,3);
[y,A]               = chebdif(Ny,2);
[~,W]               = clencurt(Ny-1);
Wcol                = repmat(W.',Nx*Nz,1);
W                   = diag(W);
bf(:,:,:,1)         = 1; % Pressure gradient value (after normalisation)

% Compute operator corresponding to pressure gradient
for i=1:nmodes
    if exist('phiH','var') % For balanced modes
        F(i)        = inprod(phiH(:,:,:,:,i),bf,dx,dz,W);
    else                   % For controllability modes
       F(i)         = inprod(phi(:,:,:,:,i),bf,dx,dz,W); 
    end
end
F                   = (abs(F)>1e-6).*F;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run ROM simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Laminar solution
qlam                = -(L/Retau)\F;
Ulam                = zeros(Ny,1);

pos                 = find(kvec(:,1)==0&kvec(:,2)==0);

for i=pos.'
    Ulam            = Ulam+squeeze(phi(:,1,1,1,i))*qlam(i);
end

Ubulklam            = -trapz(y,Ulam)/2;


% Initial condition (laminar solution + fluctuations)
q0 = 0.3*qlam + 2.0*(rand(nmodes,1)-0.5);

tspan               = 0:0.1:200;
options             = odeset('RelTol',1e-6);
disp(['Nmodes = ' int2str(nmodes)])
disp('Simulation');
% Time integration
tic;
[t,q]               = ode45(@(t,X) galerkinsys(t,X,L,QQ,F,Retau),tspan,q0,options);
toc

% Plot mode coefficients
figure
plot(t,q)
xlabel('$t$')
ylabel('$a_i$')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reconstruct velocity field
% u(x,y,z,t) =\sum q(t)\phi(x,y,z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nt                        = length(t); 
posend                    = floor(nt/2):1:nt;
posstats                  = floor(nt/2):1:nt;
qfield                    = zeros(size(phi,1),size(phi,2),size(phi,3),size(phi,4),length(posend));
for j=1:length(posstats)
    disp(['Statistics, ' int2str(j) '/' int2str(length(posstats))]);
    qfield(:,:,:,:,j)     = phi0;
    for i=1:nmodes
        qfield(:,:,:,:,j) = qfield(:,:,:,:,j)+squeeze(phi(:,:,:,:,i))*q(posend(j),i);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute flow statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U0                        = phi0(:,1,1,1);
for i=pos.'
    U0                    = U0+squeeze(phi(:,1,1,1,i))*mean(q(posend,i));
end
Ubulk                     = diag(W).'*U0 / 2; %bulk velocity


qprime                    = qfield - mean(qfield,5);

qrms                      = squeeze(sqrt(mean(mean(mean(qprime.^2,5),3),2)));

% Plot stats
yplus = (y+1)*Retau;
figure
semilogx(yplus,U0)
hold on
xlim([1 Retau])

figure
semilogx(yplus,qrms(:,1))
hold on
semilogx(yplus,qrms(:,2))
semilogx(yplus,qrms(:,3))
xlim([1 Retau])
