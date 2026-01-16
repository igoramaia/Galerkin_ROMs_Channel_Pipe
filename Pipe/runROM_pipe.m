%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Galerkin system of ODE's and reconstruct velocity field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
load('ROM_Galerkin_Cilindro_4_alfas_7_ms_alfa_max_3_m_max_3_Ncontrol_48_Nr_64_Nx_14_Nm_14_Nmodos_1200.mat');
toc

nmodes               = size(phi,5);
[x,Dx2]              = fourdif(Nx,2);
[x,Dx]               = fourdif(Nx,1);
[z,Dz2]              = fourdif(Nm,2);
[z,Dz]               = fourdif(Nm,1);

Lx = max(x)+x(2)-x(1);
Lz = max(theta)+theta(2)-theta(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set operators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


facx        = Lx/(2*pi); x = x*facx; Dx = Dx/facx; Dx2 = Dx2/facx^2;
facz        = Lz/(2*pi); z = z*facz; Dz = Dz/facz; Dz2 = Dz2/facz^2;
dx          = x(2)-x(1);
dz          = z(2)-z(1);

% Initialise pressure gradient term

F          = zeros(Nmodos,1);
bf         = zeros(Nr,Nx,Nm,3);
[y,A]      = chebdif(Nr,2);
[~, W]     = clencurt(2*Nr - 1);
W          = W(1:Nr).*r';
Wcol       = repmat(W.',Nx*Nm,1);
W          = diag(W);
bf(:,:,:,1)= 2; % Pressure gradient value (after normalisation)

% Compute operator corresponding to pressure gradient
for i=1:Nmodos
    if exist('phiH','var')
        F(i) = inprod(phiH(:,:,:,:,i),bf,dx,dz,W);
    else
        F(i) = inprod(Phi(:,:,:,:,i),bf,dx,dz,W); 
    end
end

F           = (abs(F)>1e-6).*F;

Q_aux       = full(Q);
Qn          = reshape(Q_aux,[Nmodos Nmodos Nmodos]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run ROM simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Retau       = 180;

% Laminar solution (Q = 0)
qlam = -(L/Retau)\F;

% Time-integration
tspan           = 0:0.1:100;
q0              = 0.1*qlam + 0.5*(rand(Nmodos, 1) - 0.5); % Initial condition
%options        = odeset('RelTol',1e-4,'AbsTol',1e-6,'OutputFcn',@odeplot);
%options        = odeset('RelTol',1e-4);
options          = odeset('RelTol',1e-5,'AbsTol',1e-6,'Jacobian',@(t,X) jacob(t,X,L,Qn,Retau),'JPattern',@(t,X) jacob_pat(t,X,L,Qn,Retau));
disp(['Nmodes = ' int2str(Nmodos)])
disp('Simulation');
tic;
[t, q]          = ode15s(@(t, X) galerkinsys(t, X, L, Q, F, Retau), tspan, q0, options);
toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reconstruct velocity field
% u(x,y,z,t) =\sum q(t)\phi(x,y,z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nt              = length(t);
posend          = floor(nt/4):nt;
posstats        = floor(nt/4):1:nt;
qfield          = zeros(size(Phi,1),size(Phi,2),size(Phi,3),size(Phi,4),length(posstats));
for j=1:length(posstats)
    disp(['Statistics, ' int2str(j) '/' int2str(length(posstats))]);
    for i=1:Nmodos
        qfield(:,:,:,:,j) = qfield(:,:,:,:,j)+squeeze(Phi(:,:,:,:,i))*q(posstats(j),i);
    end
end
qprime          = qfield - mean(qfield,5);
qrms            = squeeze(sqrt(mean(mean(mean(qprime.^2,5),3),2)));

Umean           = mean(qfield,5);
Umean           = mean(mean(Umean,3),2);
Umean           = squeeze(Umean);

% Plot stats

yplus = (1 - r)*Retau;
figure
semilogx(yplus,Umean)
hold on
xlim([1 Retau])

figure
semilogx(yplus,qrms(:,1))
hold on
semilogx(yplus,qrms(:,2))
semilogx(yplus,qrms(:,3))
xlim([1 Retau])
