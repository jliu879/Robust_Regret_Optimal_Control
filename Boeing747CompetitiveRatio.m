%% Boeing 747: Competitve Ratio
% This example is from reference 1. It was used in references 2-4 below.  
% References 2-3 considered (additive) regret with a nominal model.  
% Reference 4 considered (multiplicative) competitive ratio also with the 
% nominal model.  This file uses the method in the robust regret paper
% to compute the optimal competitve ratio controller.  The results match
% those given in Reference 4 as expected.
%
% Refs
% [1] J. Hong, N. Moehle, and S. Boyd, lecture notes in "Introduction to 
% Matrix Methods: LQR", https://web.stanford.edu/class/engr108.
% [2] Sabag, Goel, Lale, and Hassibi, Regret-Optimal Full-Information
% Control, ACC, p.4777-4782 2021.
% [3] Sabag, Goel, Lale, and Hassibi, Regret-Optimal Full-Information
% Control, arXiv, 2021.
% [4] Goel and Hassibi, Competitive Control, IEEE TAC, 2022.


%% Plant and Cost Data
% Dynamics: x(t+1) = A x(t) + Bd d(t) + Bu u(t)
% Per-step cost:     x(t)' Q x(t) + u(t)' R u(t) 
A = [0.99 0.03 -0.02 -0.32; 0.01 0.47 4.7 0; ...
    0.02 -0.06 0.4 0; 0.01 -0.04 0.72 0.99];
Bd = eye(4);
Bu = [0.01 0.99; -3.44 1.66; -0.83 0.44; -0.47 0.25];

[Nx,Nu] = size(Bu);

Q=eye(Nx);
R=eye(Nu);

%% Formulate general plant
% The general plant typically used in H2/Hinf synthesis is:
%  [e; y] = P [d;u]
% In this particular problem e and y are:
%  e = [sqrt(Q) x; sqrt(R) u]
%  y = [x;w] for casual controllers
Ce = [sqrtm(Q); zeros(Nu,Nx)];

Ded = [zeros(Nx,Nx); zeros(Nu,Nx)];
Deu = [zeros(Nx,Nu); sqrtm(R)];
De = [Ded Deu];

Cy = [eye(Nx); zeros(Nx)];
Dy = [zeros(Nx,Nx+Nu); eye(Nx) zeros(Nx,Nu)];
Ny = 2*Nx;

Pic = ss(A,[Bd Bu],[Ce;Cy], [De;Dy],1);

%% Solve H2/Hinf problems

% K2 is a dynamic controller. The code solves for the optimal
% state-feedback/observer assuming this is a general problem even if we
% set y=x or y=[x,w]. K2 will be "almost" a static controller if y=x or
% y=[x,w].  The code also returns the exact static gain for the full
%information case y=[x,w] in the field INFO2.KFI. This full-information
% gain and is equal to [-Kx -Kw] where Kx and Kw are computed below.
[K2,CL2,GAM2,INFO2] = h2syn(Pic,Ny,Nu);
fprintf('\nOptimal H2:  ||CL||_2 = %4.3f \t ||CL||_inf = %4.3f',...
    norm(CL2,2), norm(CL2,inf) )

% Kinf is also a dynamic controller. We can use HINFFI to solve the
% full information problem corresponding to y=[x;w].
[Kinf,CLinf,GAMinf,INFOinf] = hinfsyn(Pic,Ny,Nu);
fprintf('\nOptimal Hinf:  ||CL||_2 = %4.3f \t ||CL||_inf = %4.3f',...
    norm(CLinf,2), norm(CLinf,inf) )

% HINFI solves the full information Hinf problem with y=[x;w].  This
% returns the static full information gain Kfi.
Pfi = ss(A,[Bd Bu],Ce, De,1);
opt = hinfsynOptions('RelTol',1e-4);
[Kfi,CLfi,GAMfi] = hinffi(Pfi,Nu,opt);
fprintf('\nOptimal Hinf FI:  ||CL||_2 = %4.3f \t ||CL||_inf = %4.3f\n',...
    norm(CLfi,2), norm(CLfi,inf) )
CLinf = CLfi;

%% (Nominal) Optimal Competitive Ratio Controller
% This section of the code solves for the optimal competitive ratio 
% control.  This matches the results in the '23 TAC paper "Competitive 
% Control" by Goel and Hassibi.

% Optimal non-causal controller
[Knc,CLnc] = ncsyn(Pfi,Nu);

% Compute optimal competitive ratio
gamCInt = [0, 100];
gamd  = 0;
[KC,CLcr,GAM] = regretsyn(Pic,Ny,Nu,{gamd,gamCInt},[1e-3, 1e-4]);
gamC = GAM(2);
fprintf('\nCompetitive control gamma=%4.3f \t gamma^2 = %4.3f\n',...
    gamC,gamC^2);

%% Plot Results

% Plot: Maximum singular value vs. frequency
% This should match Fig 1 in the '23 TAC paper by Goel and Hassibi
Nfreq = 1e3;
w=linspace(0,2*pi,Nfreq);
SV2 = sigma(CL2,w);
SVinf = sigma(CLinf,w);
SVcr = sigma(CLcr,w);
SVnc = sqrt(sigma(CLnc'*CLnc,w));

figure(1)
plot(w,SV2(1,:),'g',w,SVinf(1,:),'r',w,SVcr(1,:),'b--',w,SVnc(1,:),'k');
xlim([0, 2*pi]);
xlabel('Freq, rad/sec');
ylabel('Max SV');
grid on;
legend('H2','Hinf','Comp.','NC')
if exist('garyfyFigure','file'), garyfyFigure, end

% Plot: Competitive Ratio vs. frequency
% This should match Fig 2 in the '23 TAC paper by Goel and Hassibi
PI = CLnc'*CLnc;
CRnc = ones(Nfreq,1);


CRinf = zeros(Nfreq,1);
tmp=freqresp(CLinf'*CLinf/PI, w);
for i=1:Nfreq
    CRinf(i) = max(abs(eig( tmp(:,:,i) )));
end

CR2 = zeros(Nfreq,1);
tmp=freqresp(CL2'*CL2/PI, w);
for i=1:Nfreq
    CR2(i) = max(abs(eig( tmp(:,:,i) )));
end

CRcr = zeros(Nfreq,1);
tmp=freqresp(CLcr'*CLcr/PI, w);
for i=1:Nfreq
    CRcr(i) = max(abs(eig( tmp(:,:,i) )));
end

figure(2);
plot(w,CR2,'g',w,CRinf,'r',w,CRcr,'b--',w,CRnc,'k');
xlim([0, 2*pi]);
xlabel('Freq, rad/sec');
ylabel('Max SV');
grid on;
legend('H2','Hinf','Comp.','NC')
if exist('garyfyFigure','file'), garyfyFigure, end
