%% Boeing 747: Robust Regret
% This example is from reference 1. It was used in references 2-4 below.  
% References 2-3 considered (additive) regret with a nominal model.  
% Reference 4 considered (multiplicative) competitive ratio also with the 
% nominal model.  This file uses the method in the robust regret paper
% to the robust regret
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
Nd = size(Bd,2);

Q=eye(Nx);
R=eye(Nu);

%% Formulate general plant
% The general plant typically used in H2/Hinf synthesis is:
%  [e; y] = P [d;u]
% In this particular problem e and y are:
%  e = [sqrt(Q) x; sqrt(R) u]
%  y = [x;d] for casual controllers
Ce = [sqrtm(Q); zeros(Nu,Nx)];
Ne = size(Ce,1);

Ded = [zeros(Nx,Nx); zeros(Nu,Nx)];
Deu = [zeros(Nx,Nu); sqrtm(R)];
De = [Ded Deu];

Cy = [eye(Nx); zeros(Nx)];
Dy = [zeros(Nx,Nx+Nu); eye(Nx) zeros(Nx,Nu)];
Ny = 2*Nx;

Pic = ss(A,[Bd Bu],[Ce;Cy], [De;Dy],1);

Pfi = ss(A,[Bd Bu],Ce, De,1);

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

%% (Nominal) Optimal Competitive Ratio  and (Additive) Regret Controllers
% This section of the code solves for the optimal non-causal, competitive 
% ratio, and standard (additive) regret controllers.  The competitive
% ratio controller matches the results in the '23 TAC paper "Competitive 
% Control" by Goel and Hassibi.

% Bisection tolerances
AbsTol = 1e-3;
RelTol = 1e-4;
Tol = [AbsTol, RelTol];

% Optimal non-causal controller
[Knc,CLnc] = ncsyn(Pfi,Nu);

% Compute optimal competitive ratio
gamd  = 0;
gamJInt = [0, 100];
[Kcr,CLcr,GAM] = regretsyn(Pic,Ny,Nu,{gamd,gamJInt},Tol);
gamC = GAM(2);
fprintf('\nCompetitive control gamma=%4.3f \t gamma^2 = %4.3f',...
    gamC,gamC^2);

% Compute optimal (additive) regret controller
gamdInt = [0, 100];
gamJ  = 1;
[Kar,CLar,GAM] = regretsyn(Pic,Ny,Nu,{gamdInt,gamJ},Tol);
gamR = GAM(1);
fprintf('\nAdditive regret control gamma=%4.3f \t gamma^2 = %4.3f\n',...
    gamR,gamR^2);

%% Nominal Pareto Optimal Regret Controllers
% This section of the code computes the Pareto optimal front of nominal
% regret controllers (gamd,gamJ).

N = 20; % N=10;
gamdData = GAMinf*linspace(0.001,0.999,N);
gamJData = zeros(1,N);

tic
gamJInt = [0, 1.1*gamC];
for i=1:N
    [~,~,GAM] = regretsyn(Pic,Ny,Nu,{gamdData(i),gamJInt},Tol);
    gamJData(i) = GAM(2);
end
tNom = toc


%% Robust Pareto Optimal Regret Controllers
% This section of the code computes the Pareto optimal front of robust
% regret controllers (gamd,gamJ).  Delta below satisfies ||Delta||inf<=1.
% The plant input on the control channel is multiplied by (I+beta*Delta)
% where beta scales the magnitude of the uncertainty.


% Add uncertainty at the plant input on the control channel.
Delta = ultidyn('Delta',[Nu Nu]);
if true
    % Constant uncertainty bound
    beta = 0.6;                 % Uncertainty magnitude
    Punc = Pic*blkdiag( eye(Nd), eye(Nu)+beta*Delta );
else
    % Frequency-dependent uncertainty bound
    Wdel = makeweight(0.1,[2 1],1.1,1);
    sigma( ss(A,Bu,eye(4),0,1),'b', Wdel,'r--' )
    Punc = Pic*blkdiag( eye(Nd), eye(Nu)+Delta*Wdel );
end
robgamJData = zeros(1,N);

tic
robgamJInt = [0, 10*gamC];
for i=1:N
    i
    [~,~,GAM] = robregretsyn(Punc,Ny,Nu,{gamdData(i),robgamJInt},Tol);
    robgamJData(i) = GAM(2);
end
tRob = toc


%% Plot results
figure(1)
plot(gamdData,gamJData,'b-.x'); 
hold on;
plot(gamdData,robgamJData,'r:+'); 
ph = plot(0,gamC,'gs',gamR,1,'g^',GAMinf,0,'go');
hold off
xlabel('\gamma_d');
ylabel('\gamma_J');
grid on;
legend('Nom Reg','Rob Reg','Comp. Ratio','Add. Reg','Hinf');
if exist('garyfyFigure','file'), garyfyFigure, end

% Save figure
exportgraphics(gcf, 'Boeing747RobReg.pdf');

figure(2)
plot(gamJData,gamdData,'b-.x'); 
hold on;
plot(robgamJData,gamdData,'r:+'); 
ph = plot(gamC,0,'gs',1,gamR,'g^',0,GAMinf,'go');
hold off
xlabel('\gamma_J');
ylabel('\gamma_d');
grid on;
legend('Nom Reg','Rob Reg','Comp. Ratio','Add. Reg','Hinf');
if exist('garyfyFigure','file'), garyfyFigure, end
