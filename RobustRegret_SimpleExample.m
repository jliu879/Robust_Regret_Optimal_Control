%% Robust Regret For a Simple System
% This example studies different nominal and robust regret controllers
% on a simple SISO system.

%% Bisection tolerances
AbsTol = 1e-4;
RelTol = 1e-4;
Tol = [AbsTol, RelTol];
opt = hinfsynOptions('RelTol',RelTol,'AbsTol',AbsTol);

%% Sample Time
% This example will creates continuous-time models for the plant and
% weights. These are disretized. The sample time is chosen
% sufficient fast to obtain a good approximation via discretization.
Ts = 1e-3;    % Sample time, sec

%% Uncertain Plant Model

% Nominal Dynamics
G = tf(16,[1 5.6]);
G.InputName = 'v';
G.OutputName = 'y';
G = c2d(G,Ts);

% Uncertain Actuator Model
ActNom = tf(20,[1 20]);
Wunc = makeweight(0.20,8,3);

unc = ultidyn('unc',[1 1],'SampleStateDim',5);

ActNom = c2d(ActNom,Ts);
Wunc = c2d(Wunc,Ts);
unc.Ts = Ts;

Act = ActNom*(1 + Wunc*unc); 
Act.InputName = 'u';
Act.OutputName = 'v';

%% Design Setup

% Performance weight on tracking error
e = 0.01;
M = 2;
wS = 8;
We = makeweight(1/e,wS,1/M);
We.InputName = 'e'; We.OutputName = 'etil';
We = c2d(We,Ts);

% Performance weight on control effort
Wu = makeweight(0.1,wS,1000);
Wu.InputName = 'u'; Wu.OutputName = 'util';
Wu = c2d(Wu,Ts);

% Reference weight: Model that reference is a low frequency signal.
Wd = tf(wS,[1 wS]);  
Wd.InputName = 'dtil'; Wd.OutputName = 'r';
Wd = c2d(Wd,Ts);

% Design interconnection: 
% Punc: Output Feedback with Uncertainty
% Pic: Output Feedback without Uncertainty (i.e. nominal model).
emeas  = sumblk('e = r-y');
ICinputs = {'dtil';'u'};
ICoutputs = {'etil';'util';'e'};
Punc = connect(G,Act,Wu,We,Wd,emeas,ICinputs,ICoutputs);
Pic = minreal(Punc.Nominal,[],false);
Ny = 1;
Nu = 1;
Ne = 2;  
Nd = 1;

% Design interconnection: Full Information
Pfi = Pic(1:Ne,:);

%% Solve (standard) H2/Hinf problems

% H2 synthesis
[K2,CL2,GAM2,INFO2] = h2syn(Pic,Ny,Nu);
fprintf('\nOptimal H2:  ||CL||_2 = %4.3f \t ||CL||_inf = %4.3f',...
    norm(CL2,2), norm(CL2,inf) )

% Hinf synthesis
opt.RelTol = 1e-5;
[Kinf,CLinf,GAMinf,INFOinf] = hinfsyn(Pic,Ny,Nu,opt);
fprintf('\nOptimal Hinf:  ||CL||_2 = %4.3f \t ||CL||_inf = %4.3f',...
    norm(CLinf,2), norm(CLinf,inf) )

Linf = G*Act.Nominal*Kinf;
Sinf = feedback(1,Linf);
Tinf = feedback(Linf,1);

%% (Nominal) Optimal Competitive Ratio  and (Additive) Regret Controllers
% This section of the code solves for the optimal non-causal, competitive 
% ratio, and standard (additive) regret controllers.  The competitive
% ratio controller matches the results in the '23 TAC paper "Competitive 
% Control" by Goel and Hassibi.

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

N = 20; 
gamdData = GAMinf*linspace(0.001,0.999,N);
gamJData = zeros(1,N);

tic
gamJInt = [0, 1.1*gamC];
Knr = cell(N,1);
for i=1:N
    [Knr{i},~,GAM] = regretsyn(Pic,Ny,Nu,{gamdData(i),gamJInt},Tol);
    gamJData(i) = GAM(2);
end
tNom = toc

%% Robust Pareto Optimal Regret Controllers
% This section of the code computes the Pareto optimal front of robust
% regret controllers (gamd,gamJ). 

Nrob=N;
robgamdData = GAMinf*linspace(0.001,1.5,Nrob);
robgamJData = zeros(1,Nrob);

tic
robgamJInt = [0, 10*gamC];
Krr = cell(Nrob,1);
for i=1:Nrob
    i
    [Krr{i},~,GAM] = robregretsyn(Punc,Ny,Nu,{robgamdData(i),robgamJInt},Tol);
    robgamJData(i) = GAM(2);
end
tRob = toc

%% Plot results

% Pareto Front
figure(1)
plot(gamdData,gamJData,'b-.x'); 
hold on;
plot(robgamdData,robgamJData,'r:+'); 
ph = plot(0,gamC,'gs',gamR,1,'g^',GAMinf,0,'go');
hold off
xlabel('\gamma_d');
ylabel('\gamma_J');
grid on;
legend('Nom Reg','Rob Reg','Comp. Ratio  (C)','Add. Reg (A)','Hinf');
if exist('garyfyFigure','file'), garyfyFigure, end

% Closed-loop responses
% Note: The closed-loop CL is 2-by-1. There is only one non-zero singular
% value at each frequency and it is given by SV(w) = sqrt(CL(w)' CL(w)).
Nw = 1e3;
w = logspace(-2,log10(pi/Ts), Nw);

SVnc = sigma(CLnc,w);   SVnc = SVnc(:);
% SVnc = freqresp(CLnc'*CLnc,w); SVnc = sqrt(real(SVnc(:)));
SVinf = sigma(CLinf,w);   SVinf = SVinf(:);
CVcr = sigma(CLcr,w);   CVcr = CVcr(:);

% Plot max singular value of CLnc vs. frequency
figure(2)
semilogx(w, 20*log10(SVnc),'k');
xlabel('Frequency, rad/sec')
ylabel('Magnitude, dB')
xlim( w([1 end]) )
set(gca,'Xtick',[0.01, 0.1, 1, 10, 100 1000]);
grid on;
if exist('garyfyFigure','file'), garyfyFigure, end

% Nominal controllers along the Pareto front
figure(3)
for i=1:N
    m = bode(Knr{i},w);   
    if i==1
        ph3(i) = semilogx(w, 20*log10( m(:) ) , 'b'); % CR
    elseif i==N
        ph3(i) = semilogx(w, 20*log10( m(:) ) , 'r-.'); % Hinf
    else
        ph3(i) = semilogx(w, 20*log10( m(:) ) , 'k:');
    end
    hold on;
end
hold off
xlabel('Frequency, rad/sec')
ylabel('Magnitude, dB')
axis( [w([1 end]) -60 50] )
grid on;
legend(ph3([end 1]),'$K_\infty$','$K_{C}$',...
    'Location','Southwest','Interpreter','Latex');
if exist('garyfyFigure','file'), garyfyFigure, end
set(ph3(2:end-1),'LineWidth',1.5)

% Plot closed-loop Hinf and CR costs along with their bounds
figure(4);
ph4(1) = semilogx(w, 20*log10(SVinf),'r-.');
hold on;
ph4(2) = semilogx(w,20*log10(CVcr),'b');
%ph4(3) = semilogx( w([1 end]), 20*log10(GAMinf*[1 1]),'c:');
%ph4(4) = semilogx( w, 20*log10(gamC*SVnc),'k:');
hold off;
xlabel('Frequency, rad/sec')
ylabel('Magnitude, dB')
xlim( w([1 end]) )
grid on;
legend('$T^\infty$','$T^C$','Location','South','Interpreter','Latex');
%legend('$H_\infty$','CR','$\gamma_\infty$','$\gamma_C\times$ NC',...
%    'Location','South','Interpreter','Latex');
if exist('garyfyFigure','file'), garyfyFigure, end
%set(ph4(3:4),'LineWidth',2)

% Compare robust and nominal CR controllers 
figure(5)
m = bode(Knr{1},w);  % Nominal CR
semilogx(w, 20*log10( m(:) ) , 'b');
hold on;
m = bode(Krr{1},w);  % Robust CR
semilogx(w, 20*log10( m(:) ) , 'r-.'); 
hold off
xlabel('Frequency, rad/sec')
ylabel('Magnitude, dB')
axis( [w([1 end]) -60 50] )
grid on;
legend('$K_{C,Nom}$','$K_{C,Rob}$','Location','Southwest','Interpreter','Latex');
if exist('garyfyFigure','file'), garyfyFigure, end

% Compare samples of closed-loop costs with nominal & robust CR controller
rng(0);   % Set seed to make plot repeatable
Nsamp = 20;
figure(6)
for i=1:Nsamp
    SVnom = sigma( usample(lft(Punc,Knr{1})) ,w);
    subplot(2,1,1)
    ph6nom(i) = semilogx(w, 20*log10(SVnom),'r-.');
    hold on;

    subplot(2,1,2)
    SVrob = sigma( usample(lft(Punc,Krr{1})) ,w);
    ph6rob(i) = semilogx(w, 20*log10(SVrob),'r-.');
    hold on;
end

subplot(2,1,1)
SVnom = sigma( lft(Punc.Nominal,Knr{1}) ,w);
ph6nom(Nsamp+1) = semilogx(w, 20*log10(SVnom),'b');
%xlabel('Frequency, rad/sec')
ylabel('Magnitude, dB')
axis( [w([1 end]) -30 30] )
title('Cost with Nominal CR Controller')
legend(ph6nom([Nsamp+1, 1]),'Nom. Plant','Unc. Plant','Location','NorthEast')
grid on;
hold off;

subplot(2,1,2)
SVrob = sigma( lft(Punc.Nominal,Krr{1}) ,w);
ph6rob(Nsamp+1) = semilogx(w, 20*log10(SVrob),'b');
xlabel('Frequency, rad/sec')
ylabel('Magnitude, dB')
axis( [w([1 end]) -30 30] )
title('Cost with Robust CR Controller')
legend(ph6rob([Nsamp+1, 1]),'Nom. Plant','Unc. Plant','Location','NorthEast')
grid on;
hold off;

if exist('garyfyFigure','file'), garyfyFigure, end
set([ph6nom; ph6rob],'LineWidth',1.5)

% Robust stability margins
fprintf('\nRobust stability margin for nominal CR\n');
SMcr_nom = robstab( lft(Punc,Knr{1}) )
allmargin(G*ActNom*Knr{1})

fprintf('\nRobust stability margin for robust CR\n');
SMcr_rob = robstab( lft(Punc,Krr{1}) )
allmargin(G*ActNom*Krr{1})

% Compare robust and nominal Hinf controllers 
% Note that the robust Hinf controller is a variant of mu-synthesis
% where the worst-case gain is minimized.
figure(7)
m = bode(Knr{end},w);  % Nominal CR
semilogx(w, 20*log10( m(:) ) , 'b');
hold on;
m = bode(Krr{end},w);  % Robust CR
semilogx(w, 20*log10( m(:) ) , 'r-.'); 
hold off
xlabel('Frequency, rad/sec')
ylabel('Magnitude, dB')
axis( [w([1 end]) -60 50] )
grid on;
legend('$K_{\infty,Nom}$','$K_{\infty,Rob}$','Location','Southwest','Interpreter','Latex');
if exist('garyfyFigure','file'), garyfyFigure, end


%% Save Plots
figure(1); exportgraphics(gcf, 'SimpleEx_Pareto.pdf');
figure(2); exportgraphics(gcf, 'SimpleEx_NoncausalCost.pdf');
figure(3); exportgraphics(gcf, 'SimpleEx_NomControllers.pdf');
figure(4); exportgraphics(gcf, 'SimpleEx_NomCosts.pdf');
figure(5); exportgraphics(gcf, 'SimpleEx_CRControllers.pdf');
figure(6); exportgraphics(gcf, 'SimpleEx_CRCosts.pdf');
