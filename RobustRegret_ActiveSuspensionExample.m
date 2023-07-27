%% Robust Regret Control of Active Suspension
%
% This example is based on Matlab's "Robust Control of Active Suspension"
% demo and the corresponding demo livescript ActiveSuspensionExample.mlx. 
% The code below to formulate the problem is mainly from the Matlab demo 
% livscript (with some simplifications). The livescript has been modified 
% to compare against controllers synthesized using robust regret. See the 
% demo and livescript for details. There is also a Matlab Tech Talk video 
% by Brian Douglas that provides additional context on the control designs.

%% Sample Time
% The original Matlab demo used continuous-time models for the design.
% This example will use discrete-time models. The sample time is chosen
% sufficient fast to obtain a good approximation via discretization.
Ts = 2e-3;    % Sample time, sec

%% Quarter-Car Suspension Model
% Construct a state-space model |qcar|.

% Physical parameters
mb = 300;    % kg
mw = 60;     % kg
bs = 1000;   % N/m/s
ks = 16000 ; % N/m
kt = 190000; % N/m

% State matrices
A = [ 0 1 0 0; [-ks -bs ks bs]/mb ; ...
      0 0 0 1; [ks bs -ks-kt -bs]/mw];
B = [ 0 0; 0 1e3/mb ; 0 0 ; [kt -1e3]/mw];
C = [1 0 0 0; 1 0 -1 0; A(2,:)];
D = [0 0; 0 0; B(2,:)];

qcar = ss(A,B,C,D);
qcar.StateName = {'body travel (m)';'body vel (m/s)';...
          'wheel travel (m)';'wheel vel (m/s)'};
qcar.InputName = {'r';'fs'};
qcar.OutputName = {'xb';'sd';'ab'};

% Discretize model
qcar = c2d(qcar,Ts);

%% Uncertain Actuator Model
% The hydraulic actuator used for active suspension control is connected between 
% the body mass $m_b$ and the wheel assembly mass $m_w$. Construct an
% uncertain model for the actuator dynamics.

ActNom = tf(1,[1/60 1]);
Wunc = makeweight(0.40,15,3);
unc = ultidyn('unc',[1 1],'SampleStateDim',5);

ActNom = c2d(ActNom,Ts);
Wunc = c2d(Wunc,Ts);
unc.Ts = Ts;

Act = ActNom*(1 + Wunc*unc);
Act.InputName = 'u';
Act.OutputName = 'fs';


%% Design Setup

% beta changes weights between "comfort" (beta near 0) and handling (beta
% near 1).  A balanced objective is obtained with beta = 0.5.
beta = 0.5;  

Wroad = ss(0.07);  Wroad.u = 'd1';   Wroad.y = 'r';
Wd2 = ss(0.01);  Wd2.u = 'd2';   Wd2.y = 'Wd2';
Wd3 = ss(0.5);   Wd3.u = 'd3';   Wd3.y = 'Wd3';

Wact = 0.8*tf([1 50],[1 500]);  
Wact.u = 'u';  Wact.y = 'e1';
Wact = c2d(Wact,Ts);

HandlingTarget = 0.04 * tf([1/8 1],[1/80 1]);
ComfortTarget = 0.4 * tf([1/0.45 1],[1/150 1]);
Wsd = beta / HandlingTarget;
Wsd.u = 'sd';  Wsd.y = 'e3';
Wsd = c2d(Wsd,Ts);

Wab = (1-beta) / ComfortTarget;
Wab.u = 'ab';  Wab.y = 'e2';
Wab = c2d(Wab,Ts);

sdmeas  = sumblk('y1 = sd+Wd2');
abmeas = sumblk('y2 = ab+Wd3');
ICinputs = {'d1';'d2';'d3';'u'};
ICoutputs = {'e1';'e2';'e3';'y1';'y2'};
qcaric = connect(qcar(2:3,:),Act,Wroad,Wact,Wab,Wsd,Wd2,Wd3,...
                 sdmeas,abmeas,ICinputs,ICoutputs);

%% Nominal Designs
% Use |hinfsyn| to compute an $H_\infty$ controller. Compare against
% the optimal non-causal controller and a nominal additive regret design.

% H-infinity
ncont = 1; % one control signal, u
nmeas = 2; % two measurement signals, sd and ab
[Kinf,CLinfW,gaminf] = hinfsyn(qcaric,nmeas,ncont);
gaminf

% Optimal non-causal controller
Ne = size(qcaric,1) - nmeas;
Pfi = qcaric(1:Ne,:);
Pfi = minreal(Pfi.Nominal,[],false);
[Knc,CLncW] = ncsyn(Pfi,ncont);

% Additive Regret
AbsTol = 1e-3; RelTol = 1e-4; Tol = [AbsTol, RelTol];
gamdInt = [0, 10];
gamJ  = 1;
[Kar,CLarW,GAM] = regretsyn(qcaric,nmeas,ncont,{gamdInt,gamJ},Tol);
gamar = GAM(1)

% Closed-loop models (without design weights)
% (Note: Then non-causal construction below simply inverts out the
% weights on the road (first input) and body accel (second output). This
% is a non-minimal realization. It is only used as a short-cut to
% plot the closed-loop response.)
Kinf.u = {'sd','ab'};  Kinf.y = 'u';
CLinf = connect(qcar,Act,Kinf,'r',{'xb';'sd';'ab'; 'fs'});
Kar.u = {'sd','ab'};  Kar.y = 'u';
CLar = connect(qcar,Act,Kar,'r',{'xb';'sd';'ab'; 'fs'});

CLncab = (1/Wab)*CLncW(2,1)/Wroad;  
CLncsd = (1/Wsd)*CLncW(3,1)/Wroad;  

% Plots of controllers, closed-loop with weights, closed-loop w/out weights
figure(1);
bodemag(Kinf,'b',Kar,'r-.')
legend('Hinf','AR')
grid on;
if exist('garyfyFigure','file'); garyfyFigure; end

figure(2)
sigma(CLinfW,'b',CLncW,'g--',CLarW.Nominal,'r-.',{1,140});
grid on;
legend('Hinf','NC','AR','location','SouthEast')
ylim([-40 0]);
title('Weighted closed-loops')
if exist('garyfyFigure','file'); garyfyFigure; end

figure(3)
subplot(2,1,1)
bodemag(qcar('ab','r'),'k:', CLinf('ab','r').Nominal,'b',CLncab,'g--',...
    CLar('ab','r').Nominal,'r-.',{1,140});
grid on;
title('') %title('Road to body acceleration')
if exist('garyfyFigure','file'); garyfyFigure; end

subplot(2,1,2)
bodemag(qcar('sd','r'),'k:', CLinf('sd','r').Nominal,'b',CLncsd,'g--',...
    CLar('sd','r').Nominal,'r-.',{1,140});
grid on;
legend('OL','CLinf','CLnc','CLar','location','SouthEast')
title(''); %title('Road to susp. deflection')
if exist('garyfyFigure','file'); garyfyFigure; end

%% Robust Mu Design
% Compute a robust additive regret controller and compare performance.

% Robust additive regret
[Krob,CLrobW,GAM] = robregretsyn(qcaric,nmeas,ncont,{gamdInt,gamJ},Tol); 
gamrob = GAM(1)

% Closed-loop model (without design weights)
Krob.u = {'sd','ab'};  Krob.y = 'u';
CLrob = connect(qcar,Act,Krob,'r',{'xb';'sd';'ab'; 'fs'});

% Compare controllers and closed-loops
figure(4);
bodemag(Kinf,'b',Kar,'r-.',Krob,'m:')
legend('Hinf','AR','Rob')
grid on;
if exist('garyfyFigure','file'); garyfyFigure; end

figure(5)
bodemag(qcar('ab','r'),'k:', CLinf('ab','r').Nominal,'b',CLncab,'g--',...
    CLar('ab','r').Nominal,'r-.',CLrob('ab','r').Nominal,'m:',{1,140});
grid on;
legend('OL','CLinf','CLnc','CLar','CLrob','location','SouthEast')
title('Road to body acceleration')
if exist('garyfyFigure','file'); garyfyFigure; end


% Time-Domain Evaluation
% Perform time-domain simulations using a road disturbance signal $r(t)$ 
% representing a road bump of height 5 cm. Compare nominal and robust
% additiv regret designs.
t = 0:Ts:1;
roaddist = zeros(size(t));
roaddist(1:101) = 0.025*(1-cos(8*pi*t(1:101)));

figure(6);
nsamp =  50;
lsim(usample(CLar,nsamp),'b',CLar.Nominal,'r',roaddist,t)
title('Kar'); grid on;
legend('Perturbed','Nominal','location','SouthEast')
if exist('garyfyFigure','file'); garyfyFigure; end

figure(7);
lsim(usample(CLrob,nsamp),'b',CLrob.Nominal,'r',roaddist,t)
title('Krob'); grid on;
legend('Perturbed','Nominal','location','SouthEast')
if exist('garyfyFigure','file'); garyfyFigure; end

%% Paper Plots

% Nominal frequency domain response
w = linspace(1,140,1e3);
figure(21)
clf
sys = {qcar('ab','r'), CLar('ab','r').Nominal, ...
    CLrob('ab','r').Nominal};
sysc = {'k:','b','r-.'};
for i=1:3
    m = bode(sys{i},w);
    mdB = 20*log10(m(:));
    semilogx(w,mdB,sysc{i}); 
    hold on;
end
hold off
grid on;
xlabel('Frequency, rad/sec');
ylabel('Magnitude, dB');
legend('OL','CL-AR','CL-Rob','location','SouthEast')
title('Road to body acceleration')
if exist('garyfyFigure','file'); garyfyFigure; end
exportgraphics(gcf, 'ActiveSusp_FreqDomain.pdf');

% Time Domain: Nominal AR controller 
figure(22); 
clf
t = 0:Ts:1;
roaddist = zeros(size(t));
roaddist(1:101) = 0.025*(1-cos(8*pi*t(1:101)));
nsamp =  50;
for i=1:nsamp    
    [y,t] = lsim(usample(CLar),roaddist,t);
    subplot(2,1,1)
    ph1 = plot(t,y(:,3),'b-.'); hold on;
    subplot(2,1,2)
    plot(t,y(:,2),'b-.'); hold on;
end
[y,t] = lsim(CLar.Nominal,roaddist,t);
subplot(2,1,1)
ph2 = plot(t,y(:,3),'c'); hold on;
legend([ph1; ph2],'Perturbed','Nominal','location','Northeast')
grid on; ylabel('a_b, m/s^2');
axis([t(1) t(end) -10 11.5]);

subplot(2,1,2)
plot(t,y(:,2),'c'); hold on;
grid on; xlabel('Time, sec'); ylabel('s_d, m'); 
axis([t(1) t(end) -0.05 0.05]);
if exist('garyfyFigure','file'); garyfyFigure; end
exportgraphics(gcf, 'ActiveSusp_NominalTimeDomain.pdf');

% Robust frequency domain response
% figure(23)
% clf
% sys = {qcar('ab','r'), CLar('ab','r').Nominal, CLrob('ab','r').Nominal};
% sysc = {'k:','r-.','c'};
% for i=1:3
%     m = bode(sys{i},w);
%     mdB = 20*log10(m(:));
%     semilogx(w,mdB,sysc{i}); 
%     hold on;
% end
% hold off
% grid on;
% xlabel('Frequency, rad/sec');
% ylabel('Magnitude, dB');
% legend('OL','CL-AR','CL-Rob','location','SouthEast')
% title('Road to body acceleration')
% if exist('garyfyFigure','file'); garyfyFigure; end
% exportgraphics(gcf, 'ActiveSusp_RobustFreqDomain.pdf');

% Time Domain: Robust AR controller
figure(24)
clf
for i=1:nsamp    
    [y,t] = lsim(usample(CLrob),roaddist,t);
    subplot(2,1,1)
    ph1 = plot(t,y(:,3),'r-.'); hold on;
    subplot(2,1,2)
    plot(t,y(:,2),'r-.'); hold on;
end
[y,t] = lsim(CLrob.Nominal,roaddist,t);
subplot(2,1,1)
ph2 = plot(t,y(:,3),'g'); hold on;
legend([ph1; ph2],'Perturbed','Nominal','location','Northeast')
grid on; ylabel('a_b, m/s^2');
axis([t(1) t(end) -10 11.5]);

subplot(2,1,2)
plot(t,y(:,2),'g'); hold on;
grid on; xlabel('Time, sec'); ylabel('s_d, m'); 
axis([t(1) t(end) -0.05 0.05]);
if exist('garyfyFigure','file'); garyfyFigure; end
exportgraphics(gcf, 'ActiveSusp_RobustTimeDomain.pdf');
