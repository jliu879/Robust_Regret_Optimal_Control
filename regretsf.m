function F = regretsf(P,Nu,gam)
% REGRETSF Spectral factorization for regret-based design
%
%  F = REGRETSF(P,NCON,GAM) calculates the spectral factor F associated
%  with the regret bound defined by GAM = [gamd, gamJ]. The spectral factor
%  is square, stable, and with a stable inverse. It also satisfies:
%     F'F =  gamd^2 I + gamJ^2 CL' CL
%  where CL is the optimal closed-loop from disturbance to error with
%  the optimal non-causal controller.  P is a discrete-time LTI plant with
%  state-space equations where dx=x(t+1):
%       x' =  A x +  Bd d +  Bu u
%        e = Ce x         + Deu u
%  NCON specifies the number of controls u where the inputs of P are
%  ordered as [d; u].
%

% Dimensions
[Ne,Ndu] = size(P);
Nd = Ndu-Nu;
Nx = size(P.A,1);

% Extract regret bound constants
gamd = gam(1);
gamJ = gam(2);

% Compute closed-loop from d to e with optimal non-causal controller
[K,CL] = ncsyn(P,Nu);

% Create state-space matrices for CLhat = [gamJ*CL; gamd*I]
% Note that CLhat'CLhat = gamd^2 I + gamJ^2 CL' CL
[Ahat,Bhat,Ccl,Dcl] = ssdata(CL);
Chat = [gamJ*Ccl; zeros(Nd,2*Nx)];
Dhat = [zeros(Ne,Nd); gamd*eye(Nd)];

% Construct spectral factor
Ts = P.Ts;
if gamJ == 0
    F = ss( gamd*eye(Nd) );
    F.Ts = Ts;
% elseif 1    
%     % Use built-in spectral factorization code from Matlab.
%     % This will not yield a minimal (order Nx) realization, in general.
%     % However, it seems to be a better conditioned numerical method.
%     [GG,SS] = spectralfact(gamd^2*eye(Nd)+gamJ^2*CL'*CL); 
%     F = sqrtm(SS)*GG;
else
    % Solve first DARE. This factorizes the zeros.
    Qhat = Chat'*Chat;
    Rhat = Dhat'*Dhat;
    Shat = Chat'*Dhat;
    [Xhat,Kxhat] = idare(Ahat,Bhat,Qhat,Rhat,Shat);
    Hhat = Rhat+Bhat'*Xhat*Bhat;

    % Solve second DARE. This factorizes the poles.
    [Yhat,Khaty] = idare(Ahat',Kxhat',0,inv(Hhat));
    What = inv(Hhat)+Kxhat*Yhat*Kxhat';

    % Define spectral factor
    AF = Ahat-Khaty'*Kxhat;
    BF = Bhat-Khaty';
    DF = sqrtm( inv(What) );
    CF = DF*Kxhat;

    % Remove unobservable states 
    % Note: The minimal realization can be computed directly but it
    % requires a few more steps.
    AF2 = AF(Nx+1:end,Nx+1:end);
    BF2 = BF(Nx+1:end,:);
    CF2 = CF(:,Nx+1:end);
    F = ss(AF2,BF2,CF2,DF,Ts);
end

% Check spectral factorization
% norm(F'*F - (gamd^2*eye(Nd)+gamJ^2*CL'*CL) ,inf)

