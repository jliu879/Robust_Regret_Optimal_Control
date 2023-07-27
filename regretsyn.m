function [K,CL,GAM] = regretsyn(P,Ny,Nu,gamtry,tol)
% REGRETSYN Synthesis of full-information, non-causal controller
%
%  [K,CL] = REGRETSYN(P,NMEAS,NCON,GAMTRY) calculates a regret-based
%  controller u = K y for the discrete-time, LTI plant P with state-space
%  equations where dx=x(t+1):
%        dx =  A x +  Bd d +  Bu u
%         e = Ce x         + Deu u
%         y = Cy x + Dyd d +
%  NMEAS and NCON specify the numbers of measurements y and controls u
%  (u and y must be the last inputs and outputs of P). The controller K
%  achieves the regret bound defined by GAMTRY = [gamd, gamJ]. K is
%  returned as empty if no controller exists that achieves the bound.
%  CL = LFT(P,K) is the closed-loop transfer function from disturbance
%  d to error e.
%
%  [K,CL,GAM] = REGRETSYN(P,NMEAS,NCON, {gamdInt, gamJ}) will bisect along
%  gamd with gamJ held fixed.  gamdInt is a 1-by-2 array specifying the
%  initial bisection interval [gamdInt(1) gamdInt(2)]. GAM = [gamd, gamJ]
%  returns the best regret bound obtained via bisection. If no feasible
%  point is found then GAM = [inf, inf] and K, CL are returned as empty.
%
%  [K,CL,GAM] = REGRETSYN(P,NMEAS,NCON, {gamd, gamJInt}) will bisect along
%  gamJ with gamd fixed.
%
%  [K,CL,GAM] = REGRETSYN(P,NMEAS,NCON,GAMTRY,TOL) specifies the (optional)
%  absolute and relative bisection tolerances, TOL = [AbsTol RelTol].
%  TOL is optional and default values are AbsTol=1e-2 and RelTol=1e-3.
%
% See also: bisection, hinfsyn

% Check inputs
narginchk(4,5);

if iscell(gamtry)
    if nargin==4
        tol = [1e-2, 1e-3];   % [AbsTol, RelTol]
    end
    gamd = gamtry{1};
    gamJ = gamtry{2};
    if isscalar(gamd)
        % Bisect on gamJ
        gamJInt = gamJ;
        FUN = @(gamJ) LOCALregfeas(P,Ny,Nu,[gamd,gamJ]);
        [FinalInt,BestData] = bisection(FUN,gamJInt,tol);
        if isempty(BestData)
            BestData = {[],[],[inf inf]};
        else
            BestData{3} = [gamd FinalInt(2)];
        end
    elseif isscalar(gamJ)
        % Bisect on gamd
        gamdInt = gamd;
        FUN = @(gamd) LOCALregfeas(P,Ny,Nu,[gamd,gamJ]);
        [FinalInt,BestData] = bisection(FUN,gamdInt,tol);
        if isempty(BestData)
            BestData = {[],[],[inf inf]};
        else
            BestData{3} = [FinalInt(2) gamJ];
        end
    else
        error('If GAMTRY is a 1-by-2 cell then one entry must be a scalar.');
    end
else
    [~,BestData] = LOCALregfeas(P,Ny,Nu,gamtry);
end

% Store outputs
K = BestData{1};
CL = BestData{2};
GAM = BestData{3};


%% LOCAL Function:  Regret Feasbility Problem
function [FLAG,DATA] = LOCALregfeas(P,Ny,Nu,gamtry)

% Full-information plant (removing measurement channels)
Ne = size(P,1) - Ny;
Pfi = P(1:Ne,:);

% Compute spectral factor associated with regret bound
F = regretsf(Pfi,Nu,gamtry);

% Solve H-infinty problem scaled by F
Pf = P*blkdiag( inv(F), eye(Nu) );
gaminf = 1;
[K0,~,GAM0] = hinfsyn(Pf,Ny,Nu,gaminf);
if isinf(GAM0)
    FLAG = false;  % Not feasible
    DATA = {[],[],[inf inf]};
else
    FLAG = true;  % Feasible
    DATA = {K0, lft(P,K0), gamtry};
end