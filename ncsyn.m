function [K,CL,X] = ncsyn(P,Nu) 
% NCSYN Synthesis of full-information, non-causal controller
%
%  [K,CL,X] = NCSYN(P,NCON) calculates the optimal full-information,
%  non-causal controller for the discrete-time LTI plant P with 
%  state-space equations where dx=x(t+1):
%       x' =  A x +  Bd d +  Bu u
%        e = Ce x         + Deu u
%  NCON specifies the number of controls u where the inputs of P are
%  ordered as [d; u]. The full-information, non-causal controller K 
%  is assumed to have access to the state as well as the full (past, 
%  current, and future) disturbance sequence. 
%
%  K is the non-causal controller from [x;d] to u with v'=v(t+1):
%      v' = Ak v + Bk [x;d]
%       u = Ck v + Dk [x;d]
%  The boundary condition is v(infty)=0 and the non-causal controller is 
%  iterated backward in time, i.e.
%     v = inv(Ak) ( v' - Bk [x;d] )
%  The closed-loop from d to e is the non-causal system CL.  X is the
%  solution of a DARE used to construct K and CL.


% Get plant data
[A,B,Ce,De,Ts] = ssdata(P);

Nx = size(A,1);
Ndu = size(B,2);
Nd = Ndu-Nu;
Ne = size(Ce,1);

Bd = B(:,1:Nd);
Bu = B(:,Nd+1:end);
Ded = De(:,1:Nd);
Deu = De(:,Nd+1:end);
if norm(Ded)>0 
    error('This function assumes Ded=0 and Dyd=0');    
end
if Ts==0
    error('This function assumes the plant is discrete-time')
end

% Solve DARE for X and compute stabilizing gain Kx
Q = Ce'*Ce;
S = Ce'*Deu;
R = Deu'*Deu;
[X,Kx] = idare(A,Bu,Q,R,S);

% Compute other gains in the optimal non-causal controller
Kv = inv(R+Bu'*X*Bu)*Bu';
%Kd = Kv*X*Bd;

% Form optimal non-causal controller with inputs [x;d] and output u
% The optimal non-causal controller can be expressed in standard form:
%     v' = Ak v + Bk [x;d]
%     u = Ck v + Dk [x;d]
% where Dk = [-Kx,  -Kd+Kv*X*Bd] = [-Kx 0]
Ahat11 = A-Bu*Kx;
Ak = inv(Ahat11)';
Bk = [zeros(Nx) -X*Bd];
Ck = -Kv*Ak;
Dk = [-Kx, zeros(Nu,Nd)];
K = ss(Ak,Bk,Ck,Dk,Ts);

% Form closed-loop with optimal non-causal controller 
Acl = [Ahat11 -Bu*Kv*Ak; zeros(Nx) Ak];
Bcl = [Bd; -X*Bd];
Ccl = [Ce-Deu*Kx -Deu*Kv*Ak];
Dcl = zeros(Ne,Nd);
CL = ss(Acl,Bcl,Ccl,Dcl,Ts);
