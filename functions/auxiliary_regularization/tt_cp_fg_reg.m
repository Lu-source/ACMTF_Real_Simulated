function [f,G] = tt_cp_fg_reg(Z,A,Znormsqr,lamda1)%add the Tikhonov regularization
%TT_CP_FG Computes function and gradient of the CP function.
%
%   [F,G] = TT_CP_FG(Z,A) calculates F = (1/2) ||Z - ktensor(A)||^2 where
%   Z is an N-way tensor and A is a ktensor or a cell array with N
%   factor matrices. It also calculates the gradient of the CP fit
%   function where Z is an N-way tensor and A is a ktensor or a
%   cell array with N factor matrices. The result is also a cell
%   array with N factor matrices corresponding to the gradients; in
%   other words, G{n}(:,r) is the partial derivative of the fit
%   function with respect to A{n}(:,r). 
%
%   [F,G] = TT_CP_FG(Z,A,NORMZSQR) also passes in the pre-computed
%   norm of Z, which makes the computations faster. 
%
%   See also CP_OPT, TT_CP_FUN.
%
%MATLAB Tensor Toolbox. Copyright 2018, Sandia Corporation.
%lamda is the parameter for regularization


%% Set-up
% if ~isa(Z,'tensor') && ~isa(Z,'sptensor')
%     error('Z must be a tensor or a sptensor');
% end
N = ndims(Z);

if ~iscell(A) && ~isa(A,'ktensor');
    error('A must be a cell array or ktensor');
end

if isa(A,'ktensor')
    A = tocell(A);
end
R = size(A{1},2);

%% Upsilon and Gamma
Upsilon = cell(N,1);
for n = 1:N
    Upsilon{n} = A{n}'*A{n};
end

Gamma = cell(N,1);
for n = 1:N
    Gamma{n} = ones(R,R);
    for m = [1:n-1,n+1:N]
        Gamma{n} = Gamma{n} .* Upsilon{m};
    end
end


%% Calculation

%F1
if exist('Znormsqr','var')
    f_1 = Znormsqr;
else
    f_1 = norm(Z)^2;
end

%% Calculate gradient and F2
G = cell(N,1);
U = mttkrp(Z,A,1);
V = A{1} .* U;
f_2 = sum(V(:));
% G{1} = -U + A{1}*Gamma{1};
G{1} = -U + A{1}*Gamma{1}+lamda1*A{1};%change by Lu
for n = 2:N
    U = mttkrp(Z,A,n);
    G{n} = -U + A{n}*Gamma{n}+lamda1*A{n};%change by Lu
end

%F3
W = Gamma{1} .* Upsilon{1};
f_3 = sum(W(:));

%change by Lu
f_4=norm(A{1},'fro')^2;
for i=2:N
    f_4=f_4+norm(A{i},'fro')^2;
end    
%change by Lu

    

%SUM
% f = 0.5 * f_1 - f_2 + 0.5 * f_3;
f = 0.5 * f_1 - f_2 + 0.5 * f_3+lamda1/2*f_4;%change by Lu


