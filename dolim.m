function [L,Q,B,Gtau,C0,Ctau,ualpha,valpha,galpha,tau_decay_alpha,T_mode_oscil,varC0,varCtau] = dolim(x,tau_0)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHORS:
%
%   Christopher Cox (NOAA) christopher.j.cox@noaa.gov
%
%   Appreciate guidance from Dr. Cecile Penland and Python codes from 
%   Meg Fowler (https://github.com/NOAA-PSL/Linear_Inverse_Modeling)
%
% PURPOSE:
%
%   Builds a Linear Inverse Model. 
%
% INPUT:
%
% x is a NxD array where 
%   N = length of the time series
%   D = number of stations 
%
% tau_0 = covariance lag specified in units of time
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% center the time series
x = x-nanmean(x);

% x is transposed as
xT = transpose(x);

%   and has dimensions:
D = size(x,2); % number of positions
N = size(x,1); % actual length of the data set

% check that the time dimension exceeds the state vector dimension to prove
% that the input time series is correctly oriented
disp(' ');
if N < D 
    disp('The input data needs to be transposed. STOP.');
    return
else
    disp('Dimensions look good.');
end

% From Penland (1989) Sect. 5 The POP analysis: recipe

% 1) The theory has been presented...
% Do LIM:
% - The Fokker-Plank Equation describes the probability distribution of a
%   stochastic processes and the time-evolution of the moments can be derived
%   from it.
% - The Langevin equation is closely related.
%   generalized, dx/dt = -gv + F(t)
%   discretized, dxi/dt = sum(j) Lij*xj + sum(a)Sia*Ea, where L and S are
%   systematic and noise contributions and the expectation is the the time
%   scales of variability in L >> that of S.

% 2) Form the contemporaneous covariance matrix of the process...Eq (3)
C0 = NaN(D,D);  % the contemporaneous covariance matrix 
for k1 = 1:D

    for k2 = 1:D

        a = x(:,k1);
        b = transpose(xT(k2,:));         
        ind = find(isnan(a) | isnan(b));
        a(ind) = [];
        b(ind) = [];
        C0(k1,k2) = sum(a.*b) / (length(a)-1);
        varC0(k1,k2) = var(a);
        
    end

end

% 3) Project the centered data onto EOFs
% skipping this...

% 4) Calculate sample estimate of tau-lagged correlation matrix ...
Ctau = NaN(D,D); % the lagged covariance matrix
for k1 = 1:D

    for k2 = 1:D

%         a = x(tau_0+1:end,k1);
%         b = xT(k2,1:N-tau_0)';
%         ind = find(isnan(a) | isnan(b));
%         a(ind) = [];
%         b(ind) = [];
%         Ctau(k1,k2) = sum(a.*b) / (length(a)-tau_0-1);
%         varCtau(k1,k2) = var(a);

        a = x(:,k1);
        b = transpose(xT(k2,:));
        a = a(tau_0+1:end);
        b = b(1:end-tau_0);
        Ctau(k1,k2) = nansum(a.*b) / (length(find(~isnan(a) & ~isnan(b)))-tau_0-1);
        varCtau(k1,k2) = nanvar(a);


    end

end
% ...and construct the Green function at lag tau, G(tau)
%    G is related to L in the Langevin equation as follows:
%    G(tau) = exp(L*tau)
%    Thus, Gtau is X in the equation AX = B: i.e., we treat the stations
%    as a system of linear equations relating the contemporaneous and
%    lagged covariances matrices. 
Gtau = Ctau * C0^-1;

% The Cayley-Hamilton Theorem allows us to put the eigenvectors,
% eigenvalues, and adjoints into standard linear algebra theory,
%    Av = lv where A is a square matrix, l are the eigenvalues, and v is the eigenvector
% i.e., (all substripted alpha, the generlized form of specific case, tau)
%    L*u = u*B and L'*v = v*B
% and
%    g = exp(B*tau)

% get eigenvalues, g, and eigenvectors, u.
[ualpha,galpha] = eig(Gtau);
galpha = diag(galpha);
[~,galpha_sort_ind] = sort(galpha,'descend');
galpha = galpha(galpha_sort_ind);
ualpha = ualpha(:,galpha_sort_ind);

% get the adjoints, v
[valpha,galpha_adj] = eig(transpose(Gtau));
[~,valpha_sort_ind] = sort(diag(galpha_adj),'descend');
valpha = valpha(:,valpha_sort_ind);

% resid = ualpha * valpha'; % if already normalized, this will be 1
% valpha = valpha / resid;

% Compute beta then normalize by tau to generalize
%   Note that Balpha all have negative real psrts and like modes ualpha and
%   adjoints valpha, are either real or complex conjugates.
Btau = log(galpha);
Balpha = Btau / tau_0;

% Calculate the decay time 
tau_decay_alpha = -1./real(Balpha);

% Order the modes, descending 
[tau_decay_alpha,ind_descend] = sort(tau_decay_alpha,'descend');
ualpha = ualpha(:,ind_descend);
valpha = valpha(:,ind_descend);
galpha = galpha(ind_descend);
Balpha = Balpha(ind_descend);

% Make diagonal array of beta
B = diag(Balpha); % the complex values matter for transpose...something got mixed up above?

% normalize modes and adjoints so they are 1 
% set alpha = beta do the sum on C-H slide: should be 1, if not use the result to scale one of them
resid = transpose(ualpha) * valpha; % if already normalized, this will be 1
un = ualpha / resid;
% check that the normalization is correct, benchamrking again ~single prec roundoffs
if ~any(abs(real(diag(transpose(un)*valpha))-1) < 10^-8)
    disp('u,v are NOT normalized correctly.');
else
    disp('u,v are properly normalized.');
end    
ualpha = un; % set u to normalized u

% recalculate galpha using B
%galpha = exp(diag(B)*tau_0);

% Calculate the oscillation period
T_mode_oscil = diag(2*pi ./ imag(B));

% Now we can calculate a system (alpha) L 
L = ualpha * B * transpose(valpha);

% L should be real, but in practice there may be some roundoff in imag
% space so this checks that the imag component is smaller than single
% precision.
if ~any(imag(L) > 10^-4)
    L = real(L);
else
    disp('L is complex');
    %L = real(L); %return
end

% Eq. A9 Land paper. Fluctuation-dissipation relation.
Q = -1* (L*C0 + C0*transpose(L));

% One more check. This, by definition is 0
Resid = L*transpose(C0) + transpose(C0)*transpose(L) + Q;
if any(Resid > 10^-6)
    disp('Residual not 0. Something went wrong.');
end

disp('All done. Happy Limming :)');
disp(' ');






