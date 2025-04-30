function out = waveletwrapper(input,time)
% C. Cox
% 09-28-2012
% 04-30-2024
% 
% Runs T&C wavelet codes for Bering Sea study.

% OPTIONS -----------------------------------------------------------------
anom    = 0           ;                   % Make the time series and anomaly.
rsc     = 0           ;                   % Remove the 1st harmonic (seasonal cycle from the time series).
pad     = 1           ;                   % pad the time series with zeroes (recommended)
dj      = 0.01        ;                   % sub-octaves: increased resolution in frequency space. Really only about time, but should not be larger than 0.5 for Morlet.
po2     = 15          ;                   % Powers of 2 (max frequency) (2^po2)/8 days
mother  = 'Morlet'    ;                   % Pick your favorite wavelet
units   = 'power'     ;                   % Look at 'power', 'amplitude', or 'phase' 
sa_sc   = -1          ;                   % Scale Ave. Frequencies if desired. Else should = -1.
% -------------------------------------------------------------------------

% AUTOMATIC OPTIONS AND LOAD ----------------------------------------------
data = input'         ;

% ANOMALY OPTION ----------------------------------------------------------
if rsc == 1 && anom == 1  
    y = harmfit(time,data,365.25,3,1,0);    
    data = data-y;
elseif rsc == 0 && anom == 1
    data = data-mean(data);
end
% -------------------------------------------------------------------------

% AUTO-OPTIONS ------------------------------------------------------------
n      = length(time);
dt     = time(2)-time(1) ;
s0     = 2*dt;                 % Nyquist
j1     = po2/dj;               % do po2 Powers-of-Two with dj suboctaves each 
[Cdelta dj0 gamma upsilon] = wavelet_getconstants(mother); % Constants from Torrence and Compo (1998) table 2.
% -------------------------------------------------------------------------
disp('Data Load Complete.');


% CHECKS ------------------------------------------------------------------
if dj > dj0 
    error(['Your dj is too large for the ',mother,' wavelet.']);
end
if s0 ~= 2*dt || 2^po2 < n % only report statistics if range of scale examines should include enough to reporiduce time series vis reconstruction
    noCheck = 1;
    fprintf('I am not checking for convergence of time series reconstruction\nusing inverse filter because you are not using\na large enough po2 to encompass all low frequencies.\n')
else
    noCheck = 0;
end
% -------------------------------------------------------------------------


% THE WAVELET TRANSFORM ---------------------------------------------------
[wave,period,scale,coi,X] = wavelet(data,dt,pad,dj,s0,j1,mother,upsilon,Cdelta,noCheck,0,6);
power     = (abs(wave)).^2;               % compute wavelet power spectrum
amplitude = abs(wave);                    % compute the wavelet amplitude spectrum
phase     = atan(imag(wave)./real(wave)); % compute the wavelet phase spectrum
if strcmpi(units,'power');
    spectra = power;
elseif strcmpi(units,'amplitude');
    spectra = amplitude;
elseif strcmpi(units,'phase'); 
    spectra = phase;
end
% -------------------------------------------------------------------------

% % OUTPUT ----------------------------------------------------------------
out.spectra = spectra;
out.x = X;
out.period = period;
% -------------------------------------------------------------------------