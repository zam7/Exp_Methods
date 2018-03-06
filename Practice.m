fs = 10; %sampling frequency (Hz)
N = 16384; %length of data
T = N/fs; %total time (1/Hz = seconds) 
delf=1/T; %discretized frequency

import spectra_test(1).dat;
spectral_data = spectra_test(1);
%a=spectra_test_practice(:,1); %create a vector of data numbers
a=spectral_data-mean(spectral_data); %can i leave it like this or is it just for the check?

A=fft(a);
AA=conj(A).*A; %has to be .* to make sure its element by element


Saa=(1/(fs*N))*(AA); %spectral peaks of data (equation given in notes)

check=(sum(Saa)*delf)/(var(a,1)) %should be equal to 1





