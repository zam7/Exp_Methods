%Question 1:

%Variables
fs = 10; %sampling frequency (Hz)
N = 16384; %length of data
T = N/fs; %total time (1/Hz = seconds) 
delf=1/T; %discretized frequency

%a=spectra_test_practice(:,1); %create a vector of data numbers
%import spectral data as given from course website
a=spectral_data-mean(spectral_data); %can i leave it like this or is it just for the check?

%do fast fourier transform
A=fft(a);
%calculate the conjugate
AA=conj(A).*A; %has to be .* to make sure its element by element

Saa=(1/(fs*N))*(AA); % create a vector for spectral peaks of data (equation given in notes)

%check to make sure that all inputs are correct
check=(sum(Saa)*delf)/(var(a,1)); %should be equal to 1

f=0:delf:fs-delf; % create a vector for frequency data on x axis of graph

% up to here, the data is shown in the highest resolution because 

d1=a(1:N);
d2=a(N:2*N);


