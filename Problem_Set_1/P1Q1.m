clear all
%Question 1:

%Variables
fs=10; %sampling frequency (Hz)
N=16384; %length of data
T=N/fs; %total time (1/Hz = seconds) 
delf=1/T; %discretized frequency
fN=fs/2; %nyquist frequency

%a=spectra_test_practice(:,1); %create a vector of data numbers
%import spectral data as given from course website
load spectra_test_practice.dat;
a=spectra_test_practice-mean(spectra_test_practice); 

Len=2^12; %dictates the length of each bin, in this case each will be 4096
NN=N/Len; %number of bins that are created by the length L, in this case will be 4 bins
Saatemp=zeros(Len,1); %creates a vector of zeros, length L that will be populated by for loop
Saa_avg=zeros(Len,1); %creates a vector of zeros, length L that will be populated by for loop

for i=1:NN                  %index from 1 to the number of bin being created
    atemp=a(1+(i-1)*Len:i*Len); %for index 1, the bin will go from 1 to L. For index 2, the bin will go from L+1 to 2L 
    A=fft(atemp); %do fast fourier transform
    AA=conj(A).*A; %calculate the conjugate, has to be .* to make sure its element by element
    Saatemp=(AA)*(1/(Len*fs)); %gives discrete domains, and creates a vector for spectral peaks of data (equation given in notes)
    Saa_avg=Saa_avg+Saatemp; %write over and store the information for the fft from the for loop
end
Saa4=Saa_avg/NN; %Saa for case with 4 bins
fMaxRes=fs/Len; %max resolution frequency 
f4=linspace(0,fs-fMaxRes,Len); %makes the x-axis for frequency vectors

%--------------------------------------------------------------------------------------------------------------------------------------------------
Len=2^10; %dictates the length of each bin, in this case each will be 1024
NN=N/Len; %number of bins that are created by the length L, in this case will be 16 bins
Saatemp=zeros(Len,1); %creates a vector of zeros, length L that will be populated by for loop
Saa_avg=zeros(Len,1);

%i am doing autospectral density function, not fft
for i=1:NN                  %index from 1 to the number of bin being created
    atemp=a(1+(i-1)*Len:i*Len); %for index 1, the bin will go from 1 to L. For index 2, the bin will go from L+1 to 2L 
    A=fft(atemp); %do fast fourier transform
    AA=conj(A).*A; %calculate the conjugate, has to be .* to make sure its element by element
    Saatemp=(AA)*(1/(Len*fs)); %gives discrete domains, and creates a vector for spectral peaks of data (equation given in notes)
    Saa_avg=Saa_avg+Saatemp; %write over and store the information for the fft from the for loop
end
Saa16=Saa_avg/NN; %Saa for case with 16 bins
fMaxRes=fs/Len; %max resolution frequency 
f16=linspace(0,fs-fMaxRes,Len); %makes the x-axis for frequency vectors

%--------------------------------------------------------------------------------------------------------------------------------------------------
Len=N; %dictates the length of each bin, in this case each will be 16384
NN=N/Len; %number of bins that are created by the length L, in this case will be 1 bin
Saatemp=zeros(Len,1); %creates a vector of zeros, length L that will be populated by for loop
Saa_avg=zeros(Len,1);

for i=1:NN                  %index from 1 to the number of bin being created
    atemp=a(1+(i-1)*Len:i*Len); %for index 1, the bin will go from 1 to L. For index 2, the bin will go from L+1 to 2L 
    A=fft(atemp); %do fast fourier transform
    AA=A.*conj(A); %calculate the conjugate, has to be .* to make sure its element by element
    Saatemp=(AA)*(1/(Len*fs)); %gives discrete domains, and creates a vector for spectral peaks of data (equation given in notes)
    Saa_avg=Saa_avg+Saatemp; %write over and store the information for the fft from the for loop
end
Saa1=Saa_avg/NN; %Saa for case with 1 bin
fMaxRes=fs/Len; %max resolution frequency 
f1=linspace(0,fs-fMaxRes,Len); %makes the x-axis for frequency vectors

%--------------------------------------------------------------------------------------------------------------------------------------------------%--------------------------------------------------------------------------------------------------------------------------------------------------

%check to make sure that all inputs are correct
check=(sum(Saa_avg)*delf)/(var(a,1)); %should be equal to 1

figure()
plot(f1,Saa1,'DisplayName','Original data')
hold on
plot(f16,Saa16,'red','DisplayName','Ensemble averaged data')
title('Ensemble average and original data')
xlabel('f')
ylabel('Saa')
legend('show')
%hold on
%plot (f1,Saa1,'yellow');
xlim ([0 fN/3])

%this doesn't work because it only returns 3 unique values (because the graph is symmetric)
%to extract 
%want to find max of red because its max peaks are the only true ones, so choose Saa2
% [yaxis,xaxis]=sort(Saa16,'descend'); %sort values from Saa in descending order
% maxvalues=yaxis(1:6); %yaxis is the max values on the yaxis
% indexmax=xaxis(1:6); %indexmax is where the max yaxis values are indexed
% 
% amp_desired=((Saa1(indexmax))*2)/N
% freq_desired=f1(indexmax)

%vector of uncorrected spectral autodensity 
UnA=[1934,2032,1655,1055,1215,1119];

%find actual amplitudes from spectral autodensity
AmpReal=2.*(sqrt((UnA.*fs.*16384)))./16384
