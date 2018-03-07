clear all
%Question 1:

%Variables
fs=10; %sampling frequency (Hz)
N=16384; %length of data
T=N/fs; %total time (1/Hz = seconds) 
delf=1/T; %discretized frequency
fN=fs/2; %niquist frequency

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
    Saatemp=Saatemp+(AA)*(1/(Len*fs)); %gives discrete domains, and creates a vector for spectral peaks of data (equation given in notes)
    Saa_avg=Saa_avg+Saatemp; %write over and store the information for the fft from the for loop
end
Saa4=Saa_avg/NN; %Saa for case with 4 bins
fMaxRes=fs/Len; %max resolution frequency 
f4=linspace(0,fs-fMaxRes,Len); %makes the x-axis for frequency vectors

%--------------------------------------------------------------------------------------------------------------------------------------------------
Len=2^10; %dictates the length of each bin, in this case each will be 4096
NN=N/Len; %number of bins that are created by the length L, in this case will be 4 bins
Saatemp=zeros(Len,1); %creates a vector of zeros, length L that will be populated by for loop
Saa_avg=zeros(Len,1);
for i=1:NN                  %index from 1 to the number of bin being created
    atemp=a(1+(i-1)*Len:i*Len); %for index 1, the bin will go from 1 to L. For index 2, the bin will go from L+1 to 2L 
    A=fft(atemp); %do fast fourier transform
    AA=conj(A).*A; %calculate the conjugate, has to be .* to make sure its element by element
    Saatemp=(AA)*(1/(Len*fs)); %gives discrete domains, and creates a vector for spectral peaks of data (equation given in notes)
    Saa_avg=Saa_avg+Saatemp; %write over and store the information for the fft from the for loop
end
Saa2=Saa_avg/NN; %Saa for case with 2 bins
fMaxRes=fs/Len; %max resolution frequency 
f2=linspace(0,fs-fMaxRes,Len); %makes the x-axis for frequency vectors

%--------------------------------------------------------------------------------------------------------------------------------------------------

%check to make sure that all inputs are correct
check=(sum(Saatemp)*delf)/(var(a,1)); %should be equal to 1

figure()
loglog(f4,Saa4)
hold on
loglog(f2,Saa2, 'magenta')
xlim ([0 fN])