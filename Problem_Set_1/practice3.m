%Question 3:
%Part i--------------------------------------------------------------------
load jet1.dat
JD=jet1;
JDcalc=jet1-mean(jet1);
load Gaussian.txt
Gauss=Gaussian;
Gausscalc=Gaussian-mean(Gaussian);

%Variables
fs=100; %Hz
T=194.56; %seconds
N=19456; %length of data set
delf=1/T; %discretized frequency
fN=fs/2; %nyquist frequency


%Find mean, variance, and root mean squared for LDV Jet Data 
M=mean(JD);
V=var(JD,1);
RMS=sqrt(V); %rms is the same as standard deviation

%Create a histogram
% H=hist(JD);
% hist(Gaussian,500)
% hold on
% hist(JD,500) %number of bins = 500

%Create Spectral data
Len=512; %dictates the length of each bin, in this case each will be 512
NN=N/Len; %number of bins that are created by the length L, in this case will be 1 bin
Saatemp=zeros(Len,1); %creates a vector of zeros, length L that will be populated by for loop
Saa_avg=zeros(Len,1);
meanNN=zeros(1,NN);
vNN=zeros(1,NN);
RMSNN=zeros(1,NN);
for i=1:NN                  %index from 1 to the number of bin being created
    JDtemp=JDcalc(1+(i-1)*Len:i*Len); %for index 1, the bin will go from 1 to L. For index 2, the bin will go from L+1 to 2L 
    meanNN(i)=mean(JDtemp); %finds mean of first segment of data on first pass and will report numbers that populate a 1x38 array
    vNN(i)=var(JDtemp,1); %finds variance of data along each segment 
    RMSNN(i)=sqrt(vNN(i));
    A=fft(JDtemp); %do fast fourier transform
    AA=conj(A).*A; %calculate the conjugate, has to be .* to make sure its element by element
    Saatemp=(AA)*(1/(Len*fs)); %gives discrete domains, and creates a vector for spectral peaks of data (equation given in notes)
    Saa_avg=Saa_avg+Saatemp; %write over and store the information for the fft from the for loop
        %Autocorrelation function
        for j=1:Len
           Raa1(i,j)=0; %autocorrelation of sensor 1 to sensor 1 
           for t=1:Len
              k=t+(j-1);
              if k<=Len
                 k=k;
              else
                 k=k-Len;
              end
              Raa1(i,j)=Raa1(i,j)+JDtemp(t)*JDtemp(k);
           end
        end
end
Saa1=Saa_avg/NN; %Saa for case with 38 bins
fMaxRes=fs/Len; %max resolution frequency 
f=linspace(0,fs-fMaxRes,Len); %makes the x-axis for frequency vectors

Raa1=Raa1./Len; %normalize the Raa vector by length N
Raaavg=mean(Raa1,1);
Raaavg=Raaavg./max(Raaavg);
f1=linspace(0,Len/fs,Len); 

%check to make sure that all inputs are correct
check=(sum(Saa_avg)*delf)/(var(JDcalc,1)); %should be equal to 1




%Part iii---------------------------------------------------------

%try with data set of different length
%Create Spectral data
Len=1024; %dictates the length of each bin, in this case each will be 512
NN=N/Len; %number of bins that are created by the length L, in this case will be 1 bin
Saatemp=zeros(Len,1); %creates a vector of zeros, length L that will be populated by for loop
Saa_avg=zeros(Len,1);
meanNN=zeros(1,NN);
vNN=zeros(1,NN);
RMSNN=zeros(1,NN);
for i=1:NN                  %index from 1 to the number of bin being created
    JDtemp=JDcalc(1+(i-1)*Len:i*Len); %for index 1, the bin will go from 1 to L. For index 2, the bin will go from L+1 to 2L 
    meanNN(i)=mean(JDtemp); %finds mean of first segment of data on first pass and will report numbers that populate a 1x38 array
    vNN(i)=var(JDtemp,1); %finds variance of data along each segment 
    RMSNN(i)=sqrt(vNN(i));
    A=fft(JDtemp); %do fast fourier transform
    AA=conj(A).*A; %calculate the conjugate, has to be .* to make sure its element by element
    Saatemp=(AA)*(1/(Len*fs)); %gives discrete domains, and creates a vector for spectral peaks of data (equation given in notes)
    Saa_avg=Saa_avg+Saatemp; %write over and store the information for the fft from the for loop
        %Autocorrelation function
        for j=1:Len
           Raa2(i,j)=0; %autocorrelation of sensor 1 to sensor 1 
           for t=1:Len
              k=t+(j-1);
              if k<=Len
                 k=k;
              else
                 k=k-Len;
              end
              Raa2(i,j)=Raa2(i,j)+JDtemp(t)*JDtemp(k);
           end
        end
end
Saa2=Saa_avg/NN; %Saa for case with 19 bins
fMaxRes=fs/Len; %max resolution frequency 
f2=linspace(0,fs-fMaxRes,Len); %makes the x-axis for frequency vectors

Raa2=Raa2./Len; %normalize the Raa vector by length N
Raaavg2=mean(Raa2,1);
Raaavg2=Raaavg2./max(Raaavg2);
f3=linspace(0,Len/fs,Len);

%compare a and b-----------------------------------------------------------

%lot data with the mean (i.e. don't subtract it out)
load jet1.dat
JD=jet1;
JDcalc=jet1;
load Gaussian.txt
Gauss=Gaussian;
Gausscalc=Gaussian;

%Variables
fs=100; %Hz
T=194.56; %seconds
N=19456; %length of data set
delf=1/T; %discretized frequency
fN=fs/2; %nyquist frequency


%Find mean, variance, and root mean squared for LDV Jet Data 
Mm=mean(JD);
Vm=var(JD,1);
RMSm=sqrt(V); %rms is the same as standard deviation

%Create a histogram
% H=hist(JD);
% hist(Gaussian,500)
% hold on
% hist(JD,500) %number of bins = 500

%Create Spectral data
Len=512; %dictates the length of each bin, in this case each will be 512
NN=N/Len; %number of bins that are created by the length L, in this case will be 1 bin
Saatemp=zeros(Len,1); %creates a vector of zeros, length L that will be populated by for loop
Saa_avg=zeros(Len,1);
meanNN=zeros(1,NN);
vNN=zeros(1,NN);
RMSNN=zeros(1,NN);
for i=1:NN                  %index from 1 to the number of bin being created
    JDtemp=JDcalc(1+(i-1)*Len:i*Len); %for index 1, the bin will go from 1 to L. For index 2, the bin will go from L+1 to 2L 
    meanNN(i)=mean(JDtemp); %finds mean of first segment of data on first pass and will report numbers that populate a 1x38 array
    vNN(i)=var(JDtemp,1); %finds variance of data along each segment 
    RMSNN(i)=sqrt(vNN(i));
    A=fft(JDtemp); %do fast fourier transform
    AA=conj(A).*A; %calculate the conjugate, has to be .* to make sure its element by element
    Saatemp=(AA)*(1/(Len*fs)); %gives discrete domains, and creates a vector for spectral peaks of data (equation given in notes)
    Saa_avg=Saa_avg+Saatemp; %write over and store the information for the fft from the for loop
        %Autocorrelation function
        for j=1:Len
           Raam1(i,j)=0; %autocorrelation of sensor 1 to sensor 1 
           for t=1:Len
              k=t+(j-1);
              if k<=Len
                 k=k;
              else
                 k=k-Len;
              end
              Raam1(i,j)=Raam1(i,j)+JDtemp(t)*JDtemp(k);
           end
        end
end
Saam1=Saa_avg/NN; %Saa for case with 38 bins
fMaxRes=fs/Len; %max resolution frequency 
fm=linspace(0,fs-fMaxRes,Len); %makes the x-axis for frequency vectors

Raam1=Raam1./Len; %normalize the Raa vector by length N
Raamavg=mean(Raam1,1);
Raamavg=Raamavg./max(Raamavg);
fm1=linspace(0,Len/fs,Len); 


%All the plots-------------------------------------------------------------
figure()
hist(Gaussian,500) %make histogram of normally distributed data to compare to our data set
hold on
H=hist(JD); %make histogram of data from our data set
hist(JD,500) %number of bins = 500
title('Histogram')
xlabel('u (m/s)')
ylabel('Counts')

figure()
loglog(f,Saa1,'g','DisplayName','512 data length') %plot for spectral data with 512 length
title('Loglog autospectral density data')
xlabel('f (Hz)')
ylabel('Suu (m^2/s^3)')
hold on
loglog(f2,Saa2,'r','DisplayName','1024 data length') %plot for spectral data with 1024 length
legend('show')
hold off

figure()
plot(f1,Raaavg,'DisplayName','512 data length') %plot for autocorrelation with 512 length
title('Autocorrelation')
xlim([0 2.56]);
hold on
plot(f3,Raaavg2,'r','DisplayName','1024 data length') %plot for autocorrelation data with 1024 length
legend('show')
hold off

figure() %comparing spectral data without mean to data with mean
loglog(f,Saa1,'g','DisplayName','Data without mean') %plot for spectral data with 512 length without mean
title('Loglog autospectral density data')
xlabel('f (Hz)')
ylabel('Suu (m^2/s^3)')
hold on
loglog(fm,Saam1,'m','DisplayName','Data with mean') %plot for spectral data with 512 length with mean
legend('show')
hold off
