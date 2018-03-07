% Lab 1
% Zoe Maisel

% Columns for data
%  1   Ensemble counter                 (1-16777216)
%  2   Status
%  3   Velocity (Beam1|X)               (m/s)
%  4   Velocity (Beam2|Y)               (m/s)
%  5   Velocity (Beam3|Z)               (m/s)
%  6   Velocity (Beam4|Z2)              (m/s) %data is more gaussian so
%  only work with Z2
%  7   Amplitude (Beam1)                (counts)
%  8   Amplitude (Beam2)                (counts)
%  9   Amplitude (Beam3)                (counts)
% 10   Amplitude (Beam4)                (counts)
% 11   SNR (Beam1)                      (dB)
% 12   SNR (Beam2)                      (dB)
% 13   SNR (Beam3)                      (dB)
% 14   SNR (Beam4)                      (dB)
% 15   Correlation (Beam1)              (%)
% 16   Correlation (Beam2)              (%)
% 17   Correlation (Beam3)              (%)
% 18   Correlation (Beam4)              (%)


%Reynolds number for turbulent jet
v=0.010034; %cm2/s at 20 degrees Celcius
Djet=0.595; %cm 
Vjet=322.25; %cm/s, found by Q=VA of jet
Re=(Vjet*Djet)/v;1   

%row, column, dimesion
tic();
data(:,:,1) = load('30.9556cm.dat');0.595; %cm
data(:,:,2) = load('31.2604.dat');
data(:,:,3) = load('31.5652cm.dat'); %looks most like centerline where we started
data(:,:,4) = load('31.87cm.dat');  
data(:,:,5) = load('32.1748cm.dat');
data(:,:,6) = load('32.4796cm.dat');
data(:,:,7) = load('32.7844.dat');
data(:,:,8) = load('33.394cm.dat');
data(:,:,9) = load('34.0036cm.dat');
data(:,:,10) = load('34.6132cm.dat');
data(:,:,11) = load('35.2228cm.dat');
data(:,:,12) = load('35.8324cm.dat');
data(:,:,13) = load('36.442cm.dat');
data(:,:,14) = load('37.0516cm.dat'); %has symmetric bumps at end from resonance bouncing back
data(:,:,15) = load('37.6612cm.dat'); %has symmetric bumps at end from resonance bouncing back
data(:,:,16) = load('38.2708cm.dat'); %has symmetric bumps at end from resonance bouncing back
data(:,:,17) = load('38.8804cm.dat');
% heights that data was gathered at 31.5652,31.87,32.1748,32.4796,32.7844,33.394,34.0036,34.6132,35.2228,35.8324,36.442,37.0516,37.6612,38.2708,38.8804];
D=17; %dimension of data is 17 because there are 17 data sets
LD=length(data(:,:,1)); %vector to make length of data
toc(); %loading the data takes ~7 sec

%find mean of unfiltered data
for i = 1:D
    xmean(i)= mean(data(:,3,i));
    ymean(i)= mean(data(:,4,i));
    zmean1(i) = mean(data(:,5,i));
    zmean2(i) = mean(data(:,6,i));
    zmean(i) = mean(data(:,6,i));
    corrx(i) = mean(data(:,15,i));
end

%find perturbations of unfiltered data

for i = 1:D 
    for j = 1:LD
     xprime(j,1,i) = data(j,3,i)-xmean(i);
     yprime(j,1,i) = data(j,4,i)-ymean(i);
     zprime1(j,1,i) = data(j,5,i) -zmean1(i);
     zprime2(j,1,i) = data(j,6,i) -zmean2(i);
     zprime(j,1,i) = zprime2(j,1,i);
    end
end

U0=xmean(3); %based on the graph, we think that the third data point that we have is the actual midpoint
Uhalf=U0/2;
%from gaussian data

r12=3.0811; %radial distance of r(1/2) from centerline
H = [30.9556, 31.2604, 31.5652,31.87,32.1748,32.4796,32.7844,33.394,34.0036,34.6132,35.2228,35.8324,36.442,37.0516,37.6612,38.2708,38.8804]; % z value heights that all measurements were taken at
center = 31.5652; %cm centerline of jet, chosen from gaussian plot of maximum
Z = H-center; %convert height distances from bottom of bed to radial from centerline
ZoR= Z./r12; %z/r1/2

%the 3rd, 10th, and 15th data sets are the ones of interest for centerline,
%z/r(1/2)=1 and z/r(1/2)=2
xpZoR0=xprime(:,1,3); %ZoR1=0 occurs during the 3rd data set
xpZoR1=xprime(:,1,10); %ZoR1=1 occurs during the 10th data set
xpZoR2=xprime(:,1,15); %ZoR1=2 occurs during the 15th data set

ypZoR0=yprime(:,1,3); %ZoR1=0 occurs during the 3rd data set
ypZoR1=yprime(:,1,10); %ZoR1=1 occurs during the 10th data set
ypZoR2=yprime(:,1,15); %ZoR1=2 occurs during the 15th data set

zpZoR0=zprime(:,1,3); %ZoR1=0 occurs during the 3rd data set
zpZoR1=zprime(:,1,10); %ZoR1=1 occurs during the 10th data set
zpZoR2=zprime(:,1,15); %ZoR1=2 occurs during the 15th data set

%-----------------------Spectral Data--------------------------------------
% worked with original data because the unfiltered data has the correct
% time scale

%----------At Centerline----------- 
% Creating spectra data for xprime %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs=100; %Hz
Nxc=length(xpZoR0); % N is the length of the xprime vector
Lenxc=2^10; %dictates the length of each bin, in this case each will be 1024
NNxc=Nxc/Lenxc; %number of bins that are created by the length L
Saatempxc=zeros(Lenxc,1); %creates a vector of zeros, length L that will be populated by for loop
Saa_avgxc=zeros(Lenxc,1);

% autospectral density function 
for i=1:NNxc                  %index from 1 to the number of bin being created
    atemp=xpZoR0(1+(i-1)*Lenxc:i*Lenxc); %for index 1, the bin will go from 1 to L. For index 2, the bin will go from L+1 to 2L 
    Axc=fft(atemp); %do fast fourier transform
    AAxc=conj(Axc).*Axc; %calculate the conjugate, has to be .* to make sure its element by element
    Saatempxc=(AAxc)*(1/(Lenxc*fs)); %gives discrete domains, and creates a vector for spectral peaks of data (equation given in notes)
    Saa_avgxc=Saa_avgxc+Saatempxc; %write over and store the information for the fft from the for loop
end
Saaxc=Saa_avgxc/NNxc; %Saa for case with 16 bins
fMaxResxc=fs/Lenxc; %max resolution frequency 
fxc=linspace(0,fs-fMaxResxc,Lenxc); %makes the x-axis for frequency vectors

% Creating spectra data for yprime %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs=100; %Hz
Nyc=length(ypZoR0); % N is the length of the yprime vector
Lenyc=2^10; %dictates the length of each bin, in this case each will be 1024
NNyc=Nyc/Lenyc; %number of bins that are created by the length L
Saatempyc=zeros(Lenyc,1); %creates a vector of zeros, length L that will be populated by for loop
Saa_avgyc=zeros(Lenyc,1);

% autospectral density function 
for i=1:NNyc                  %index from 1 to the number of bin being created
    atempyc=ypZoR0(1+(i-1)*Lenyc:i*Lenyc); %for index 1, the bin will go from 1 to L. For index 2, the bin will go from L+1 to 2L 
    Ayc=fft(atempyc); %do fast fourier transform
    AAyc=conj(Ayc).*Ayc; %calculate the conjugate, has to be .* to make sure its element by element
    Saatempyc=(AAyc)*(1/(Lenyc*fs)); %gives discrete domains, and creates a vector for spectral peaks of data (equation given in notes)
    Saa_avgyc=Saa_avgyc+Saatempyc; %write over and store the information for the fft from the for loop
end
Saayc=Saa_avgyc/NNyc; %Saa for case with NNy bins
fMaxResyc=fs/Lenyc; %max resolution frequency 
fyc=linspace(0,fs-fMaxResyc,Lenyc); %makes the y-axis for frequency vectors

% Creating spectra data for zprime %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs=100; %Hz
Nzc=length(zpZoR0); % N is the length of the zprime vector
Lenzc=2^10; %dictates the length of each bin, in this case each will be 1024
NNzc=Nzc/Lenzc; %number of bins that are created by the length L
Saatempzc=zeros(Lenzc,1); %creates a vector of zeros, length L that will be populated by for loop
Saa_avgzc=zeros(Lenzc,1);

% autospectral density function 
for i=1:NNzc                  %index from 1 to the number of bin being created
    atempzc=zpZoR0(1+(i-1)*Lenzc:i*Lenzc); %for index 1, the bin will go from 1 to L. For index 2, the bin will go from L+1 to 2L 
    Azc=fft(atempzc); %do fast fourier transform
    AAzc=conj(Azc).*Azc; %calculate the conjugate, has to be .* to make sure its element by element
    Saatempzc=(AAzc)*(1/(Lenzc*fs)); %gives discrete domains, and creates a vector for spectral peaks of data (equation given in notes)
    Saa_avgzc=Saa_avgzc+Saatempzc; %write over and store the information for the fft from the for loop
end
Saazc=Saa_avgzc/NNzc; %Saa for case with NNz bins
fMaxReszc=fs/Lenzc; %max resolution frequency 
fzc=linspace(0,fs-fMaxReszc,Lenzc); %makes the z-axis for frequency vectors

%----------At Z/r(1/2)=1----------- 
% Creating spectra data for xprime %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs=100; %Hz
Nx1=length(xpZoR1); % N is the length of the xprime vector
Lenx1=2^10; %dictates the length of each bin, in this case each will be 1024
NNx1=Nx1/Lenx1; %number of bins that are created by the length L
Saatempx1=zeros(Lenx1,1); %creates a vector of zeros, length L that will be populated by for loop
Saa_avgx1=zeros(Lenx1,1);

% autospectral density function 
for i=1:NNx1                  %index from 1 to the number of bin being created
    atempx1=xpZoR1(1+(i-1)*Lenx1:i*Lenx1); %for index 1, the bin will go from 1 to L. For index 2, the bin will go from L+1 to 2L 
    Ax1=fft(atempx1); %do fast fourier transform
    AAx1=conj(Ax1).*Ax1; %calculate the conjugate, has to be .* to make sure its element by element
    Saatempx1=(AAx1)*(1/(Lenx1*fs)); %gives discrete domains, and creates a vector for spectral peaks of data (equation given in notes)
    Saa_avgx1=Saa_avgx1+Saatempx1; %write over and store the information for the fft from the for loop
end
Saax1=Saa_avgx1/NNx1; %Saa for case with 16 bins
fMaxResx1=fs/Lenx1; %max resolution frequency 
fx1=linspace(0,fs-fMaxResx1,Lenx1); %makes the x-axis for frequency vectors

% Creating spectra data for yprime %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs=100; %Hz
Ny1=length(ypZoR1); % N is the length of the yprime vector
Leny1=2^10; %dictates the length of each bin, in this case each will be 1024
NNy1=Ny1/Leny1; %number of bins that are created by the length L
Saatempy1=zeros(Leny1,1); %creates a vector of zeros, length L that will be populated by for loop
Saa_avgy1=zeros(Leny1,1);

% autospectral density function 
for i=1:NNy1                  %index from 1 to the number of bin being created
    atempy1=ypZoR1(1+(i-1)*Leny1:i*Leny1); %for index 1, the bin will go from 1 to L. For index 2, the bin will go from L+1 to 2L 
    Ay1=fft(atempy1); %do fast fourier transform
    AAy1=conj(Ay1).*Ay1; %calculate the conjugate, has to be .* to make sure its element by element
    Saatempy1=(AAy1)*(1/(Leny1*fs)); %gives discrete domains, and creates a vector for spectral peaks of data (equation given in notes)
    Saa_avgy1=Saa_avgy1+Saatempy1; %write over and store the information for the fft from the for loop
end
Saay1=Saa_avgy1/NNy1; %Saa for case with NNy bins
fMaxResy1=fs/Leny1; %max resolution frequency 
fy1=linspace(0,fs-fMaxResy1,Leny1); %makes the y-axis for frequency vectors

% Creating spectra data for zprime %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs=100; %Hz
Nz1=length(zpZoR1); % N is the length of the zprime vector
Lenz1=2^10; %dictates the length of each bin, in this case each will be 1024
NNz1=Nz1/Lenz1; %number of bins that are created by the length L
Saatempz1=zeros(Lenz1,1); %creates a vector of zeros, length L that will be populated by for loop
Saa_avgz1=zeros(Lenz1,1);

% autospectral density function 
for i=1:NNz1                  %index from 1 to the number of bin being created
    atempz1=zpZoR1(1+(i-1)*Lenz1:i*Lenz1); %for index 1, the bin will go from 1 to L. For index 2, the bin will go from L+1 to 2L 
    Az1=fft(atempz1); %do fast fourier transform
    AAz1=conj(Az1).*Az1; %calculate the conjugate, has to be .* to make sure its element by element
    Saatempz1=(AAz1)*(1/(Lenz1*fs)); %gives discrete domains, and creates a vector for spectral peaks of data (equation given in notes)
    Saa_avgz1=Saa_avgz1+Saatempz1; %write over and store the information for the fft from the for loop
end
Saaz1=Saa_avgz1/NNz1; %Saa for case with NNz bins
fMaxResz1=fs/Lenz1; %max resolution frequency 
fz1=linspace(0,fs-fMaxResz1,Lenz1); %makes the z-axis for frequency vectors

%----------At z/r(1/2)=2-----------
% Creating spectra data for xprime %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs=100; %Hz
Nx2=length(xpZoR2); % N is the length of the xprime vector
Lenx2=2^10; %dictates the length of each bin, in this case each will be 1024
NNx2=Nx2/Lenx2; %number of bins that are created by the length L
Saatempx2=zeros(Lenx2,1); %creates a vector of zeros, length L that will be populated by for loop
Saa_avgx2=zeros(Lenx2,1);

% autospectral density function 
for i=1:NNx2                  %index from 1 to the number of bin being created
    atempx2=xpZoR2(1+(i-1)*Lenx2:i*Lenx2); %for index 1, the bin will go from 1 to L. For index 2, the bin will go from L+1 to 2L 
    Ax2=fft(atempx2); %do fast fourier transform
    AAx2=conj(Ax2).*Ax2; %calculate the conjugate, has to be .* to make sure its element by element
    Saatempx2=(AAx2)*(1/(Lenx2*fs)); %gives discrete domains, and creates a vector for spectral peaks of data (equation given in notes)
    Saa_avgx2=Saa_avgx2+Saatempx2; %write over and store the information for the fft from the for loop
end
Saax2=Saa_avgx2/NNx2; %Saa for case with 16 bins
fMaxResx2=fs/Lenx2; %max resolution frequency 
fx2=linspace(0,fs-fMaxResx2,Lenx2); %makes the x-axis for frequency vectors

% Creating spectra data for yprime %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs=100; %Hz
Ny2=length(ypZoR2); % N is the length of the yprime vector
Leny2=2^10; %dictates the length of each bin, in this case each will be 1024
NNy2=Ny2/Leny2; %number of bins that are created by the length L
Saatempy2=zeros(Leny2,1); %creates a vector of zeros, length L that will be populated by for loop
Saa_avgy2=zeros(Leny2,1);

% autospectral density function 
for i=1:NNy2                  %index from 1 to the number of bin being created
    atempy2=xpZoR2(1+(i-1)*Leny2:i*Leny2); %for index 1, the bin will go from 1 to L. For index 2, the bin will go from L+1 to 2L 
    Ay2=fft(atempy2); %do fast fourier transform
    AAy2=conj(Ay2).*Ay2; %calculate the conjugate, has to be .* to make sure its element by element
    Saatempy2=(AAy2)*(1/(Leny2*fs)); %gives discrete domains, and creates a vector for spectral peaks of data (equation given in notes)
    Saa_avgy2=Saa_avgy2+Saatempy2; %write over and store the information for the fft from the for loop
end
Saay2=Saa_avgy2/NNy2; %Saa for case with NNy bins
fMaxResy2=fs/Leny2; %max resolution frequency 
fy2=linspace(0,fs-fMaxResy2,Leny2); %makes the y-axis for frequency vectors

% Creating spectra data for zprime %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs=100; %Hz
Nz2=length(zpZoR2); % N is the length of the zprime vector
Lenz2=2^10; %dictates the length of each bin, in this case each will be 1024
NNz2=Nz2/Lenz2; %number of bins that are created by the length L
Saatempz2=zeros(Lenz2,1); %creates a vector of zeros, length L that will be populated by for loop
Saa_avgz2=zeros(Lenz2,1);

% autospectral density function 
for i=1:NNz2                  %index from 1 to the number of bin being created
    atempz2=zpZoR2(1+(i-1)*Lenz2:i*Lenz2); %for index 1, the bin will go from 1 to L. For index 2, the bin will go from L+1 to 2L 
    Az2=fft(atempz2); %do fast fourier transform
    AAz2=conj(Az2).*Az2; %calculate the conjugate, has to be .* to make sure its element by element
    Saatempz2=(AAz2)*(1/(Lenz2*fs)); %gives discrete domains, and creates a vector for spectral peaks of data (equation given in notes)
    Saa_avgz2=Saa_avgz2+Saatempz2; %write over and store the information for the fft from the for loop
end
Saaz2=Saa_avgz2/NNz2; %Saa for case with NNz bins
fMaxResz2=fs/Lenz2; %max resolution frequency 
fz2=linspace(0,fs-fMaxResz2,Lenz2); %makes the z-axis for frequency vectors

figure()
axc = subplot(3,2,1);
plot(axc,fxc,Saaxc)
title(axc,'uprime at centerline')
ylabel(axc,'Saa')
axclog = subplot(3,2,2);
loglog(axclog,fxc,Saaxc)
title(axclog,'loglog uprime at centerline')
ylabel(axclog,'Saa')
ayc = subplot(3,2,3);
plot(ayc,fyc,Saayc)
title(ayc,'yprime at centerline')
ylabel(ayc,'Saa')
ayclog = subplot(3,2,4);
loglog(ayclog,fyc,Saayc)
title(ayclog,'loglog yprime at centerline')
ylabel(ayclog,'Saa')
azc = subplot(3,2,5);
plot(azc,fzc,Saazc)
title(azc,'zprime at centerline')
ylabel(azc,'Saa')
azclog = subplot(3,2,6);
loglog(azclog,fzc,Saazc)
title(azclog,'loglog zprime at centerline')
ylabel(azclog,'Saa')
xlabel('f')

figure()
ax1 = subplot(3,2,1);
plot(ax1,fx1,Saax1)
title(ax1,'uprime at z/rhalf=1')
ylabel(ax1,'Saa')
ax1log = subplot(3,2,2);
loglog(ax1log,fx1,Saax1)
title(ax1log,'loglog uprime at z/rhalf=1')
ay1 = subplot(3,2,3);
plot(ay1,fy1,Saay1)
title(ay1,'yprime at z/rhalf=1')
ylabel(ay1,'Saa')
ay1log = subplot(3,2,4);
loglog(ay1log,fy1,Saay1)
title(ay1log,'loglog yprime at z/rhalf=1')
az1 = subplot(3,2,5);
plot(az1,fz1,Saaz1)
title(az1,'zprime at z/rhalf=1')
ylabel(az1,'Saa')
az1log = subplot(3,2,6);
loglog(az1log,fz1,Saaz1)
title(az1log,'zprime at z/rhalf=1')
xlabel('f')

figure()
ax2 = subplot(3,2,1);
plot(ax2,fx2,Saax2)
title(ax2,'uprime at z/rhalf=2')
ylabel(ax2,'Saa')
ax2log = subplot(3,2,2);
loglog(ax2log,fx2,Saax2)
title(ax2log,'loglog uprime at z/rhalf=2')
ay2 = subplot(3,2,3);
plot(ay2,fy2,Saay2)
title(ay2,'yprime at z/rhalf=2')
ylabel(ay2,'Saa')
ay2log = subplot(3,2,4);
loglog(ay2log,fy2,Saay2)
title(ay2log,'loglog yprime at z/rhalf=2')
az2 = subplot(3,2,5);
plot(az2,fz2,Saaz2)
title(az2,'zprime at z/rhalf=2')
ylabel(az2,'Saa')
az2log = subplot(3,2,6);
loglog(az2log,fz2,Saaz2)
title(az2log,'loglog zprime at z/rhalf=2')
xlabel('f')

%------------Filter script from Cowen--------------------------------------
%first data filter 
filtertd1=data(:,3:6,1); %data that we are going to test filter with
    %take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the first data set
count1=data(:,1,1); %count of data all along length
filtertdMin=-1; %RANDOM (shouldbe centered around 0)
filtertdMax=1; %RANDOM (largest value in this set is 1.01)

[Data1,Time1]=agw_filter(filtertd1, count1, filtertdMin, filtertdMax);

%second data filter 
filtertd2=data(:,3:6,2); %data that we are going to test filter with
    %take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the second data set
count2=data(:,1,2); %count of data all along length

[Data2,Time2]=agw_filter(filtertd2, count2, filtertdMin, filtertdMax);

%third data filter ------ where we think the midpoint is
filtertd3=data(:,3:6,3); %data that we are going to test filter with
%take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the third data set
count3=data(:,1,3); %count of data all along length

[Data3,Time3]=agw_filter(filtertd3, count3, filtertdMin, filtertdMax);

%fourth data filter
filtertd4=data(:,3:6,4); %data that we are going to test filter with
%take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the fourth data set
count4=data(:,1,4); %count of data all along length

[Data4,Time4]=agw_filter(filtertd4, count4, filtertdMin, filtertdMax);

%fifth data filter
filtertd5=data(:,3:6,5); %data that we are going to test filter with
%take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the fifth data set
count5=data(:,1,5); %count of data all along length

[Data5,Time5]=agw_filter(filtertd5, count5, filtertdMin, filtertdMax);

%sixth data filter
filtertd6=data(:,3:6,6); %data that we are going to test filter with
%take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the sixth data set
count6=data(:,1,6); %count of data all along length

[Data6,Time6]=agw_filter(filtertd6, count6, filtertdMin, filtertdMax);

%seventh data filter
filtertd7=data(:,3:6,7); %data that we are going to test filter with
%take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the seventh data set
count7=data(:,1,7); %count of data all along length

[Data7,Time7]=agw_filter(filtertd7, count7, filtertdMin, filtertdMax);

%eigth data filter
filtertd8=data(:,3:6,8); %data that we are going to test filter with
%take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the eigth data set
count8=data(:,1,8); %count of data all along length

[Data8,Time8]=agw_filter(filtertd8, count8, filtertdMin, filtertdMax);

%ninth data filter
filtertd9=data(:,3:6,9); %data that we are going to test filter with
%take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the ninth data set
count9=data(:,1,9); %count of data all along length

[Data9,Time9]=agw_filter(filtertd9, count9, filtertdMin, filtertdMax);

%tenth data filter
filtertd10=data(:,3:6,10); %data that we are going to test filter with
%take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the tenth data set
count10=data(:,1,10); %count of data all along length

[Data10,Time10]=agw_filter(filtertd10, count10, filtertdMin, filtertdMax);

%eleventh data filter
filtertd11=data(:,3:6,11); %data that we are going to test filter with
%take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the eleventh data set
count11=data(:,1,11); %count of data all along length

[Data11,Time11]=agw_filter(filtertd11, count11, filtertdMin, filtertdMax);

%twelfth data filter
filtertd12=data(:,3:6,12); %data that we are going to test filter with
%take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the twelfth data set
count12=data(:,1,12); %count of data all along length

[Data12,Time12]=agw_filter(filtertd12, count12, filtertdMin, filtertdMax);

%thirteenth data filter
filtertd13=data(:,3:6,13); %data that we are going to test filter with
%take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the thirteenth data set
count13=data(:,1,13); %count of data all along length

[Data13,Time13]=agw_filter(filtertd13, count13, filtertdMin, filtertdMax);

%fourteenth data filter
filtertd14=data(:,3:6,14); %data that we are going to test filter with
%take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the fourteenth data set
count14=data(:,1,14); %count of data all along length

[Data14,Time14]=agw_filter(filtertd14, count14, filtertdMin, filtertdMax);

%fifteenth data filter
filtertd15=data(:,3:6,15); %data that we are going to test filter with
%take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the fifthteenth data set
count15=data(:,1,15); %count of data all along length

[Data15,Time15]=agw_filter(filtertd15, count15, filtertdMin, filtertdMax);

%sixteenth data filter
filtertd16=data(:,3:6,16); %data that we are going to test filter with
%take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the sixteenth data set
count16=data(:,1,16); %count of data all along length

[Data16,Time16]=agw_filter(filtertd16, count16, filtertdMin, filtertdMax);

%seventeenth data filter
filtertd17=data(:,3:6,17); %data that we are going to test filter with
%take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the seventeenth data set
count17=data(:,1,17); %count of data all along length

[Data17,Time17]=agw_filter(filtertd17, count17, filtertdMin, filtertdMax);

%create an empty cell of 1x17 so that we can put the first Data 1 in the first space, Data 2 in the second space, and so on
c = cell(1,17); 
cc = cell(1,17);
%populate it with all of the data arrays we have
c{1}=Data1;
c{2}=Data2;
c{3}=Data3;
c{4}=Data4;
c{5}=Data5;
c{6}=Data6;
c{7}=Data7;
c{8}=Data8;
c{9}=Data9;
c{10}=Data10;
c{11}=Data11;
c{12}=Data12;
c{13}=Data13;
c{14}=Data14;
c{15}=Data15;
c{16}=Data16;
c{17}=Data17;

%------------Prepare mean, prime, prime^2----------------------------------

Nca=length(c); %length of the cell array, which is 17 because we are using 17 data sets that are put into a cell array

%make the cells readable as a matrix so use cell2mat function
%find the means (ubar)
for i=1:Nca
   tempdata=cell2mat(c(i));
   umean(i)=mean(tempdata(1,:));
   vmean(i)=mean(tempdata(2,:));
   w1mean(i)=mean(tempdata(3,:)); 
   w2mean(i)=mean(tempdata(4,:));
   wavgmean(i)=w2mean(i);
end

% find the perturbations (u')
for i=1:Nca
   tempdata=cell2mat(c(i));
   uprime=(tempdata(1,:)-mean(tempdata(1,:)));
   vprime=(tempdata(2,:)-mean(tempdata(2,:)));
   w1prime=(tempdata(3,:)-mean(tempdata(3,:)));
   w2prime=(tempdata(4,:)-mean(tempdata(4,:)));
   wavgprime=w2prime;
   cc{i}=[uprime; vprime; w1prime; w2prime; wavgprime];
   clear uprime vprime w1prime w2prime wavgprime tempdata 
end 

% find the perturbations squared (u'^2)
for i=1:Nca
   tempdata=cell2mat(c(i));
   uprime2=(tempdata(1,:)-mean(tempdata(1,:))).^2;
   vprime2=(tempdata(2,:)-mean(tempdata(2,:))).^2;
   w1prime2=(tempdata(3,:)-mean(tempdata(3,:))).^2;
   w2prime2=(tempdata(4,:)-mean(tempdata(4,:))).^2;
   wavgprime2=w2prime2;
   cd{i}=[uprime2; vprime2; w1prime2; w2prime2; wavgprime2];
   clear uprime2 vprime2 w1prime2 w2prime2 wavgprime2 tempdata 
end 

%----------------Graph [mean(z/rhalf)]/U0----------------------------------

meanunorm=(umean./U0); %ubar/U0
meanvnorm=(vmean./U0); %vbar/U0
meanwnorm=(wavgmean./U0); %wbar/U0

figure()
subplot(3,1,1);
plot(ZoR,meanunorm)
title('normalized u mean')
ylabel('velocity')
xlim([0 2])
subplot(3,1,2);
plot(ZoR,meanvnorm)
title('normalized v mean')
ylabel('velocity')
xlim([0 2])
subplot(3,1,3);
plot(ZoR,meanwnorm)
title('normalized w mean')
ylabel('velocity')
xlabel('z/rhalf')
xlim([0 2])

%----------------Graph [sqrt(prime squared mean)(z/rhalf)]/U0----------------------------------
for i=1:Nca
    temp=cell2mat(cc(1));
uprime0=cc{i};
vprime0=cc{i};
wprime0=cc{i};

xfluc(i)=(sqrt(mean((uprime0(1,:).^2))))./U0; %sqrt(mean((uprime)^2))/U0
yfluc(i)=(sqrt(mean((vprime0(2,:).^2))))./U0; %sqrt(mean((vprime)^2))/U0
zfluc(i)=(sqrt(mean((wprime0(3,:).^2))))./U0; %sqrt(mean((wprime)^2))/U0
end

figure()
subplot(3,1,1);
plot(ZoR,xfluc)
title('u turbulent fluctuations')
ylabel('velocity')

subplot(3,1,2);
plot(ZoR,yfluc)
title('v turbulent fluctuations')
ylabel('velocity')

subplot(3,1,3);
plot(ZoR,zfluc)
title('w turbulent fluctuations')
ylabel('velocity')
xlabel('z/rhalf')

%----------------Graph average Re stress [mean(prime*prime)(z/rhalf)]/U0^2----------------------------------

for i=1:Nca
uprime=cc{i};
vprime=cc{i};
wprime=cc{i};

xRe(i)=(mean(uprime(1,:).*vprime(2,:)))./U0; %mean(uprime*vprime)(z/rhalf)/U0^2
yRe(i)=(mean(uprime(1,:).*wprime(3,:)))./U0; %mean(uprime*wprime)(z/rhalf)/U0^2
zRe(i)=(mean(vprime(2,:).*wprime(3,:)))./U0; %mean(vprime*wprime)(z/rhalf)/U0^2
end

figure()
subplot(3,1,1);
plot(ZoR,xRe)
title('u Re stress')
ylabel('velocity')
subplot(3,1,2);
plot(ZoR,yfluc)
title('v Re stress')
ylabel('velocity')
subplot(3,1,3);
plot(ZoR,zfluc)
title('w Re stress')
ylabel('velocity')
xlabel('z/rhalf')


%---------------------------------PDF--------------------------------------

%----------At Centerline----------- 
uprime0=cc{3};
vprime0=cc{3};
wprime0=cc{3};

uprimenorm0=uprime0(1,:)./U0; %uprime/U0
sortupn0=sort(uprimenorm0); %sort the data so that when the pdf is made 
vprimenorm0=vprime0(2,:)./U0; %vprime/U0
sortvpn0=sort(vprimenorm0); %sort the data so that when the pdf is made 
wprimenorm0=wprime0(3,:)./U0; %wprime/U0
sortwpn0=sort(wprimenorm0); %sort the data so that when the pdf is made 

% Creating PDF for uprime %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
updf0=pdf('Normal',sortupn0,mean(sortupn0),std(sortupn0));

% Creating PDF for vprime %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vpdf0=pdf('Normal',sortvpn0,mean(sortvpn0),std(sortvpn0));

% Creating PDF for wprime %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wpdf0=pdf('Normal',sortwpn0,mean(sortwpn0),std(sortwpn0));

figure()
subplot(3,1,1);
plot(sortupn0,updf0) 
title('PDF of uprime at centerline')
ylabel('velocity')
subplot(3,1,2);
plot(sortvpn0,vpdf0) 
title('PDF of vprime at centerline')
ylabel('velocity')
subplot(3,1,3);
plot(sortwpn0,wpdf0) 
title('PDF of wprime at centerline')
ylabel('velocity')
xlabel('z/rhalf')

%----------At Z/r(1/2)=1----------- 
uprime1=cc{10};
vprime1=cc{10};
wprime1=cc{10};

uprimenorm1=uprime1(1,:)./U0; %uprime/U0
sortupn1=sort(uprimenorm1); %sort the data so that when the pdf is made 
vprimenorm1=vprime1(2,:)./U0; %vprime/U0
sortvpn1=sort(vprimenorm1); %sort the data so that when the pdf is made 
wprimenorm1=wprime1(3,:)./U0; %wprime/U0
sortwpn1=sort(wprimenorm1); %sort the data so that when the pdf is made

% Creating PDF for uprime %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
updf1=pdf('Normal',sortupn1,mean(sortupn1),std(sortupn1,1));

% Creating PDF for vprime %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vpdf1=pdf('Normal',sortvpn1,mean(sortvpn1),std(sortvpn1,1));

% Creating PDF for wprime %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wpdf1=pdf('Normal',sortwpn1,mean(sortwpn1),std(sortwpn1,1));

figure()
subplot(3,1,1);
plot(sortupn1,updf1) 
title('PDF of uprime at z/rhalf=1')
ylabel('velocity')
subplot(3,1,2);
plot(sortvpn1,vpdf1) 
title('PDF of vprime at z/rhalf=1')
ylabel('velocity')
subplot(3,1,3);
plot(sortwpn1,wpdf1) 
title('PDF of wprime at z/rhalf=1')
ylabel('velocity')
xlabel('z/rhalf')

%----------At Z/r(1/2)=2----------- 
uprime2=cc{15};
vprime2=cc{15};
wprime2=cc{15};

uprimenorm2=uprime2(1,:)./U0; %uprime/U0
sortupn2=sort(uprimenorm2); %sort the data so that when the pdf is made 
vprimenorm2=vprime2(2,:)./U0; %vprime/U0
sortvpn2=sort(vprimenorm2); %sort the data so that when the pdf is made 
wprimenorm2=wprime2(3,:)./U0; %wprime/U0
sortwpn2=sort(wprimenorm2); %sort the data so that when the pdf is made

% Creating PDF for uprime %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
updf2=pdf('Normal',sortupn2,mean(sortupn2),std(sortupn2,1));

% Creating PDF for vprime %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vpdf2=pdf('Normal',sortvpn2,mean(sortvpn2),std(sortvpn2,1));

% Creating PDF for wprime %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wpdf2=pdf('Normal',sortwpn2,mean(sortwpn2),std(sortwpn2,1));

figure()
subplot(3,1,1);
plot(sortupn2,updf2) 
title('PDF of uprime at z/rhalf=2')
ylabel('velocity')
subplot(3,1,2);
plot(sortvpn2,vpdf2) 
title('PDF of vprime at z/rhalf=2')
ylabel('velocity')
subplot(3,1,3);
plot(sortwpn2,wpdf2) 
title('PDF of wprime at z/rhalf=2')
ylabel('velocity')
xlabel('z/rhalf')

%----------------------------Scatter Plots---------------------------------
%----------At Centerline----------- 
uprime0=cc{3};
vprime0=cc{3};
wprime0=cc{3};

up0=uprime0(1,:)./(U0);
vp0=vprime0(2,:)./(U0);
wp0=wprime0(3,:)./(U0);

figure()
subplot(3,1,1);
scatter(up0,vp0) 
title('scatter plot of (uprime*vprime)/(U0^2) at centerline')
ylabel('velocity')
subplot(3,1,2);
scatter(vp0,wp0) 
title('scatter plot of (vprime*wprime)/(U0^2) at centerline')
ylabel('velocity')
subplot(3,1,3);
scatter(up0,wp0) 
title('scatter plot of (uprime*wprime)/(U0^2) at centerline')
ylabel('velocity')
xlabel('z/rhalf')

%----------At Z/r(1/2)=1----------- 
uprime1=cc{10};
vprime1=cc{10};
wprime1=cc{10};

up1=(uprime1(1,:))./((U0));
vp1=(vprime1(2,:))./((U0));
wp1=(wprime1(3,:))./((U0));

figure()
subplot(3,1,1);
scatter(up1,vp1) 
title('scatter plot of (uprime*vprime)/(U0^2) at z/rhalf=1')
ylabel('velocity')
subplot(3,1,2);
scatter(vp1,wp1) 
title('scatter plot of (vprime*wprime)/(U0^2) at z/rhalf=1')
ylabel('velocity')
subplot(3,1,3);
scatter(up1,wp1) 
title('scatter plot of (uprime*wprime)/(U0^2) at z/rhalf=1')
ylabel('velocity')
xlabel('z/rhalf')

%----------At Z/r(1/2)=2----------- 
uprime2=cc{15};
vprime2=cc{15};
wprime2=cc{15};

up2=(uprime2(1,:))./((U0));
vp2=(vprime2(2,:))./((U0));
wp2=(wprime2(3,:))./((U0));

figure()
subplot(3,1,1);
scatter(up2,vp2)
title('scatter plot of (uprime*vprime)/(U0^2) at z/rhalf=2')
ylabel('velocity')
subplot(3,1,2);
scatter(vp2,wp2) 
title('scatter plot of (vprime*wprime)/(U0^2) at z/rhalf=2')
ylabel('velocity')
subplot(3,1,3);
scatter(up2,wp2) 
title('scatter plot of (uprime*wprime)/(U0^2) at z/rhalf=2')
ylabel('velocity')
xlabel('z/rhalf')
