% Zoe Maisel
% CEE 4370
% Lab 2

% Calibration: 33 pixels/mm for images
% Delta t = 2.5 milliseconds for 9 pixel displacement
% fs = 5 Hz
% temp = 22 celcius

% PIV algorithm with subroutines from Todd Cowen
% ARRAY ORGANIZATION:
%                      X,J ----->
%          XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%          X         |                     X
%      -Y, X         | NyStart             X
%       I  X---------xxxxxxxxxxxx          X
%       |  X NxStart x          x          X  Nrows
%       |  X         x          x NySize   X
%       V  X         x          x          X
%          X         xxxxxxxxxxxx          X
%          X           NxSize              X
%          XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%                     Ncols
%       Note:  Positve vectors produced for direction shown
%       Note:  Y is positive up and I is positive down!
%
% Modified March 16, 2013 to process single TIFF images from the Phantom
% Camera for the CEE6370 Boundary Layer Experiment Lab
%
% Flow in this lab is assumed Right to Left in the images with the flume
% bed located near the image bottom (around row=1150)

% Note this version currently set up for file names wwith a fixed prefix
% followed by three digigits where a zero prececedes numbers less than 100
% and two zeros preced numbers less than 10 (i.e., 034, 008).  If you have
% more or less numerical characters there are three lcoations below that
% must be changed and three inavgMin.m - all are marked with comments.

% Note this version currently set up with camera rotated 90 degrees such
% that the wall normal coordinate is x and the wall parallel coordiante is
% y.  The V velocity (y direciton) will be negative given the rotation for
% flows in the down-flume direction.  You will want to multiply those
% velocities by -1.  The u velocity (in the wall normal directiub) will be
% positive upwards so you are all set in that coordinate.

% The camera used produces an image size of 1200 x 1632

% Last updated 03/15/17 by EAC

clear all
close all

% Parameters for processing
Nx=32; % Size of interrogation subwindow in x direction
Ny=32; % Size of interrogation subwindow in y direction
iXshift=0;  % Initial displacement estimate in X direction
iYshift=9;  % Initial displacement estimate in Y direction - from delta t
cal = 33; %pixels/mm
delt = 2.5; %milliseconds

% Files to process
iFileStart=1; % First image to analyze in sequence
iFileStop=1200; % Last image to analyze in sequence
iPairSep=2; % Location of second image in pair with respect to first (e.g., 2= image pair, 3 = image triple)
sPrefix=''; % Constant prefix of image file names
sSuffix='.tif'; % Constant suffix of image file names
sDataSuffix='.dat'; % Constant suffix of ASCII data file names
sImagePath='/Users/Zoeannem/Documents/MATLAB/Fluids_Lab_2/zeb/'; % path to image directory
sOutputPath='/Users/Zoeannem/Documents/MATLAB/Fluids_Lab_2/zeb/'; % path to results directory - NOTE THIS DIRECTORY MUST EXIST

% Region to process
NxStart=27; % left edge of start of interrogation region
NyStart=25; % upper edge of start of interrogation region
NxSize=1536; % horizontal length of interrogation region
NySize=1152; % vertical length of interrogation region

%Overlap processing (/1 = no overlap, /2=50% overlap, /4=75% overlap)
NxStep=Nx/1;  % horizontal distance between subwindows
NyStep=Ny/1;  % vertical distance between subwindows

%Parameters for figures
Cmin=4;     % Min color level to show in pseudocolored images
Cmax=50;    % Max color level to show in pseudocolored images

% Check to see if NxSize <= Nx and NySize <= Ny
sFile1=num2str(iFileStart,'%04d'); %note this sets the numeric value in filename
%if your filename has a number of nueric values different then 3 change
%this to reflect that numbher of numeric values.  You also must cange it in
%in AvgMin
fname=[sImagePath sPrefix sFile1 sSuffix];
disp(fname)
I1=double(imread(fname)); 
[Nrows Ncols]=size(I1)
if NySize+NyStart-1 > Nrows
    disp(' NySize+NyStart-1 must be <= Nrows')
    return
elseif NxSize+NxStart-1 > Ncols
    disp(' NxSize+NxStart-1 must be <= Ncols')
    return
end           

% Determine number of subwindows within interrogation region

NsubX=floor((NxSize-Nx)/NxStep)+1;
NsubY=floor((NySize-Ny)/NyStep)+1;

% Determine interrogation grid nodes
IgX=zeros(NsubY,NsubX); %Initialize arrays as this leads to faster execution in Matlab
IgY=IgX;
IgXc=IgX;
IgYc=IgX;
for i=1:NsubY
   for j=1:NsubX
     IgX(i,j)= NxStart+(j-1)*NxStep-1;  % one left of subwindow start
     IgY(i,j)= NyStart+(i-1)*NyStep-1;  % one above subwindow start
     IgXc(i,j)= NxStart+(j-1)*NxStep+Nx/2-0.5;   % centers of subwindows
     IgYc(i,j)= NyStart+(i-1)*NyStep+Ny/2-0.5;   % centers of subwindows
   end
end

% Loop over image pairs

% Determine Pre-process images
disp(sprintf('Pre-processing images ...\n'))
[Bkgrnd1,Bkgrnd2,iMean]=avgMin(iFileStart,iPairSep,iFileStop,sImagePath,sPrefix,sSuffix);
save([sOutputPath 'backgroundImg1'], 'Bkgrnd1') %save background images to .mat files
save([sOutputPath 'backgroundImg2'], 'Bkgrnd2') %save background images to .mat files

% Plot images Bkgrnd1 and Bkgrnd2 and curves iMean
FigAvgMin(Bkgrnd1,Bkgrnd2,iMean)
   
% Plot mean intensity along rows to find bed
figure
plot(mean(Bkgrnd1))
xlabel('Row (pixel)')
ylabel('Minimum intensity (counts)')
disp(sprintf('\n'))
sKeyboard = input('Hit return to being image pair analysis.');

for iFile=iFileStart:iPairSep:iFileStop

% Read in Images into arrays I1 and I2
    sFile1=num2str(iFile,'%04d');  %must change this if number of numeric 
%values in filename is differnt than 3
    fname=[sImagePath sPrefix sFile1 sSuffix];
    I1=double(imread(fname)); %loads image into a double precision real variable matrix  
  
    sFile1=num2str(iFile+1,'%04d');  %must change this if number of numeric 
%values in filename is differnt than 3
    fname=[sImagePath sPrefix sFile1 sSuffix];
    I2=double(imread(fname)); %loads image into a double precision real variable matrix  
    disp(['Reading image number: ' sFile1])
   
% Subtract pre-process images from original images
   I1=I1-Bkgrnd1;
   I2=I2-Bkgrnd2;
 
% Plot images I1 and I2 and difference (I2-I1)
   FigImagesAndDifferences(I1,I2, Cmin, Cmax)
     
% Loop over subwindows
   t0=clock;
   for nsY=1:NsubY
      for nsX=1:NsubX     
% Load subwindows into F and G
         for i=1:Ny
            for j=1:Nx
               F(i,j)=I1(IgY(nsY,nsX)+i,IgX(nsY,nsX)+j);
               G(i,j)=I2(IgY(nsY,nsX)+i+iYshift,IgX(nsY,nsX)+j+iXshift);
            end
         end
% Fast PIV - periodic
         H=ifft2(fft2(F).*conj(fft2(G)));
         H(1:Ny,2:Nx)=fliplr(H(1:Ny,2:Nx)); %correct for matlab oddly stored 2-D fft step-1
         H(2:Ny,1:Nx)=flipud(H(2:Ny,1:Nx)); %correct for matlab oddly stored 2-D fft step-2
      
% Call integer diplacement finding subroutine and return vector of
% magnitudes at nearest neighbor locations
         [a,nXdis(nsY,nsX),nYdis(nsY,nsX)] = Pixel_displacement(H,Nx,Ny);

% Call sub-pixel accuracy subroutine
         [Xdis(nsY,nsX),Ydis(nsY,nsX)]=Subpixel_displacement(a);
      end
   end
   t_elapse=etime(clock,t0);
   disp(sprintf('Elapsed time for one image pair: %6.2f\n',t_elapse))
   
% Set Xdis to sum of integer and sub-pixel parts
   Xdis=nXdis+Xdis+iXshift;
   Ydis=nYdis-Ydis-iYshift;
   
% Plot figure showing displacements
   FigDisplacements(Xdis,Ydis,Nx,Ny)
   pause(0.01)
   
% Output Raw Vectors - note directory in sOutputPath MUST exist prior to running this code   
   %outputRawVectorsASCII(NsubY,IgXc,IgYc,Xdis,Ydis,sOutputPath,sFrame,sDataSuffix); %4 ASCII files  
   outputRawVectorsMat(IgXc,IgYc,Xdis,Ydis,sOutputPath,sFile1); %1 Mat file
end

%Read in processed data and arrange in 3-D array of U and V velocities
binCount=0;
for iFile=iFileStart+iPairSep-1:iPairSep:iFileStop
   binCount=binCount+1; %track number of data points at each gridpoint (bin)
   sFrame=sprintf('%04d',iFile); %note this sets the numeric value in filename
%if your filename has a number of nueric values different then 3 change
   sFilename=[sOutputPath 'xyM_' sFrame];
   eval(['load ' sFilename ])
% Store data in 3-D matrix with time running in the first index
   U(binCount,:,:)=Xdis;
   V(binCount,:,:)=Ydis;
end

%save results to UV.mat
save([sOutputPath 'UV'], 'U','V') %save background images to .mat files

%% Code from Diego to do an AGW filter%%%%%%%%%%%%%%%%%%
dims=size(U);
%%% Filter U
Vel=zeros(size(U,1)*size(U,2),2,size(U,3));
VelBack=Vel;

for i=1:size(U,3);
   Vel(:,1,i)=reshape(U(:,:,i),[],1);
   %Vel(:,2,i)=reshape(U(:,:,i),[],1);
end

t=[0:1:length(Vel)-1];

for i=1:size(U,3);
  [uT, tf] = agw_filter(Vel(:,:,i),t,-inf,inf);
  data=interp1(tf',uT',t','linear','extrap');
  VelBack(:,:,i)=data;
end
%%% Filter V
VelV=zeros(size(V,1)*size(V,2),2,size(V,3));
VelBackV=VelV;

for i=1:size(V,3);
   VelV(:,1,i)=reshape(V(:,:,i),[],1);
   %VelV(:,2,i)=reshape(V(:,:,i),[],1);
end

t=[0:1:length(VelV)-1];

for i=1:size(V,3);
  [vT, tf] = agw_filter(VelV(:,:,i),t,-inf,inf);
  data=interp1(tf',vT',t','linear', 'extrap');
  VelBackV(:,:,i)=data;
end

U=reshape(VelBack(:,1,:),dims); %reshape vectors so that they are the correct size by U
V=-1.*(reshape(VelBackV(:,1,:),dims));

%%
Ubar=squeeze(mean(U,1));  % find the temporal mean of the results and store in 2D array
Vbar=squeeze(mean(V,1)); % pixel displacement (mm)

Umed=squeeze(median(U,1));  % find the temporal median of the results and store in 2D array
Vmed=squeeze(median(V,1)); % pixel displacement (mm)

Udisp=(U./cal); % new displacement vector for U (mm)
Vdisp=(V./cal); % new displacement vector for V (mm)

Uvel=(U./cal)./delt; % new velocity vector for U (m/s)
Vvel=(V./cal)./delt; % new velocity vector for V (m/s)

Uvelmean=squeeze(mean(Uvel)); % mean velocity vector for U (m/s)
Vvelmean=squeeze(mean(Vvel)); % mean velocity vector for V (m/s)

Vvelmeanvector=mean(Vvelmean,1); % mean velocity vector for V as a 1x48 vector (m/s)
Vvelmax=max(mean(Vvelmean,1)); % max velocity for V as a (m/s)

[D1 D2 D3]=size(U); %D1 is the number of image pairs, D2 is the NsubY, D3 is the NsubX

%% Plots

%plot mean vector field 
figure
quiver(IgXc,IgYc,Ubar,Vbar) % multiply velocity by negative 1 to correct for negative error
axis ij
axis image

%plot median velocity profile assuming horizontally homogenous
% - have not removed spurrious vectors yet so use median
figure
z=1:NsubX;
plot(median(Ubar),z,'b',median(Vbar),z,'r') % multiply velocity by negative 1 to correct for negative error
legend('U median','V median')
title('Median of horizontally and temporally averaged data')
xlabel('U, V (mm)')
ylabel('Interrogation grid point')

Uprime=U-repmat(mean(U,1),[D1 1 1]);  %find perturbations based on Median - uprime
Vprime=V-repmat(mean(V,1),[D1 1 1]);  %find perturbations based on Median - vprime

UUprime=Uprime.^2;
VVprime = Vprime.^2;

%keeping time
Uvelprimetime = squeeze(((Uprime))./cal./delt); % working vector of velocity primes, squeezed into simple 36x48 matrix
Vvelprimetime = squeeze(((Vprime))./cal./delt);

%getting rid of time
Uvelprime = squeeze(sqrt(mean(UUprime))./cal./delt); % working vector of velocity primes, squeezed into simple 36x48 matrix
Vvelprime = squeeze(sqrt(mean(VVprime))./cal./delt);

uu=sort(Uprime,1); %sort the perturbations at each i,j (i.i., in time)
vv=sort(Vprime,1); %sort the perturbations at each i,j (i.i., in time)

%normal distribution +/-1 standard deviation is at about 16% and 84%

ind84=round(D1*0.84);

%determine the robust estimate of the variance of u and w based on the 84% index
upIQR=squeeze(uu(ind84,:,:));
wpIQR=squeeze(vv(ind84,:,:));

%% Velocity Field
% plot IQR determined estimate of RMS velocity fluctuation profile 
%assuming horizontally homogenous.
% - have not removed spurrious vectors yet so use IQR - a robust estimator of the RMS)
figure
plot(median(upIQR)/cal/delt,z,'b',median(wpIQR)/cal/delt,z,'r') %divide by calibration constant to convert from pixels to mm, divide by delt (m/s)
legend('(u''^2)^{0.5}','(v''^2)^{0.5}')
xlabel('(u''^2)^{0.5}, (v''^2)^{0.5} (m/\s)')
ylabel('Interrogation grid point')
title('RMS velocity fluctions of horizontally and temporally averaged data')

%% Determining log layer
loglayer=10;

%% Histograms 

%of x and z displacements
figure()
temp=squeeze(Udisp(:,:,1));
subplot(3,1,1); hist(temp(:),100); title('histogram of displacement normal to wall at bed') % histogram of displacement at bed
temp=squeeze(Udisp(:,:,loglayer));
subplot(3,1,2); hist(temp(:),100); title('histogram of displacement normal to wall at log-layer') %histogram of displacement at log-layer
temp=squeeze(Udisp(:,:,23));
subplot(3,1,3); hist(temp(:),100); title('histogram of displacement normal to wall at free stream') %histogram of displacement at free stream

figure()
temp=squeeze(Vdisp(:,:,1));
subplot(3,1,1); hist(temp(:),100); title('histogram of displacement parallel to wall at bed')% histogram of displacement at bed
temp=squeeze(Vdisp(:,:,loglayer));
subplot(3,1,2); hist(temp(:),100); title('histogram of displacement parallel to wall at log-layer') %histogram of displacement at log-layer
temp=squeeze(Vdisp(:,:,23));
subplot(3,1,3); hist(temp(:),100); title('histogram of displacement parallel to wall at free stream')%histogram of displacement at free stream

% of x' and w' 
figure()
temp=squeeze(Uprime(:,:,1));
subplot(3,1,1); hist(temp(:),100); title('histogram of perturbations normal to wall at bed') % histogram of wprime at bed
temp=squeeze(Uprime(:,:,loglayer));
subplot(3,1,2); hist(temp(:),100); title('histogram of perturbations normal to wall at log-layer') %histogram of wprime at log-layer
temp=squeeze(Uprime(:,:,23));
subplot(3,1,3); hist(temp(:),100); title('histogram of perturbations normal to wall at free stream') %histogram of wprime at free stream

figure()
temp=squeeze(Vprime(:,:,1));
subplot(3,1,1); hist(temp(:),100); title('histogram of perturbations parallel to wall at bed')% histogram of uprime at bed
temp=squeeze(Vprime(:,:,loglayer));
subplot(3,1,2); hist(temp(:),100); title('histogram of perturbations parallel to wall at log-layer') %histogram of uprime at log-layer
temp=squeeze(Vprime(:,:,23));
subplot(3,1,3); hist(temp(:),100); title('histogram of perturbations parallel to wall at free stream')%histogram of uprime at free stream

%Display color gradient of flow 
figure()
imagesc(Ubar')
colorbar
title('flow normal to wall')
Vdisp = imrotate(Vbar,90);
figure()
imagesc(Vdisp)
colorbar 
title('flow parallel to wall')

%% Scatter Plots 
% Uprime and Vprime are the U' and V' vectors

Ubedprime=reshape((Uvelprimetime(:,:,1)),[21600 1]); % Uprime velocity (m/s)
Ulogprime=reshape((Uvelprimetime(:,:,loglayer)),[21600 1]);
Ufstrmprime=reshape((Uvelprimetime(:,:,23)),[21600 1]);

Vbedprime=reshape((Vvelprimetime(:,:,1)),[21600 1]);
Vlogprime=reshape((Vvelprimetime(:,:,loglayer)),[21600 1]);
Vfstrmprime=reshape((Vvelprimetime(:,:,23)),[21600 1]);

figure()
subplot(3,1,1);
scatter(Vbedprime,Ubedprime) 
title('scatter plot of uprime*vprime at near-bed')
xlabel('u'' (m/s)')
ylabel('w'' (m/s)')
axis equal

subplot(3,1,2);
scatter(Vlogprime,Ulogprime) 
title('scatter plot of uprime*vprime at log-layer')
xlabel('u'' (m/s)')
ylabel('w'' (m/s)')
axis equal

subplot(3,1,3);
scatter(Vfstrmprime,Ufstrmprime) 
title('scatter plot of uprime*vprime at free stream')
xlabel('u'' (m/s)')
ylabel('w'' (m/s)')
axis equal


%% Determining v*
 
nu = 0.9565*10^-6; %visosity for temperature at 22C, m^2/s
delx = 0.0009697; %m
znorm=z*delx;

% figure();
% scatter(log(znorm),mean(Vbar))
% % make a linear fit    
% hold on    
% lsline;
% slope values: point 1 [x=3.403, Y=8.363]
              % point 2 [x=3.06,  Y=7.982]
slope = (8.363-7.982)/(3.403-3.06)/100; %found from graph hardcoded in

Kappa = 0.41; % von Karman's constant
vstar = slope*Kappa; 
              
%% Determining Vplus(zplus)
% vplus = vbar/vstar
Vplus = mean(Vvelmean./vstar); % equation from notes, want to compress Vplu so it is along only one dimension
% should be about 5% of maximum velocity 

% zplus = (z)*(v*)/nu
zplus = (znorm.*vstar)./nu; %multiply by height of each subwindow

figure();
scatter(log(znorm),mean(Vvelmean))
title('U vs lnz')
xlabel('lnz (z in m)')
ylabel('U (m/s)')
% make a linear fit    
figure()
vstarest = 0.052*Vvelmax;
zplusest = (znorm.*vstar)./nu;
uplusest = Vvelmeanvector/vstarest; 
loglaw = 1/Kappa*log(zplusest)+5.5;
semilogx(zplusest,uplusest,'o',zplusest,loglaw)
xlabel('z+')
ylabel('U+')
title('U+ vs z+')


 %calculating vstar using two different ways shows that it is 0.0058
 %between
 %0.0046
 
% subplot(2,1,2)
% semilogx(zplus,Vplus) %divide by calibration constant to convert from pixels to mm, divide by delt (m/s)

%
%% Determining v'plus(zplus), u'plus(zplus)
Uprimeplus = mean(Uvelprime./vstar); 
Vprimeplus = mean(Vvelprime./vstar); 

figure
subplot(2,1,1); plot(zplus,Vprimeplus); title('Uprime+ vs z+'); ylabel('Uprime+') 
subplot(2,1,2); plot(zplus,Uprimeplus); title('Wprime+ vs z+'); ylabel('Wprime+') 
xlabel('z+')

% 
%% Spectral Data
% Need to plot spectral data at each z+ (multiple elevations on a single
% plot)

% ks = 2*pi/delx

fs=5; %Hz from lab
delx = 0.9697; %mm
ks = 2*pi/delx; % convert f to temporal
N=length(zplus); % N is the length of the vel vector
Len=36; %dictates the length of each bin
NN=N/Len; %number of bins that are created by the length L
 
Saa = zeros(25,36,48);
for r = 1:(iFileStop/2)
% autospectral density function for the first data row
for i=1:N                  %index from 1 to the number of bin being created
%     ParvelZbed = (V(r,:,i));
    atemp=(V(r,:,i))-mean(Vvelmean(i),1); %for index 1, the bin will go from 1 to L. For index 2, the bin will go from L+1 to 2L 
    A=abs(fft(atemp)); %do fast fourier transform
    AA=A.*conj(A); %calculate the conjugate, has to be .* to make sure its element by element
    Saatemp=(AA)*(1/(Len*ks)); %gives discrete domains, and creates a vector for spectral peaks of data (equation given in notes)
    Saa(r,:,i) = (Saatemp);
end

end 

kMaxRes=ks/Len; %max resolution frequency 
k=linspace(0,ks-kMaxRes,Len); %makes xaxis 

Saa = squeeze(mean(Saa)); %make a 36x48 matrix 
% autospectral density function for the first data row

figure()

for i=1:N  
    loglog(k,Saa(:,i),'color',rand(1,3))
    hold on
end
hold on
loglog(k,k.^(-5/3));
ylabel('Saa')
xlabel('k (rad/mm)')
title('Saa(k) at all z elevations')
hold off 

figure()
loglog(k,Saa(:,1),'color',rand(1,3))
hold on
loglog(k,Saa(:,loglayer),'color',rand(1,3))
hold on
loglog(k,Saa(:,23),'color',rand(1,3))
hold on
loglog(k,k.^(-5/3));
ylabel('Saa')
xlabel('k (rad/mm)')
title('Saa(k) at three elevations')
legend('bed','log-layer','free stream','-5/3 slope')
hold off

%% Dissipation Energy WRONG
% compensated spectra %p.s. education has all the vowels, like unkangaroolike(y)

kaxis = (k.^(5/3));

figure()
Saa1 = zeros(36,48);
for i=1:N  
    for j=1:36
        Saa1(j,i) = Saa(j,i).*kaxis(j);
        loglog(k,Saa1(:,i),'color',rand(1,3))
    hold on
    end
end
hold on
ylabel('Saa')
xlabel('k (rad/mm)')
xlim([0 10])
hold off 

figure()
loglog(k,Saa1(:,1),'color',rand(1,3))
hold on
loglog(k,Saa1(:,loglayer),'color',rand(1,3))
hold on
loglog(k,Saa1(:,23),'color',rand(1,3))
ylabel('Saa')
xlabel('k (rad/mm)')
title('Compensated Saa(k) at three elevations')
legend('bed','log-layer','free stream')
hold off

