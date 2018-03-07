% PIV algorithm with subroutines
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
% velocities by -1.  The u velocity (in the wall nornak directiub) will be
% positive upwards so you are all set in that coordinate.

% The camera used prodoces and image size of 1200 x 1632

% Last updated 03/15/17 by EAC

clear all
close all

% Parameters for processing
Nx=64; % Size of interrogation subwindow in x direction
Ny=64; % Size of interrogation subwindow in y direction
iXshift=0;  % Initial displacement estimate in X direction
iYshift=12;  % Initial displacement estimate in Y direction

% Files to process
iFileStart=1; % First image to analyze in sequence
iFileStop=10; % Last imgage to analyze in sequence
iPairSep=2; % Location of second image in pair with respect to first (e.g., 2= image pair, 3 = image triple)
sPrefix='Dt3ms'; % Constant prefix of image file names
sSuffix='.tif'; % Constant suffix of image file names
sDataSuffix='.dat'; % Constant suffix of ASCII data file names
sImagePath='/Users/eac/classes/Cee637/PIVLab/Lab02_2017_Diego/'; % path to image directory
sOutputPath='/Users/eac/classes/Cee637/PIVLab/Results/' % path to results directory - NOTE THIS DIRECTORY MUS EXIST

% Region to process
NxStart=6; % left edge of start of interrogation region
NyStart=25; % upper edge of start of interrogation region
NxSize=1536; % horizontal length of interrogation region
NySize=1152; % horizontal length of interrogation region

%Overlap processing (/1 = no overlap, /2=50% overlap, /4=75% overlap)
NxStep=Nx/1;  % horizontal distance between subwindows
NyStep=Ny/1;  % vertical distance between subwindows

%Parameters for figures
Cmin=4;     % Min color level to show in pseudocolored images
Cmax=50;    % Max color level to show in pseudocolored images

% Check to see if NxSize <= Nx and NySize <= Ny
sFile1=num2str(iFileStart,'%03d'); %note this sets the numeric value in filename
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
    sFile1=num2str(iFile,'%03d');  %must change this if number of numeric 
%values in filename is differnt than 3
    fname=[sImagePath sPrefix sFile1 sSuffix];
    I1=double(imread(fname)); %loads image into a double precision real variable matrix  
  
    sFile1=num2str(iFile+1,'%03d');  %must change this if number of numeric 
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
   
% Output Raw Vectors - note director in sOutputPath MUST exist prior to running this code   
   %outputRawVectorsASCII(NsubY,IgXc,IgYc,Xdis,Ydis,sOutputPath,sFrame,sDataSuffix); %4 ASCII files  
   outputRawVectorsMat(IgXc,IgYc,Xdis,Ydis,sOutputPath,sFile1); %1 Mat file
end

%Read in processed data and arrange in 3-D array of U and V velocities
binCount=0;
for iFile=iFileStart+iPairSep-1:iPairSep:iFileStop
   binCount=binCount+1; %track number of data points at each gridpoint (bin)
   sFrame=sprintf('%03d',iFile);; %note this sets the numeric value in filename
%if your filename has a number of nueric values different then 3 change
   sFilename=[sOutputPath 'xyM_' sFrame];
   eval(['load ' sFilename ])
% Store data in 3-D matrix with time running in the first index
   U(binCount,:,:)=Xdis;
   V(binCount,:,:)=Ydis;
end

%save results to UV.mat
save([sOutputPath 'UV'], 'U','V') %save background images to .mat files

Ubar=squeeze(mean(U,1));  % find the temporal mean of the results and store in 2D array
Vbar=squeeze(mean(V,1));

Umed=squeeze(median(U,1));  % find the temporal median of the results and store in 2D array
Vmed=squeeze(median(V,1));

[D1 D2 D3]=size(U); %D1 is the number of image pairs, D2 is the NsubY, D3 is the NsubX

%plot median vector field 
% - have not removed spurrious vectors yet so use median)
figure
quiver(IgXc,IgYc,Umed,Vmed)
axis ij
axis image

%plot median velocity profile assuming horizontally homogenous
% - have not removed spurrious vectors yet so use median)
figure
z=1:NsubX;
plot(median(Umed),z,'b',median(Vmed),z,'r')
legend('U median','V median')
title('Median of horizontally and temporally averaged data')
xlabel('U, V (pixels/\Delta t)')
ylabel('Interrogation grid point')

UUmed=U-repmat(median(U,1),[D1 1 1]);  %find perturbations based on Median
VVmed=V-repmat(median(V,1),[D1 1 1]);  %find perturbations based on Median

uu=sort(UUmed,1); %sort the perturbations at each i,j (i.i., in time)
vv=sort(VVmed,1); %sort the perturbations at each i,j (i.i., in time)

%normal distribution +/-1 standard deviation is at about 16% and 84%

ind84=round(D1*0.84);

%determine the robust estimate of the variance of u and w based on the 84% index
upIQR=squeeze(uu(ind84,:,:));
wpIQR=squeeze(vv(ind84,:,:));

%plot IQR determined estimate of RMS velocity fluctuation profile 
%assuming horizontally homogenous.
% - have not removed spurrious vectors yet so use IQR - a robust estimator of the RMS)
figure
plot(median(upIQR),z,'b',median(wpIQR),z,'r')
legend('(u''^2)^{0.5}','(v''^2)^{0.5}')
xlabel('(u''^2)^{0.5}, (v''^2)^{0.5} (pixels/\Delta t)')
ylabel('Interrogation grid point')
title('RMS velocity fluctions of horizontally and temporally averaged data')


