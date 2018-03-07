%Homework 2
% bias error arises only from calibration because then every measurement is
% wrong
% the mean is constant - how to actually get at it
% they tell us what the accuracy is to certain level probe geometry  - they
% have their own tolerance to make the claibration tolerance matrix
% (norteck vectrino to find how tight the calibraiton is)

%PIV within a few pixels on each tic, and those are the sum to find the
%total pixels that we're off (spatial accuracy is the worst case error),
%time accuracy is the timing box to control the yag is delta t (time
%accuracy 1 microsecond)

% root sum of the squares any time you have to combine errors 
clear all
close all
%%
% Load in data from Lab 1
% Error reported from +/-0.5% of measured value to +/- 1mm/s accuracy (whichever is greater)
%Reynolds number for turbulent jet
v=0.010034; %cm2/s at 20 degrees Celcius
Djet=0.595; %cm
Vjet=322.25; %cm/s, found by Q=VA of jet
Re=(Vjet*Djet)/v;   

%row, column, dimesion
tic();
datao(:,:,1) = load('30.9556cm.dat');0.595; %cm
datao(:,:,2) = load('31.2604.dat');
datao(:,:,3) = load('31.5652cm.dat'); %looks most like centerline where we started
datao(:,:,4) = load('31.87cm.dat');  
datao(:,:,5) = load('32.1748cm.dat');
datao(:,:,6) = load('32.4796cm.dat');
datao(:,:,7) = load('32.7844.dat');
datao(:,:,8) = load('33.394cm.dat');
datao(:,:,9) = load('34.0036cm.dat');
datao(:,:,10) = load('34.6132cm.dat');
datao(:,:,11) = load('35.2228cm.dat');
datao(:,:,12) = load('35.8324cm.dat');
datao(:,:,13) = load('36.442cm.dat');
datao(:,:,14) = load('37.0516cm.dat'); %has symmetric bumps at end from resonance bouncing back
datao(:,:,15) = load('37.6612cm.dat'); %has symmetric bumps at end from resonance bouncing back
datao(:,:,16) = load('38.2708cm.dat'); %has symmetric bumps at end from resonance bouncing back
datao(:,:,17) = load('38.8804cm.dat');
% heights that data was gathered at 31.5652,31.87,32.1748,32.4796,32.7844,33.394,34.0036,34.6132,35.2228,35.8324,36.442,37.0516,37.6612,38.2708,38.8804];
D=17; %dimension of data is 17 because there are 17 data sets
LD=length(datao(:,:,1)); %vector to make length of data
toc(); %loading the data takes ~7 sec

%find mean of unfiltered data
for i = 1:D
    xmean(i)= mean(datao(:,3,i));
    ymean(i)= mean(datao(:,4,i));
    zmean1(i) = mean(datao(:,5,i));
    zmean2(i) = mean(datao(:,6,i));
    zmean(i) = mean(datao(:,6,i));
    corrx(i) = mean(datao(:,15,i));
end

%find perturbations of unfiltered data

for i = 1:D 
    for j = 1:LD
     xprime(j,1,i) = datao(j,3,i)-xmean(i);
     yprime(j,1,i) = datao(j,4,i)-ymean(i);
     zprime1(j,1,i) = datao(j,5,i) -zmean1(i);
     zprime2(j,1,i) = datao(j,6,i) -zmean2(i);
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

%first data filter 
filtertd1=datao(:,3:6,1); %data that we are going to test filter with
    %take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the first data set
count1=datao(:,1,1); %count of data all along length
filtertdMin=-1; %RANDOM (shouldbe centered around 0)
filtertdMax=1; %RANDOM (largest value in this set is 1.01)

[Data1,Time1]=agw_filter(filtertd1, count1, filtertdMin, filtertdMax);

%second data filter 
filtertd2=datao(:,3:6,2); %data that we are going to test filter with
    %take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the second data set
count2=datao(:,1,2); %count of data all along length

[Data2,Time2]=agw_filter(filtertd2, count2, filtertdMin, filtertdMax);

%third data filter ------ where we think the midpoint is
filtertd3=datao(:,3:6,3); %data that we are going to test filter with
%take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the third data set
count3=datao(:,1,3); %count of data all along length

[Data3,Time3]=agw_filter(filtertd3, count3, filtertdMin, filtertdMax);

%fourth data filter
filtertd4=datao(:,3:6,4); %data that we are going to test filter with
%take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the fourth data set
count4=datao(:,1,4); %count of data all along length

[Data4,Time4]=agw_filter(filtertd4, count4, filtertdMin, filtertdMax);

%fifth data filter
filtertd5=datao(:,3:6,5); %data that we are going to test filter with
%take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the fifth data set
count5=datao(:,1,5); %count of data all along length

[Data5,Time5]=agw_filter(filtertd5, count5, filtertdMin, filtertdMax);

%sixth data filter
filtertd6=datao(:,3:6,6); %data that we are going to test filter with
%take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the sixth data set
count6=datao(:,1,6); %count of data all along length

[Data6,Time6]=agw_filter(filtertd6, count6, filtertdMin, filtertdMax);

%seventh data filter
filtertd7=datao(:,3:6,7); %data that we are going to test filter with
%take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the seventh data set
count7=datao(:,1,7); %count of data all along length

[Data7,Time7]=agw_filter(filtertd7, count7, filtertdMin, filtertdMax);

%eigth data filter
filtertd8=datao(:,3:6,8); %data that we are going to test filter with
%take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the eigth data set
count8=datao(:,1,8); %count of data all along length

[Data8,Time8]=agw_filter(filtertd8, count8, filtertdMin, filtertdMax);

%ninth data filter
filtertd9=datao(:,3:6,9); %data that we are going to test filter with
%take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the ninth data set
count9=datao(:,1,9); %count of data all along length

[Data9,Time9]=agw_filter(filtertd9, count9, filtertdMin, filtertdMax);

%tenth data filter
filtertd10=datao(:,3:6,10); %data that we are going to test filter with
%take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the tenth data set
count10=datao(:,1,10); %count of data all along length

[Data10,Time10]=agw_filter(filtertd10, count10, filtertdMin, filtertdMax);

%eleventh data filter
filtertd11=datao(:,3:6,11); %data that we are going to test filter with
%take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the eleventh data set
count11=datao(:,1,11); %count of data all along length

[Data11,Time11]=agw_filter(filtertd11, count11, filtertdMin, filtertdMax);

%twelfth data filter
filtertd12=datao(:,3:6,12); %data that we are going to test filter with
%take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the twelfth data set
count12=datao(:,1,12); %count of data all along length

[Data12,Time12]=agw_filter(filtertd12, count12, filtertdMin, filtertdMax);

%thirteenth data filter
filtertd13=datao(:,3:6,13); %data that we are going to test filter with
%take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the thirteenth data set
count13=datao(:,1,13); %count of data all along length

[Data13,Time13]=agw_filter(filtertd13, count13, filtertdMin, filtertdMax);

%fourteenth data filter
filtertd14=datao(:,3:6,14); %data that we are going to test filter with
%take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the fourteenth data set
count14=datao(:,1,14); %count of data all along length

[Data14,Time14]=agw_filter(filtertd14, count14, filtertdMin, filtertdMax);

%fifteenth data filter
filtertd15=datao(:,3:6,15); %data that we are going to test filter with
%take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the fifthteenth data set
count15=datao(:,1,15); %count of data all along length

[Data15,Time15]=agw_filter(filtertd15, count15, filtertdMin, filtertdMax);

%sixteenth data filter
filtertd16=datao(:,3:6,16); %data that we are going to test filter with
%take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the sixteenth data set
count16=datao(:,1,16); %count of data all along length

[Data16,Time16]=agw_filter(filtertd16, count16, filtertdMin, filtertdMax);

%seventeenth data filter
filtertd17=datao(:,3:6,17); %data that we are going to test filter with
%take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the seventeenth data set
count17=datao(:,1,17); %count of data all along length

[Data17,Time17]=agw_filter(filtertd17, count17, filtertdMin, filtertdMax);

%create an empty cell of 1x17 so that we can put the first Data 1 in the first space, Data 2 in the second space, and so on
c = cell(1,17); 
cc = cell(1,17);
%populate it with all of the data arrays we have
c{1}=Data1;
c{2}=Data2;
c{3}=Data3; %centerline
c{4}=Data4;
c{5}=Data5;
c{6}=Data6;
c{7}=Data7;
c{8}=Data8;
c{9}=Data9;
c{10}=Data10; %z/r1/2=1
c{11}=Data11;
c{12}=Data12;
c{13}=Data13;
c{14}=Data14;
c{15}=Data15; %z/r1/2=2
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
   uprime2=sqrt((tempdata(1,:)-mean(tempdata(1,:))).^2);
   vprime2=sqrt((tempdata(2,:)-mean(tempdata(2,:))).^2);
   w1prime2=sqrt((tempdata(3,:)-mean(tempdata(3,:))).^2);
   w2prime2=sqrt((tempdata(4,:)-mean(tempdata(4,:))).^2);
   wavgprime2=w2prime2;
   cd{i}=[uprime2; vprime2; w1prime2; w2prime2; wavgprime2];
   clear uprime2 vprime2 w1prime2 w2prime2 wavgprime2 tempdata 
end 

for i=1:Nca
    temp=cell2mat(cc(1));
uprime0=cc{i};
vprime0=cc{i};
wprime0=cc{i};

xfluc(i)=(sqrt(mean((uprime0(1,:).^2)))); %sqrt(mean((uprime)^2))
yfluc(i)=(sqrt(mean((vprime0(2,:).^2)))); %sqrt(mean((vprime)^2))
zfluc(i)=(sqrt(mean((wprime0(3,:).^2)))); %sqrt(mean((wprime)^2))
end

%%
% Load in data from Lab 2
% need to choose only three locations: boundary layer, log-law, free stream

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
 
%Code from Diego to do an AGW filter%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

Uvel=(U./cal)./delt; % new velocity vector for U (m/s)
Vvel=(V./cal)./delt; % new velocity vector for V (m/s)

Uvelmean=squeeze(mean(Uvel)); % mean velocity vector for U (m/s)
Vvelmean=squeeze(mean(Vvel)); % mean velocity vector for V (m/s)

Vvelmeanvector=mean(Vvelmean,1); % mean velocity vector for V as a 1x48 vector (m/s)
Vvelmax=max(mean(Vvelmean,1)); % max velocity for V as a (m/s)
%
[D1 D2 D3]=size(U); 
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

%%

%ADV%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% You need to first calculate the 0.5% range
% all of the velocities are in m/s
umeanuncert=0.05*umean;
    umeanuncertcent=umeanuncert(3); % bigger than 1 mm/s
    umeanuncertzr1=umeanuncert(10); % bigger than 1 mm/s
    umeanuncertzr2=umeanuncert(15); % bigger than 1 mm/s

vmeanuncert=0.05*vmean;
    vmeanuncertcent=vmeanuncert(3); % smaller than 1 mm/s
    vmeanuncertzr1=vmeanuncert(10); % smaller than 1 mm/s
    vmeanuncertzr2=vmeanuncert(15); % smaller than 1 mm/s
    
wmeanuncert=0.05*w2mean;
    wmeanuncertcent=wmeanuncert(3); % smaller than 1 mm/s
    wmeanuncertzr1=wmeanuncert(10); % smaller than 1 mm/s
    wmeanuncertzr2=wmeanuncert(15); % smaller than 1 mm/s
    
uprime2uncert=0.05*xfluc;
    uprime2uncertcent=uprime2uncert(3);
    uprime2uncertzr1=uprime2uncert(10);
    uprime2uncertzr2=uprime2uncert(15);
    
wprime2uncert=0.05*zfluc;
    wprime2uncertcent=wprime2uncert(3);
    wprime2uncertzr1=wprime2uncert(10);
    wprime2uncertzr2=wprime2uncert(15);

% prime thing    
% need to do sqrt sum of squares because it is two uncertainties that are coming together
xfluczfluc=sqrt((xfluc.^2)+(zfluc.^2));

uwprimecent=sqrt((uprime2uncertcent.^2)+(wprime2uncertcent.^2)); % bigger than 1 mm/s
uwprimezr1=sqrt((uprime2uncertzr1.^2)+(wprime2uncertzr1.^2));    % bigger than 1 mm/s
uwprimezr2=sqrt((uprime2uncertzr2.^2)+(wprime2uncertzr2.^2));    % bigger than 1 mm/s

% to calculate the % uncertainty now for the ones that are not 5%
Vcentperc=((0.0022+0.001)/0.0022)*100;
Wcentperc=((0.0148+0.001)/0.0148)*100;

Vzr1perc=((-0.0056+0.001)/-0.0056)*100;
Wzr1perc=((0.013+0.001)/0.013)*100;

Vzr2perc=((0.0123+0.001)/0.0123)*100;
Wzr2perc=((-0.0089+0.001)/-0.0089)*100;

%%
% ADV
% NOW FOR BOOTFUN
xbootcent=Data3(1,:);
ybootcent=Data3(2,:);
zbootcent=Data3(4,:);

xbootzr1=Data10(1,:);
ybootzr1=Data10(2,:);
zbootzr1=Data10(4,:);

xbootzr2=Data15(1,:);
ybootzr2=Data15(2,:);
zbootzr2=Data15(4,:);

xbootfuncent=mean(bootstrap(1000,@mean,xbootcent));
ybootfuncent=mean(bootstrap(1000,@mean,ybootcent));
zbootfuncent=mean(bootstrap(1000,@mean,zbootcent));

xbootfunzr1=mean(bootstrap(1000,@mean,xbootzr1));
ybootfunzr1=mean(bootstrap(1000,@mean,ybootzr1));
zbootfunzr1=mean(bootstrap(1000,@mean,zbootzr1));

xbootfunzr2=mean(bootstrap(1000,@mean,xbootzr2));
ybootfunzr2=mean(bootstrap(1000,@mean,ybootzr2));
zbootfunzr2=mean(bootstrap(1000,@mean,zbootzr2));

clare=cell2mat(cd(3)); %for centerline
clareucent=clare(1,:);
clarevcent=clare(2,:);
UVrandomcent=clareucent.*clarevcent;

clare1=cell2mat(cd(10)); %for z/r=1
clareu1=clare1(1,:);
clarev1=clare1(2,:);
UVrandom1=clareu1.*clarev1;

clare2=cell2mat(cd(15)); %for z/r=2
clareu2=clare2(1,:);
clarev2=clare2(2,:);
UVrandom2=clareu2.*clarev2;

bootfunuvcent=mean(bootstrap(1000,@mean,UVrandomcent));
bootfunuv1=mean(bootstrap(1000,@mean,UVrandom1));
bootfunuv2=mean(bootstrap(1000,@mean,UVrandom2));

%random error percent fun!
Ucentperc=((0.3739+xbootfuncent)/0.3739)*100;
Vcentperc=((0.0022+ybootfunzr1)/0.0022)*100;
Wcentperc=((0.0148+zbootfuncent)/0.0148)*100;
primeUVcentperc=((0.0148+bootfunuvcent)/0.0148)*100;

Uzr1perc=((-0.0056+xbootfunzr1)/-0.0056)*100;
Vzr1perc=((-0.0056+ybootfunzr1)/-0.0056)*100;
Wzr1perc=((0.013+zbootfunzr1)/0.013)*100;
primeUV1perc=((0.0148+bootfunuv1)/0.0148)*100;

Uzr2perc=((-0.0056+xbootfunzr2)/-0.0056)*100;
Vzr2perc=((0.0123+ybootfunzr2)/0.0123)*100;
Wzr2perc=((-0.0089+zbootfunzr2)/-0.0089)*100;
primeUV2perc=((0.0148+bootfunuv2)/0.0148)*100;
%%

%PIV%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Error for the time calibration on time box is 1 microsecond
timeerror=0.001; %milliseconds

% cal = 33; %pixels/mm
% delt = 2.5; %milliseconds
% U and V are the velocity vectors
% maximum PIV error is when displacement is bigger than expected and time
% is shorter than expected

% maximum time error
PIVterr=(timeerror/delt)*100;

%maximum space error
% Spatial error in ruler calibration
% Original calibration gave 33 pixels/mm
% Later calibration attempt gives 33.7 pixels/mm
% [670 255], [688 930]

% percent error space
PIVsperr=((33.7-33)/33)*100; %gives the percent error between the two measurements
% percent error total
PIVtoterr=sqrt((PIVterr^2)+(PIVsperr.^2));
%total error from primes
PIVprimetoterr=sqrt((PIVtoterr^2)+(PIVtoterr^2));

%Uvelocity 
Umeanreport=[mean(Uvelmean(:,1)) mean(Uvelmean(:,10)) mean(Uvelmean(:,23))];

%Vvelocity
Vmeanreport=[mean(Vvelmean(:,1)) mean(Vvelmean(:,10)) mean(Vvelmean(:,23))];

%UVprimes
UVvelprimebed=mean(Uvelprime(:,1))*mean(Vvelprime(:,1));
UVvelprimelog=mean(Uvelprime(:,10))*mean(Vvelprime(:,10));
UVvelprimefree=mean(Uvelprime(:,23))*mean(Vvelprime(:,23));
UVprimereport=[UVvelprimebed UVvelprimelog UVvelprimefree];

%multiply percent error by values
Uerrorreports=Umeanreport.*PIVtoterr;
Verrorreports=Vmeanreport.*PIVtoterr;
UVprimeerrorreports=UVprimereport.*PIVprimetoterr;
%%
%bootstrap that PIV
%need to take 600*36 for the vector that will be bootstrapped
bedU=reshape(squeeze(Uvel(:,:,1)),[1 21600]);
bedV=reshape(squeeze(Vvel(:,:,1)),[1 21600]);
bedUV=reshape(squeeze(UUprime(:,:,1).*VVprime(:,:,1)),[1 21600]); 

logU=reshape(squeeze(Uvel(:,:,10)),[1 21600]);
logV=reshape(squeeze(Vvel(:,:,10)),[1 21600]);
logUV=reshape(squeeze(UUprime(:,:,10).*VVprime(:,:,10)),[1 21600]); 

freeU=reshape(squeeze(Uvel(:,:,15)),[1 21600]);
freeV=reshape(squeeze(Vvel(:,:,15)),[1 21600]);
freeUV=reshape(squeeze(UUprime(:,:,15).*VVprime(:,:,15)),[1 21600]); 

bootlyUbed=mean(bootstrap(1000,@mean,bedU));
bootlyVbed=mean(bootstrap(1000,@mean,bedV));
bootlyUVbed=mean(bootstrap(1000,@mean,bedUV));

bootlyUlog=mean(bootstrap(1000,@mean,logU));
bootlyVlog=mean(bootstrap(1000,@mean,logV));
bootlyUVlog=mean(bootstrap(1000,@mean,logUV));

bootlyUfree=mean(bootstrap(1000,@mean,freeU));
bootlyVfree=mean(bootstrap(1000,@mean,freeV));
bootlyUVfree=mean(bootstrap(1000,@mean,freeUV));

%random error percent fun!
PIVUcentperc=((0.0145+bootlyUbed)/0.0145)*100;
PIVVcentperc=(((3.98*10^-4)+bootlyVbed)/(3.98*10^-4))*100;
PIVprimeUVcentperc=(((1.91*10^-5)+bootlyUVbed)/(1.91*10^-5))*100;

PIVUzr1perc=((0.0863+bootlyUlog)/0.0863)*100;
PIVVzr1perc=(((-6.44*10^-6)+bootlyVlog)/(-6.44*10^-6))*100;
PIVprimeUV1perc=(((5.47*10^-5)+bootlyUVlog)/(5.47*10^-5))*100;

PIVUzr2perc=((0.0978+bootlyUfree)/0.0978)*100;
PIVVzr2perc=(((1.47*10^-4)+bootlyVfree)/(1.47*10^-4))*100;
PIVprimeUV2perc=(((5.23*10^-5)+bootlyUVfree)/(5.23*10^-5))*100;
