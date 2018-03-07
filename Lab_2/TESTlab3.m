%Lab 3

% Notes on LIF dataset 2017
% 
% Images were collected at 5 Hz. We used a 50 mm lens with an aperture of
% f1.4. A Number 16 filter was placed on the lens to remove particles from
% the image.
% 
% Images follow the naming convention ?LIF_####.tif" and are
% 16 bit uncompressed TIFF images coming from a 14bit CMOS sensor.
% 
% To read the images into Matlab's workspace, try the following commands
% from the image directory:
% 
% image1 = double( imread( ?LIF_####.tif' ) ); imagesc( image1 ) set(gca,
% 'clim', [0 255] )
% 
% If you have questions or problems reading the image, please email Prof.
% Cowen (eac20) or the TA Diego Muriel (dfm99).
% 
% Calibration data is provided within the same folder. The
% calibration levels start with a blank (i.e. 0 concentration flourecein)
% and are incremented at 10 ppb. The number after calibration in the
% filename specifies how many increments above the blank have been added 
% as follows: calib_00_##.tif 0ppb; calib_01_##.tif 10ppb; calib_02_##.tif 
% 20ppb; calib_03_##.tif 30ppb; calib_04_##.tif 40ppb; calib_05_##.tif 
% 50ppb; calib_06_##.tif 70ppb; calib_07_##.tif 100ppb; calib_08_##.tif 
% 140ppb. calib_newBlank_##.tif is also provided if you want to account 
% for the small amount of tracer spilled when removing the fish tank. 
% 
% For spatial calibration use calib_ruler.tif (ruler in vertical position),
%  and for the location of the jet origin use calib_position_##.tif and 
% calib_xd_50.tif. The vertical string to the left of calib_xd_50.tif 
% corresponds to x/d=50. Use this to confirm your estimate of the origin.
% 
% 
% The source concentration for the jet was ~400 ppb. You will need to verify 
% this.

clear all
close all

concentration = [0 10 20 30 40 50 70 100 140];

iFileStart=1;
iFileEnd=300;
sPath='/Users/Zoeannem/Documents/MATLAB/Fluids_Lab_3/';
sPrefix='LIF_';
sSuffix='.tif';

Nrows=1200;
Ncols=1632;
Sum=zeros(Nrows,Ncols);

N=0;
for i=iFileStart:iFileEnd
   N=N+1;
   sName=num2str(i,'%04d');
   fname=[sPath sPrefix sName sSuffix]
   img=double(imread(fname));
   Sum=Sum+img;
   figure(1)
   imagesc(img,[10,30000])
   title(sName)
   pause(0.1)
end
Cbar=Sum/N;

figure(1)
set(gca,'fontsize',12)
imagesc(Cbar)
colorbar
axis image
xlabel('x (pixels)')
ylabel('y (pixels)')
%%
% load in the ruler and x/d images
sPath = '/Users/Zoeannem/Documents/MATLAB/Fluids_Lab_3/';
sPrefix='calib_';
sSuffix='.tif';
sName='ruler';

ruler=[sPath sPrefix sName sSuffix]
imgruler=double(imread(ruler));

figure()
imagesc(imgruler)

sPath = '/Users/Zoeannem/Documents/MATLAB/Fluids_Lab_3/';
sPrefix='calib_';
sSuffix='.tif';
sName='xd_50';

ruler=[sPath sPrefix sName sSuffix]
imgxd=double(imread(ruler));

figure()
imagesc(imgxd)

%calibration data:
% (901,475), (901, 521)
cal = (521-475); %pixels/cm




%%
% Do calibration tests
xaxis = 1632;
yaxis = 1200;
sPath = '/Users/Zoeannem/Documents/MATLAB/Fluids_Lab_3/';
sPrefix='calib_newBlank_';
sSuffix='.tif';
iImageStart = 1;
iImageEnd= 10;

blank = zeros(10,1200,1632);
for i=iImageStart:iImageEnd
   sName=num2str(i,'%02d');
   fname=[sPath sPrefix sName sSuffix]
   img=double(imread(fname));
   title(sName)
   blank(i,:,:) = img;
end
medmat(1,:,:) = reshape(median(blank),[1200 1632]);

%%%%%%%%%%%
sPath = '/Users/Zoeannem/Documents/MATLAB/Fluids_Lab_3/';
sPrefix='calib_01_';
sSuffix='.tif';
iImageStart = 1;
iImageEnd= 10;

one = zeros(10,1200,1632);
for i=iImageStart:iImageEnd
   sName=num2str(i,'%02d');
   fname=[sPath sPrefix sName sSuffix]
   img=double(imread(fname));
   title(sName)
   one(i,:,:) = img;
end
medmat(2,:,:) = reshape(median(one),[1200 1632]);
%%%%%%%%%%%%%

sPath = '/Users/Zoeannem/Documents/MATLAB/Fluids_Lab_3/';
sPrefix='calib_02_';
sSuffix='.tif';
iImageStart = 1;
iImageEnd= 10;

two = zeros(10,1200,1632);
for i=iImageStart:iImageEnd
   sName=num2str(i,'%02d');
   fname=[sPath sPrefix sName sSuffix]
   img=double(imread(fname));
   title(sName)
   two(i,:,:) = img;
end
medmat(3,:,:) = reshape(median(two),[1200 1632]);
%%%%%%%%%%%%

sPath = '/Users/Zoeannem/Documents/MATLAB/Fluids_Lab_3/';
sPrefix='calib_03_';
sSuffix='.tif';
iImageStart = 1;
iImageEnd= 10;

three = zeros(10,1200,1632);
for i=iImageStart:iImageEnd
   sName=num2str(i,'%02d');
   fname=[sPath sPrefix sName sSuffix]
   img=double(imread(fname));
   title(sName)
   three(i,:,:) = img;
end
medmat(4,:,:) = reshape(median(three),[1200 1632]);
%%%%%%%%%%%%

sPath = '/Users/Zoeannem/Documents/MATLAB/Fluids_Lab_3/';
sPrefix='calib_04_';
sSuffix='.tif';
iImageStart = 1;
iImageEnd= 10;

four = zeros(10,1200,1632);
for i=iImageStart:iImageEnd
   sName=num2str(i,'%02d');
   fname=[sPath sPrefix sName sSuffix]
   img=double(imread(fname));
   title(sName)
   four(i,:,:) = img;
end
medmat(5,:,:) = reshape(median(four),[1200 1632]);  
%%%%%%%%%%%%

sPath = '/Users/Zoeannem/Documents/MATLAB/Fluids_Lab_3/';
sPrefix='calib_05_';
sSuffix='.tif';
iImageStart = 1;
iImageEnd= 10;

five = zeros(10,1200,1632);
for i=iImageStart:iImageEnd
   sName=num2str(i,'%02d');
   fname=[sPath sPrefix sName sSuffix]
   img=double(imread(fname));
   title(sName)
   five(i,:,:) = img;
end
medmat(6,:,:) = reshape(median(five),[1200 1632]);
%%%%%%%%%%%%

sPath = '/Users/Zoeannem/Documents/MATLAB/Fluids_Lab_3/';
sPrefix='calib_06_';
sSuffix='.tif';
iImageStart = 1;
iImageEnd= 10;

six = zeros(10,1200,1632);
for i=iImageStart:iImageEnd
   sName=num2str(i,'%02d');
   fname=[sPath sPrefix sName sSuffix]
   img=double(imread(fname));
   title(sName)
   six(i,:,:) = img;
end
medmat(7,:,:) = reshape(median(six),[1200 1632]);
%%%%%%%%%%%%

sPath = '/Users/Zoeannem/Documents/MATLAB/Fluids_Lab_3/';
sPrefix='calib_07_';
sSuffix='.tif';
iImageStart = 1;
iImageEnd= 10;

seven = zeros(10,1200,1632);
for i=iImageStart:iImageEnd
   sName=num2str(i,'%02d');
   fname=[sPath sPrefix sName sSuffix]
   img=double(imread(fname));
   title(sName)
   seven(i,:,:) = img;
end
medmat(8,:,:) = reshape(median(seven),[1200 1632]);
%%%%%%%%%%%%

sPath = '/Users/Zoeannem/Documents/MATLAB/Fluids_Lab_3/';
sPrefix='calib_08_';
sSuffix='.tif';
iImageStart = 1;
iImageEnd= 10;

eight = zeros(10,1200,1632);
for i=iImageStart:iImageEnd
   sName=num2str(i,'%02d');
   fname=[sPath sPrefix sName sSuffix]
   img=double(imread(fname));
   title(sName)
   eight(i,:,:) = img;
end
medmat(9,:,:) = reshape(median(eight),[1200 1632]);
%%%%%%%%%%%%
%%
% Make vectors of representative pixel concentrations 
% Choose 9 pixels that can be modeled against concentrations
% [(200,200),(600,200),(1000,200),(200,800),(600,800),(1000,800),(200,1400),(600,1400),(1000,1400)]

% Create vectors of the different points across concentrations
RT = medmat(:,200,200);
RM = medmat(:,600,200);
RB = medmat(:,1000,200);
MT = medmat(:,200,800);
MM = medmat(:,600,800);
MB = medmat(:,1000,800);
LT = medmat(:,200,1400);
LM = medmat(:,600,1400);
LB = medmat(:,1000,1400);

figure()
plot(concentration,RT,'color',rand(1,3))
hold on
plot(concentration,RM,'color',rand(1,3))
hold on
plot(concentration,RB,'color',rand(1,3))
hold on
plot(concentration,MT,'color',rand(1,3))
hold on
plot(concentration,MM,'color',rand(1,3))
hold on
plot(concentration,MB,'color',rand(1,3))
hold on
plot(concentration,LT,'color',rand(1,3))
hold on
plot(concentration,LM,'color',rand(1,3))
hold on
plot(concentration,LB,'color',rand(1,3))
legend('right top','right middle','right bottom','middle top','middle middle','middle bottom','left top','left middle','left bottom')
xlabel('concentration')
hold off

conc = transpose(concentration);
%%
%Loop around every single pixel for the first 4 points of data

M = zeros(Nrows,Ncols);
B = zeros(Nrows,Ncols);

for j=1:Ncols
    for i=1:Nrows
        points = medmat(1:4,i,j);
        
        %medmat = reshape(avim(i,j,:),1,9);
        X = [ones(length(points),1),points(1:4)];
        %X = [ones(length(concentration(1:4)),1),concentration(1:4)'];
        l = (X'*X)\(X'*concentration(1:4)');
        B(i,j)=l(1); %matrix of y-intercept at every point
        M(i,j)=l(2); %matrix of slope at every point
    end
end

%%
% subtract blank from the the total average image

IM = (Cbar)-squeeze(medmat(1,:,:)); % intensity matrix with blanks subtracted

CM = IM.*M; % concentration matrix formed by multiplying slopes by average image intensity

d=0.595; %cm
% need to nondimensionalize by d
xx=linspace(0,Ncols,Ncols);
yy=linspace(0,Nrows,Nrows);

xoffset=(-1)*(200-(d*50*cal)); % find the pixel offset from start of jet to beginning of image frame

xaxis=((xx+xoffset)/cal)/d; % new x-axis where all numbers are in x/d

yaxis=((yy+xoffset)/cal)/d; % new y-axis where all numbers are in y/d
   %need to come up with a way to get the centerline as 0 

CMnorm=((CM+xoffset)/cal)/d;   
   
figure()
set(gca,'fontsize',12)
imagesc(xaxis,yaxis,CMnorm)
colorbar
axis image
xlabel('x/d')
ylabel('y/d')

%contour plot?
% figure()
% contour(CM,4)

%% spatial conversions
% x/d = 50, column = 200
% d=0.595; %cm
% cal = (521-475); %pixels/cm

% to find different x/d values
% for x/d=50, known
% xoffset=(-1)*(200-(d*50*cal)); % find the pixel offset from start of jet to beginning of image frame

% for x/d=50
xd50=(d*50*cal)-xoffset;
% for x/d=60
xd60=(d*60*cal)-xoffset;
% for x/d=70
xd70=(d*70*cal)-xoffset;
% for x/d=80
xd80=(d*80*cal)-xoffset;

%convert pixels to desired x/d values
Cxd50=CMnorm(:,round(xd50));
Cxd60=CMnorm(:,round(xd60));
Cxd70=CMnorm(:,round(xd70));
Cxd80=CMnorm(:,round(xd80));
%% make curves of all the data at 4 locations
figure()
plot(Cxd50,'color',rand(1,3))
hold on
plot(Cxd60,'color',rand(1,3))
hold on
plot(Cxd70,'color',rand(1,3))
hold on
plot(Cxd80,'color',rand(1,3))
legend('x/d=50', 'x/d=60', 'x/d=70', 'x/d=80')

%Use this to determine Co
%% find the rmss

Cp50=sqrt((Cxd50-mean(Cxd50)).^2);  %find rms based on mean 
Cp60=sqrt((Cxd60-mean(Cxd60)).^2);  %find rms based on mean 
Cp70=sqrt((Cxd70-mean(Cxd70)).^2);  %find rms based on mean 
Cp80=sqrt((Cxd80-mean(Cxd80)).^2);  %find rms based on mean 

figure()
plot(Cp50,'color',rand(1,3))
hold on
plot(Cp60,'color',rand(1,3))
hold on
plot(Cp70,'color',rand(1,3))
hold on
plot(Cp80,'color',rand(1,3))
legend('x/d=50', 'x/d=60', 'x/d=70', 'x/d=80')