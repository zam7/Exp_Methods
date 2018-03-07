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
%xaxis = 1632;
%yaxis = 1200;
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
title('dye linearity test plot')
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

cal = (521-475); %pixels/cm

d=0.595; %cm
% need to nondimensionalize by d
xx=linspace(0,Ncols,Ncols);
yy=linspace(0,Nrows,Nrows);

xoffset=(-1)*(200-(d*50*cal)); % find the pixel offset from start of jet to beginning of image frame

xaxis=((xx+xoffset)/cal)/d; % new x-axis where all numbers are in x/d

centerline=557;

ymidline=max(mean(CM,1)); %this is row 558
ytopline=0+(1200-centerline);
ybottomline=0-(1200-centerline);

yaxis=((linspace(ytopline,ybottomline,1200))/cal)/d; % new y-axis where all numbers are in y/d
   %need to come up with a way to get the centerline as 0 

figure()
set(gca,'fontsize',12)
imagesc(xaxis,yaxis,CM)
colorbar
axis image
xlabel('x/d')
ylabel('y/d')
title('mean concentration field')

%contour plot?
% figure()
% contour(CM,4)
%% find Co
concmax=zeros(1,Ncols);
for i=1:Ncols
    concmax(i)=max(CM(:,i));
end

concnorm=zeros(Nrows,Ncols);
for i=1:Ncols
    concnorm(:,i)=CM(:,i)./concmax(:,i);
end

figure()
imagesc(xaxis,yaxis,concnorm,[0,1])
colorbar
axis image
xlabel('x/d')
ylabel('y/d')
title('normalized mean concentration')
%% spatial conversions
% x/d = 50, column = 200
% d=0.595; %cm
% cal = (521-475); %pixels/cm

% to find different x/d values
% for x/d=50, known
% xoffset=(-1)*(200-(d*50*cal)); % find the pixel offset from start of jet to beginning of image frame

% for x/d=50
xd50=round((d*50*cal)-xoffset);
% for x/d=60
xd60=round((d*60*cal)-xoffset);
% for x/d=70
xd70=round((d*70*cal)-xoffset);
% for x/d=80
xd80=round((d*80*cal)-xoffset);

%convert pixels to desired x/d values
Cxd50=concnorm(:,(xd50));
Cxd60=concnorm(:,(xd60));
Cxd70=concnorm(:,(xd70));
Cxd80=concnorm(:,(xd80));

%
%% make curves of all the data at 4 locations
figure()
plot(yaxis,Cxd50,'color',rand(1,3))
hold on
plot(yaxis,Cxd60,'color',rand(1,3))
hold on
plot(yaxis,Cxd70,'color',rand(1,3))
hold on
plot(yaxis,Cxd80,'color',rand(1,3))
legend('x/d=50', 'x/d=60', 'x/d=70', 'x/d=80')
xlabel('y/d')
title('normalized mean concentrations at different x/d values')

%% find the rmss

% concrms=zeros(Nrows,Ncols);
% for i=1:Ncols
%     for j=1:Nrows
%     concrms(j,i)=sqrt(((CM(j,i)-mean(CM(j,i))).^2));
%     end
% end
concrms=zeros(Nrows,Ncols);
newSum=zeros(Nrows,Ncols);
tPrefix='LIF_';
for i=iFileStart:iFileEnd
      sName=num2str(i,'%04d');
   fname=[sPath tPrefix sName sSuffix];
   img=double(imread(fname));
    for k = 1:Ncols
        for j = 1:Nrows
   concrms(j,k) = ((M(j,k)*(img(j,k)-medmat(1,j,k))) - CM(j,k))^2;
   newSum(j,k) = newSum(j,k) + concrms(j,k);
        end
    end
    test=1;
end

rmsbar = newSum/N;
concrms = sqrt(rmsbar);

concrmsnorm=zeros(Nrows,Ncols);
for i=1:Ncols
    concrmsnorm(:,i)=concrms(:,i)./concmax(:,i);
end
%%
%not normalized
figure()
imagesc(xaxis,yaxis,concrms,[0,10])
colorbar
axis image
xlabel('x/d')
ylabel('y/d')
title('root mean squared concentration')


%
Cp50=concrmsnorm(:,xd50);  %find rms based on mean %NEED TO SUBTRACT THE MEAN OF EACH PIXEL IN SPACE
Cp60=concrmsnorm(:,xd60);  %find rms based on mean 
Cp70=concrmsnorm(:,xd70);  %find rms based on mean 
Cp80=concrmsnorm(:,xd80);  %find rms based on mean 


% Cp50=sqrt((Cxd50-mean(Cxd50)).^2);  %find rms based on mean %NEED TO SUBTRACT THE MEAN OF EACH PIXEL IN SPACE
% Cp60=sqrt((Cxd60-mean(Cxd60)).^2);  %find rms based on mean 
% Cp70=sqrt((Cxd70-mean(Cxd70)).^2);  %find rms based on mean 
% Cp80=sqrt((Cxd80-mean(Cxd80)).^2);  %find rms based on mean 

figure()
plot(yaxis,Cp50,'color',rand(1,3))
hold on
plot(yaxis,Cp60,'color',rand(1,3))
hold on
plot(yaxis,Cp70,'color',rand(1,3))
hold on
plot(yaxis,Cp80,'color',rand(1,3))
xlabel('y/d')
legend('x/d=50', 'x/d=60', 'x/d=70', 'x/d=80')

%%

rhalfconc=zeros(1,Ncols);
for i=1:Ncols
    rhalfconc(i)=max(CM(:,i))/2;
end

% 
%%

%plots of normalized
%figure()
count = 0;
xdnew = 0;
ronehalf = 0;
for i = 1:20:Ncols
   xdnew20(((i-1)/20)+1) = i;  %makes a new axis
   xddist = CM(:,i)/CM(centerline,i); %normalized concentration
   
   
    %within this loop look for the r1/2 and then plot it as a function of things
   half(((i-1)/20)+1) = .5; %length 82
    for j = centerline:Nrows %558 -1200
        test1 = xddist(j);
        if test1 <= half(((i-1)/20)+1)
            ronehalf(((i-1)/20)+1) = j-centerline; %this value is in pixels right now
            count = count+1;
            break
        end
    end
    %plot(yaxis,xddist);
   %hold on
  
end
%xlabel('y/D')
%ylabel('concentration (ppb)')
%title('C at many x/D nomalized by C_0')
%%
%to find slope of ronehalf vectore

ronehalfcm=ronehalf/cal;
xaxis20=linspace(0,length(xaxis),82);
 
% WHAT ARE MY AXES AND UNITS SUPPOSED TO BE IN
figure()
PP=polyfit(xaxis20,ronehalf,1);
plot(xaxis20,ronehalf)
hold on
plot(xaxis20,PP(1)*xaxis20+PP(2))
%xlabel('x/d')
title('spreading rate')
ylabel('ronehalf')
xlabel('pixels')
hold off

%Spreading rate=slope
slope=PP(1);

%%
figure()
PP1=polyfit(xaxis(1:20:end),ronehalf,1);
plot(xaxis(1:20:end),ronehalf)
hold on
plot(xaxis(1:20:end),PP1(1)*xaxis(1:20:end)+PP1(2))
%xlabel('x/d')
title('spreading rate to find x0')
ylabel('ronehalf')
xlabel('x/d')
hold off

%Spreading rate=slope



% to find the jet virtual origin, find what is x where y=0 and 
x0=-PP(2)/slope; %this is x/d value

centerconc=CM(centerline,:);

% B IS EQUAL TO THE VELOCITY DECAY NOT THE CONCENTRATION DECAY

%need to plot (Uj/Uo)=(1/B)(x-xo/D)
PC=polyfit(xd,yplot,1);

xd=xaxis-x0;
yplot=400./centerconc;
figure()
plot(xd,yplot)
hold on
plot(xd,PC(1)*xd+PC(2))
hold off
xlabel('x/d')
ylabel('mean concentration (ppb)')
title('centerline velocity decay')


B=1/PC(1);



%CX=linspace(0,length(centerconc),length(centerconc)); %%%%what is my xaxis
% PO=polyfit(xaxis,centerconc,1);
% 
% figure()
% plot(xaxis,centerconc)
% hold on
% plot(centerconc,centerconc*PO(1)+PO(2),'m')
% hold off

%B=;

%bad stuff
for i=1:Ncols
f = CM(:,1);
val = rhalfconc(1); %value to find
tmp = abs(f-val);
[idx1,idx1] = min(tmp); %index of closest value
closest1 = f(idx1); %closest value
end
radialidx1=abs(idx1-centerline);
% 
% for i=1:Ncols
% f = CM(:,30);
% val = rhalfconc(30); %value to find
% tmp = abs(f-val);
% [idx30,idx30] = min(tmp); %index of closest value
% closest30 = f(idx30); %closest value
% end
% radialidx30=abs(idx30-centerline);
% 
% for i=1:Ncols
% f = CM(:,200);
% val = rhalfconc(200); %value to find
% tmp = abs(f-val);
% [idx200,idx200] = min(tmp); %index of closest value
% closest200 = f(idx200); %closest value
% end
% radialidx200=abs(idx200-centerline);
% 
% rhalfvector = [radialidx1 radialidx30 radialidx200];
% rxdvector = [xaxis(1) xaxis(50) xaxis(200)];
% 
% figure()
% plot(rxdvector,rhalfvector)

% 
% for i=1:Nrows
% f(i) = CM(:,i);
% val(i) = rhalfconc(i); %value to find
% tmp(i) = abs(f(:,i)-val(i));
% [idx idx] = min(tmp(i)); %index of closest value
% closest(i) = f(idx(i)); %closest value
% end
%%
% now to find the Uo(x)/Uj=B/((x-xo)/D))

