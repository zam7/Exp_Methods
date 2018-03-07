%Basic m-file to display tiff images collected as individual tiff

clear all
close all

s_path='/Users/eac/ResearchSSD/CEE637BL/';
prefix='im289_';
suffix='.tif';
iFiles=136;
Cmax=400;
Cmin=0;

iPairSep=2; %number of images across pair (2 = adjacent images)
iStart=163; %First image
%Loop over images and display them
for i=iStart:iPairSep:iStart+iFiles-1
    sFile1=num2str(i,'%03d');
    fname=[s_path prefix sFile1 suffix];
    disp(fname)
    I1=double(imread(fname)); %loads image into a double precision real variable matrix
    sFile2=num2str(i+1,'%03d');
    fname=[s_path prefix sFile2 suffix];
    I2=double(imread(fname)); %loads image into a double precision real variable matrix
    min1=min(I1(:));
    max1=max(I1(:));
    min2=min(I2(:));
    max2=max(I2(:));
    mean1=mean(I1(:));
    mean2=mean(I2(:));
    median1=median(I1(:));
    median2=median(I2(:));
    disp(sprintf('Min1=%3d, Min2=%3d, Max1=%4d, Max2=%4d',min1,min2,max1,max2))
    disp(sprintf('Mean1=%6.2f, Mean2=%6.2f, Median1=%6.2f, Median2=%6.2f',mean1,mean2,median1,median2))
    
    figure(1)
    imagesc(I1,[Cmin Cmax])
    axis image
    colormap(gray(256))
    colorbar
    title (num2str(i))
    
    figure(2)
    imagesc(I2,[Cmin Cmax])
    axis image
    colormap(gray(256))
    colorbar
    title (num2str(i+1))
     
    figure(3)
    imagesc(I2-I1,[-Cmax Cmax])
    axis image
    colormap(gray(256))
    colorbar
    title ([sFile2 ' - ' sFile1])
    pause
end

