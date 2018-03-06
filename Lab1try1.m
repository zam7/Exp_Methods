% Lab 1
% Zoe Maisel

% what the columns of .dat file are
%  1   Ensemble counter                 (1-16777216)
%  2   Status
%  3   Velocity (Beam1|X)               (m/s)
%  4   Velocity (Beam2|Y)               (m/s)
%  5   Velocity (Beam3|Z)               (m/s)
%  6   Velocity (Beam4|Z2)              (m/s)
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

%row, column, dimesion
tic();
data(:,:,1) = load('30.9556cm.dat');
data(:,:,2) = load('31.2604.dat');
data(:,:,3) = load('31.5652cm.dat'); %assumed centerline where we started
data(:,:,4) = load('31.87cm.dat');  %based on the data, this actually looks more like centerline
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

toc(); %loading the data takes ~7 sec
% for i = 1:D
%     xvelmean(i)=mean(data(:,3,i));
%     yvelmean(i)= mean(data(:,4,i));
%     zvelmean1(i) = mean(data(:,5,i));
%     zvelmean2(i) = mean(data(:,6,i));
%     zvelmean(i) = mean((data(:,5,i) + mean(data(:,6,i))))/2;
%     corrx(i) = mean(data(:,15,i));
% end
% for i = 1:D 
%     for j = 1:R
%      xvel(j,1,i) = data(j,3,i)-xvelmean(i);
%      yvel(j,1,i) = data(j,4,i)-yvelmean(i);
%      zvel1(j,1,i) = data(j,5,i) -zvelmean1(i);
%      zvel2(j,1,i) =data(j,6,i) -zvelmean2(i);
%      zvel(j,1,i) = (zvel1(j,1,i) + zvel2(j,1,i))/2;
%     end
% end
%hists for xvel normality 2 looks really good, not that its super important...




%plot(data(:,:,1))

%%%%%%%U bar(z/r1/2)/U0%%%%%%%%%%%

%%%%%%%V bar(z/r1/2)/U0%%%%%%%%%%%

%%%%%%%W bar(z/r1/2)/U0%%%%%%%%%%%

%%%%%%%sqrt(U'2 bar((z/r1/2))/U0%%%%%%%%%%%

%%%%%%%sqrt(V'2 bar((z/r1/2))/U0%%%%%%%%%%%

%%%%%%%sqrt(W'2 bar((z/r1/2))/U0%%%%%%%%%%%

%%%%%%%(U'V' bar((z/r1/2))/U0^2%%%%%%%%%%%

%%%%%%%(U'W' bar((z/r1/2))/U0^2%%%%%%%%%%%

%%%%%%%(V'W' bar((z/r1/2))/U0^2%%%%%%%%%%%

%%%%PDF's of u'/U0, v'/U0, w'/U0 @ 3 rep locations%%%%%%%%%%

%%%%scatter plots of u'v'/U0^2, v'w'/U0^2, w'u'/U0^2 @ 3 rep locations%%%%%%%%%%

%%%%spectra of u', v', w'@ 3 rep locations%%%%%%%%%%

%------------Filter script from Cowen--------------------------------------
%first data filter 
filtertd1=data(:,3:6,1); %data that we are going to test filter with
    %take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the first data set
count1=data(:,1,1); %count of data all along length
filtertdMin=-1; %RANDOM (shouldbe centered around 0)
filtertdMax=1; %RANDOM (largest value in this set is 1.01)

[Data1,Time1]=agw_filter(filtertd1, count1, filtertdMin, filtertdMax);

%second data filter ------ where we think the midpoint is
filtertd2=data(:,3:6,2); %data that we are going to test filter with
    %take all of the data down the column, from columns 3 through 6 (X,
    %Y,Z1, Z2), of the second data set
count2=data(:,1,2); %count of data all along length

[Data2,Time2]=agw_filter(filtertd2, count2, filtertdMin, filtertdMax);

%third data filter
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

Nca=length(c); %length of the cell array, which is 17 because we are using 17 data sets that are put into a cell array

%make the cells readable as a matrix so use cell2mat function
%find the means (ubar)
for i=1:Nca
   tempdata=cell2mat(c(i));
   xmean(i)=mean(tempdata(1,:));
   ymean(i)=mean(tempdata(2,:));
   z1mean(i)=mean(tempdata(3,:));
   z2mean(i)=mean(tempdata(4,:));
   zavgmean(i)=(z1mean(i)+z2mean(i))/2;
end 

% %find the perturbations (u')
% for i=1:Nca
%    tempdata=cell2mat(c(i));
%    xprime(i,:)=(tempdata(1,:)-mean(tempdata(1,:)));
%    yprime(i,:)=(tempdata(2,:)-mean(tempdata(2,:)));
%    z1prime(i,:)=(tempdata(3,:)-mean(tempdata(3,:)));
%    z2prime(i,:)=(tempdata(4,:)-mean(tempdata(4,:)));
%    zavgprime(i,:)=(z1prime(i,:)+z2prime(i,:))/2;
%    c{i}=[xprime, yprime, z1prime, z2prime, zavgprime];
% end 
% 
% find the perturbations squared (u'^2)
% for i=1:Nca
%    tempdata=cell2mat(c(i));
%    xprime2=(tempdata(1,:)-mean(tempdata(1,:))).^2;
%    yprime2=(tempdata(2,:)-mean(tempdata(2,:))).^2;
%    z1prime2=(tempdata(3,:)-mean(tempdata(3,:))).^2;
%    z2prime2=(tempdata(4,:)-mean(tempdata(4,:))).^2;
%    zavgprime2=(z1prime2+z2prime2)/2;
% end 

U0=xmean(2); %based on the graph, we think that the second data point that we have is the actual midpoint
Uhalf=U0/2;

% xaxis = [31.5652,31.87,32.1748,32.4796,32.7844,33.394,34.0036,34.6132,35.2228,35.8324,36.442,37.0516,37.6612,38.2708,38.8804];
% f=fit(xaxis,xmean,'gauss2');
% plot(f,xaxis,xmean)

% these are the heights 
% xaxis = [31.5652,31.87,32.1748,32.4796,32.7844,33.394,34.0036,34.6132,35.2228,35.8324,36.442,37.0516,37.6612,38.2708,38.8804];
% figure();
% plot(xaxis,xmean,'c',[30 40], [0 0],'r');
% title('meanx');
% figure();
% plot(yvelmean);
% title('meany');
% figure();
% plot(zvelmean1);
% title('z1mean');
% figure();
% plot(zvelmean2);
% title('z2mean');
% figure();
% plot(corrx);
% title('corrx');
% grid;

%-----------------------Spectral Data--------------------------------------

% Creating spectra data for xprime
fs=100; %Hz
N=length(xprime); % N is the length of the xprime vector
Len=2^10; %dictates the length of each bin, in this case each will be 1024
NN=N/Len; %number of bins that are created by the length L, in this case will be 16 bins
Saatemp=zeros(Len,1); %creates a vector of zeros, length L that will be populated by for loop
Saa_avg=zeros(Len,1);

% autospectral density function for
for i=1:NN                  %index from 1 to the number of bin being created
    atemp=a(1+(i-1)*Len:i*Len); %for index 1, the bin will go from 1 to L. For index 2, the bin will go from L+1 to 2L 
    A=fft(atemp); %do fast fourier transform
    AA=conj(A).*A; %calculate the conjugate, has to be .* to make sure its element by element
    Saatemp=(AA)*(1/(Len*fs)); %gives discrete domains, and creates a vector for spectral peaks of data (equation given in notes)
    Saa_avg=Saa_avg+Saatemp; %write over and store the information for the fft from the for loop
end
Saax=Saa_avg/NN; %Saa for case with 16 bins
fMaxRes=fs/Len; %max resolution frequency 
fx=linspace(0,fs-fMaxRes,Len); %makes the x-axis for frequency vectors
