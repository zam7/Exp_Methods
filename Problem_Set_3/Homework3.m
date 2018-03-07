% Homework 3
%%
% Problem 1

% create a random array of numbers to start working with 
% 30 seconds at 0.5 Hz --> 16 samples

% need to make a moving average of the 16 samples - take an average of 5
% samples, shifting down by 1 so there is overlap

uprime=transpose([0.70,0.65,0.6,0.55,0.50,0.45,0.40,0.35,0.3,0.25,0.2,0.15,0.1,0.05,0.03,0.01]);
upsqr=(uprime.^2);

%to take moving average that will create 15 averaged data points
% first point will stay the same, second one will be averaged by 3,
% third-thirteenth will be averaged by 5, 14th by 3, 15th stays the same

temp=zeros(16,1);
for i=1:16
    if i==1
        temp(i)=upsqr(i);
    elseif i==2
        temp(i)=(upsqr(i-1)+upsqr(i)+upsqr(i+1))/3;
    elseif 2<i && i<15
        temp(i)=(upsqr(i-2)+upsqr(i-1)+upsqr(i)+upsqr(i+1)+upsqr(i+2))/5;
    elseif i==15
        temp(i)=(upsqr(i-1)+upsqr(i)+upsqr(i+1))/3;
    elseif i==16
        temp(i)=upsqr(i);
    end
end

% now need to take square root of averages
urms=sqrt(temp);

% to do a bootstrap on this analysis, we need to do a bootstrap of all of
% the different sets that we took averages on  - want to populate

btep=zeros(16,1);
for i=1:16
    if i==1
        btep(i)=upsqr(i);
        for j=1:15
        btep(i,j)=btep(i);
        end
    elseif i==2
        A=[upsqr(i-1) upsqr(i) upsqr(i+1)];
            btep(i,:)=repmat(A,1,5);
    elseif 2<i && i<15
        B=[upsqr(i-2) upsqr(i-1) upsqr(i) upsqr(i+1) upsqr(i+2)];
        for j=1:15
        btep(i,:)=repmat(B,1,3);
        end
    elseif i==15
        A1=[upsqr(i-1) upsqr(i) upsqr(i+1)];
            btep(i,:)=repmat(A1,1,5);
    elseif i==16
        btep(i)=upsqr(i);
        for j=1:15
        btep(i,j)=btep(i);
        end
    end
end

bootmat=zeros(16,1000);
for i=1:16
    bootmat(i,:)=sqrt(sort(bootstrap(1000,@mean,btep(i,:))));
end

maxline=bootmat(:,950);
q=0.1; %quiescient limit

timely=linspace(0,30,16);
figure()
plot(timely,maxline,'color',rand(1,3))
hold on
plot(timely,urms,'color',rand(1,3))
hold on
refline(0,q)
ylim([0 0.7])
legend('bootstrap','mean')
hold off

x=maxline-q;

for i=1:16
if x(i)>0
    valuehigh=x(i);
    timehigh=timely(i);
elseif x(i)<0
    valuelow=x(i);
    timelow=timely(i);
    break
end
end

tqmet=timelow; %choose more conservative time to prevent messing with quiescence

%%
% Problem 2

load final_signal.dat;

freq=1024; %data was sampled uniformly over 1 second and there are 1024 data points so it was sampled at 1024 Hz
t=(0:1/1023:1);
x=final_signal; 
y=hilbert_6370(x);

figure()
plot(t,real(y),'color',rand(1,3));
hold on
plot(t,imag(y),':','color',rand(1,3));
hold off
legend('x(t)', 'H[x(t)]')
xlabel('t')
ylabel('x(t) and H[x(t)]')

figure()
A=abs(y);
phi=(unwrap(2*atan(imag(y)./real(y))))/2;
plot(t,A)
title('Amplitude')
xlabel('t')
ylabel('A(t)')

figure()
plot(t,phi)
title('Phase Unwrap')
xlabel('y')
ylabel('theta (radians)')

%ylim([0 10])

figure()
f=1/(2*pi)*diff(phi)/(t(2)-t(1));
plot(t(1:1024/2),f(1:1024/2))
title('Frequency')
xlabel('t')
ylabel('f')


%%
% Problem 3
load final_auv.dat;

UAV=final_auv;
T=0.0516; %s
fs=1/T; %Hz
N=length(UAV);
fN=fs/2;

Len=N; %dictates the length of each bin, in this case each will be 4096
NN=N/Len; %number of bins that are created by the length L, in this case will be 4 bins
Saatemp=zeros(Len,1); %creates a vector of zeros, length L that will be populated by for loop
Saa_avg=zeros(Len,1); %creates a vector of zeros, length L that will be populated by for loop

for i=1:NN                  %index from 1 to the number of bin being created
    atemp=UAV(1+(i-1)*Len:i*Len); %for index 1, the bin will go from 1 to L. For index 2, the bin will go from L+1 to 2L 
    A=fft(atemp); %do fast fourier transform
    AA=conj(A).*A; %calculate the conjugate, has to be .* to make sure its element by element
    Saatemp=(AA)*(1/(Len*fs)); %gives discrete domains, and creates a vector for spectral peaks of data (equation given in notes)
    Saa_avg=Saa_avg+Saatemp; %write over and store the information for the fft from the for loop
end
Saa=Saa_avg/NN; %Saa for case with 4 bins
fMaxRes=fs/Len; %max resolution frequency 
f=linspace(0,fs-fMaxRes,Len); %makes the x-axis for frequency vectors

figure()
plot(f,Saa,'DisplayName','Original data')
xlabel('f')
ylabel('Saa')
legend('show')
%hold on
%plot (f1,Saa1,'yellow');
xlim ([0 fN/3])
%%
% Problem 4

%%