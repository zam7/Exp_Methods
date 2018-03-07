%Question 2:

%Variables
fs = 100; %sampling frequency (Hz)
N = 32; %length of data
T = N/fs; %total time (1/Hz = seconds) 
delf=1/T; %discretized frequency
x = 10; %distance in ft between sensors

load car.dat

s1=car(:,1); %data loaded from sensor 1
s2=car(:,2); %data loaded from sensor 2
s3=car(:,3); %data loaded from sensor 3
N=32; %length of data
for i=1:N
   Raa1(i)=0; %autocorrelation of sensor 1 to sensor 1 
   Raa2(i)=0; %autocorrelation of sensor 2 to sensor 2
   Raa3(i)=0; %autocorrelation of sensor 3 to sensor 3
   Raa12(i)=0; %cross correlation of sensor 1 to sensor 2
   Raa23(i)=0; %cross correlation of sensor 2 to sensor 3
   Raa31(i)=0; %cross correlation of sensor 3 to sensor 1
   for t=1:N
      k=t+(i-1);
      if k<=N
         k=k;
      else
         k=k-N;
      end
      Raa1(i)=Raa1(i)+s1(t)*s1(k);
      Raa2(i)=Raa2(i)+s2(t)*s2(k);
      Raa3(i)=Raa3(i)+s3(t)*s3(k);
      Raa12(i)=Raa12(i)+s1(t)*s2(k);
      Raa23(i)=Raa23(i)+s2(t)*s3(k);
      Raa31(i)=Raa31(i)+s3(t)*s1(k);
   end
end
Raa1=Raa1/N; %normalize the Raa vector by length N
Raa2=Raa2/N;
Raa3=Raa3/N;
Raa12=Raa12/N;
Raa23=Raa23/N;
Raa31=Raa31/N;

%plot(Raa1,'red'); hold on; %plot all the auto and cross correlations on same graph
%plot(Raa2,'blue')
%plot(Raa3,'green')
plot(Raa12,'black','DisplayName','Cross correlation sensor 1 to 2'); 
hold on;
title('Cross correlations between sensors')
xlabel('t')
ylabel('Raa')
plot(Raa23,'yellow','DisplayName','Cross correlation sensor 2 to 3');
plot(Raa31, 'cyan','DisplayName','Cross correlation sensor 3 to 1');
legend('show')
hold off

%need to find time (index) that corresponds to max peak in the y axis
%divide by 100 to convert from 100 Hz to seconds
%subtract 1 to get data back to 
[M1,I1]=max(Raa12(:));
t1=((I1-1)/100);
[M2,I2]=max(Raa23(:));
t2=((I2-1)/100);
[M3,I3]=max(Raa31(:));
t3=((I3-1)/100);

V1=x/t1; %to find velocity, divide travel distance by time
V2=x/t2;
V3=(2*x)/t3; %there is 20 ft between sensors 1 and 2

A=[V1 V2 V3];

Vavgfs=mean(A); %ft/seconds

V1mph=V1*(3600/5280)
V2mph=V2*(3600/5280)
V3mph=V3*(3600/5280)
Vavgmph=Vavgfs*(3600/5280)



