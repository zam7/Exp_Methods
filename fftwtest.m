clear all
close all

isprime(4194301) % 4194301 is a prime number

N=2^22+1 %Set N=4194305

a=rand(N,1); %create a vector of random numbers of length N

%imitalize the timer
tic
fft(a); %take the fft using fftw3
t(1)=toc %calculate the elapsed time since the call to tic
len(1)=length(a) %store the length of a in the vector len

%shorten the vector a by 1 value and repeat 
a=a(1:end-1);
N=length(a)
tic
fft(a);
t(2)=toc
len(2)=length(a)

%shorten the vector a by 1 value and repeat 
a=a(1:end-1);
N=length(a)
tic
fft(a);
t(3)=toc
len(3)=length(a)

%shorten the vector a by 1 value and repeat 
a=a(1:end-1);
N=length(a)
tic
fft(a);
t(4)=toc
len(4)=length(a)

%shorten the vector a by 1 value and repeat 
a=a(1:end-1);
N=length(a)
tic
fft(a);
t(5)=toc
len(5)=length(a)

%shorten the vector a by 1 value and repeat 
a=a(1:end-1);
N=length(a)
tic
fft(a);
t(6)=toc
len(6)=length(a)

%shorten the vector a by 1 value and repeat 
a=a(1:end-1);
N=length(a)
tic
fft(a);
t(7)=toc
len(7)=length(a)

%shorten the vector a by 1 value and repeat 
a=a(1:end-1);
N=length(a)
tic
fft(a);
t(8)=toc
len(8)=length(a)

%shorten the vector a by 1 value and repeat 
a=a(1:end-1);
N=length(a)
tic
fft(a);
t(9)=toc
len(9)=length(a)

%normalize time by the fastest calculation time (2^22)
time=t/min(t);

%results (subtrace 2^22 of of length so 0 = 4194304 and prime = -3
R=[len-2^22; time]