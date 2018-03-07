
Saa = zeros(25,36,48);

for r = 1:(iFileStop/2)
% autospectral density function for the first data row
for i=1:N                  %index from 1 to the number of bin being created
%     ParvelZbed = (V(r,:,i));
    atemp=(V(r,:,i))-mean(Vvelmean(i),1); %for index 1, the bin will go from 1 to L. For index 2, the bin will go from L+1 to 2L 
    A=fft(atemp); %do fast fourier transform
    AA=conj(A).*A; %calculate the conjugate, has to be .* to make sure its element by element
    Saatemp=(AA)*(1/(Len*ks)); %gives discrete domains, and creates a vector for spectral peaks of data (equation given in notes)
    Saa(r,:,i) = (Saatemp);
end

end 