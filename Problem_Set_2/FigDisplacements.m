function FigDisplacements(Xdis,Ydis,Nx,Ny)
%FigDisplacements(Xdis,Ydis) - function to plot displacement
%results from an image pair

h=findobj('type','figure');
if sum(h==10)>0
 close(10)
end
figure(10)
set(10,'position',[100 500 500 170])
 
subplot(1,2,1)
imagesc(real(Xdis),[-Nx/2 Nx/2])
axis image
colorbar
title('Xdis')
subplot(1,2,2)
imagesc(real(Ydis),[-Ny/2 Ny/2])
axis image
colorbar
title('Ydis')

%{
figure(11)
U=abs(Xdis)<(Nx/2);
U=U.*Xdis;
V=abs(Ydis)<(Ny/2);
V=V.*Ydis;
quiver(U, V)
%}




  