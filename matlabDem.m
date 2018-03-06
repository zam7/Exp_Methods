% Creating a vector of integer numbers
i=0:10

% An ending semicolon suppresses the screen dump
i=0:10;

% An increment (step size) can be used
x=0:0.1:10;

% To find the number of elements in a vector (1-D array)
N=length(x);

% Functional expressions can be written in standard form
y=sin(pi*x);

% To plot this
plot(x,y)

% To label it
xlabel('x')
ylabel('y')
title('Plot of sin(x)')

% How about a second funtional dependence on x?

% y2=0.01*x^2

% Error 'Matrix must be square' occurs because Matlab is trying to multiply the 1-D array
% (which is not square) by itself.  We want Matlab to square each individual element in x
% thus we need to use the '.^' opperator.  In Matlab a '.' before an operator indicates to
% perform the operation on individual elements in an array

y2=0.01*x.^2;

% And plotting this

plot(x,y2)

% What if I want both plots on the same figure?
% Two ways to do this, I can just type

hold

% And now I will continue to write on top of my existing plot

plot(x,y)

% And to release the plot so that it will no longer write over

hold

% A second way, and the more common way in this case, is to use

plot(x,y,x,y2)

% Note that Matlab automatically choses the colors so that the two linetypes are different
% We can control linetype and color

plot(x,y,'r--',x,y2,'g:')

% to see all the options use the Matlab help feature

help plot

% Note that help is available for all functions

% we can add a grid

grid

% And specify a legend

h=plot(x,y,'r--',x,y2,'g:')
legend(h,'y','y2')
grid

% or thicken a line

set(h(2),'linewidth',2)

% In general we can adjust a lot of parameters of the individual lines - to see them all

get(h(1))

get(h(2))

% Note that the linewidth parameter, the linestyle parameter and the color parameter are 
% different for each line

% If we want to change the axis limits we can

axis([0 5 -1 1])

% Where the bracketed term is a 4 component vector with [ xmin xmax ymin ymax]

axis([0 10 -1.5 1.5])

% In general we can build arrays by typing directly

a = [1 2 3 4 5 6]

% 2D arrays can be direclty entered as

b=[1 2 3;4 5 6]

% We can get the dimensions of the 2D array using the size command

size(a)
size(b)
[nrows ncols]=size(b)

% In 2D arrays the first index is row number and the second index is column number

% Note that matlab is case sensitive

a
%A is not equal to a, try it

% Three other standard plot features are (note use of pause)

semilogx(x,y2)
pause
semilogy(x,y2)
pause
loglog(x,y2)
grid
pause

% Note ctrl-C will interrupt any command

% Multiple plots can be put on a single page with subplot

subplot(3,1,1)
semilogx(x,y2)

subplot(3,1,2)
semilogy(x,y2)

subplot(3,1,3)
loglog(x,y2)

% .m files are essentially scripts that can be run at the command line

% A .m file to look at the grade distribution for problem set 1 is given in pr1.m