function [p,l]=von_Mises(y,k,mu)
% k=min(k,100);
p=exp(k*cos(2*pi*(y-mu)))/2/pi/besseli(0,k); % von-Mises PDF
l=sum(k*cos(2*pi*(y-mu)));
