function [s,F]=crlb_pulsar(f,x0,sigma);
% in terms of n under assumptiom that dist of n is normal with mean 
% mu=f*x0+1/2*x0^2*fdot and variance v=f^2*sigma^2
%FIM in this case is F_{i,j}=(d mu/d u_i)(1/v)^2(d mu/d u_j) +1/2trace{1/v dv/du_i 1/v dv/du_j }
% sigma = 2e-3 sec
% x0 - TOA gaps
%f - estimated f 
% u=[f fdot];
v=f^2*sigma^2;
dmu=[x0 1/2*x0.^2]; 
dv=[2*f*sigma^2 0];

F=dmu'*dmu/v+1/2*dv'*dv/v^2;
s=inv(F);
 
