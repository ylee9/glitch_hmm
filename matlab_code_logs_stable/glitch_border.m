function y=glitch_border(sigma,x,d,df,ddf,N,state)
% the function computes the glitch border as a line, tangential to the
% ellipse at Mahalanobis distance ``d'' from ``state''
% sigma - std of the fdotdot noise 
% x - TOA gap
% d - how many sigmas define the border
% df, ddf - f and fdot bin sizes
% N - number of fdot bins --- must be odd
% state -  state at previous time
% output : y - border between glitch and line modes 

k=(ddf*x)/(2*df);
m=(1:N)';
%line equation
y0=-sqrt(3)/6/df*d*sigma*x^(3/2);

y=round(k*(m-state(2))+(state(1)+y0));

m=(-K:K)';
ss=real(asin((m*ddf)/(dist*sigma*z^(1/2)))) - pi/6;
nprime=n-mun(k); %displacement in previous TOA
%line equation
%shift it to the correct position
s0=ss((K+1:2*K)+1-k);
gb=3^(1/2)*dist*sigma*z^(3/2)*sin(s0)/3/df+nprime;
