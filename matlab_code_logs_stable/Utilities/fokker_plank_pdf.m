function f=fokker_plank_pdf(x,y,mu,V)
%compute pdf given moments

c=1/(2*pi)/sqrt(det(V));
X=[x(:)-mu(1) y(:)-mu(2)];
v=sum((X/V).*X,2);
%f=c*reshape(exp(-1/2*v),size(x));
f = log(c)+reshape(-1/2*v, size(x));

