function [mu,s,pi,L,f,x,F] = gmem(X,k,mu)
%
M=length(X);
m=1;
dl=10;
pi=ones(1,k)/k;
s=ones(1,k)*1e-0;
% mu=linspace(min(X),max(X),k);

while dl>1e-6
for n=1:k
	g(:,n)=pi(n)*normpdf(X,mu(n),s(n));
end
g(isnan(g))=eps;
L(m)=nansum(log(sum(g,2)));
if m>1
	dl=abs(L(m)-L(m-1));
end
disp(L(m))

if m>1e4
	break
end

gamma=g./repmat(sum(g,2),1,k);
gamma(isnan(gamma))=0;
N=sum(gamma);

mu=1./N.*sum(gamma.*repmat(X,1,k));
pi=N/M;
for n=1:k
	s(n)=max(0.01,sqrt((gamma(:,n).*(X-mu(n)))'*(X-mu(n))/N(n)));
end
if nargout>=6
	x=min(X):0.001:max(X);
	f=zeros(size(x));
	for n=1:length(mu)
		f=f+pi(n)*normpdf(x,mu(n),s(n));
	end
	plot(x,f,X,0,'xr')
	F(m)=getframe;
	drawnow
end
m=m+1;
end
% n=find(pi>1e-6);
% mu=mu(n);s=s(n);pi=pi(n);
% if nargout==6
% 	x=min(X):0.001:max(X);
% 	f=zeros(size(x));
% 	for n=1:length(mu)
% 		f=f+pi(n)*normpdf(x,mu(n),s(n));
% 	end
% end


