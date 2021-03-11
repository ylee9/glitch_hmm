function [q,r,r0,V]=fokker_plank(sigma,z,x,y,df,dfdot,mu)
%evaluate discrete distribution over rectangular grid
M=length(y);
N=length(x);

d=10;% d std
m=max(1,min(N,ceil(sigma*z.^(3/2)/3^(1/2)*(2*d+1)/df)));
% m=max(1,min(N,ceil(sigma*sqrt(z)*(2*d+1)/df)));
if (m/2)==round(m/2), m=m-1;end
xq=x(1:m);

V=[z^3/3 z^2/2;z^2/2 z]*sigma^2;
% V=[z 1;1 1]*sigma^2;
% mu=[-z*y0+x0 y0];
if nargin==6
	mu=[xq((m+1)/2) y((M+1)/2)];
end
[xx,yy]=meshgrid(xq,y);
q=fokker_plank_pdf(xx-mu(1),yy-mu(2),[0 0],V)*df*dfdot;
% fun = @(x,y)fokker_plank_pdf(x,y,mu,V);
% 
% 
% 
% for n=1:N
% 	for m=1:M
% 		q(m,n) = integral2(fun,x(n)-df/2,x(n)+df/2,y(m)-dfdot/2,y(m)+dfdot/2);
% 	end
% end
if sum(q(:))==0
	q((M+1)/2,(m+1)/2)=1;
end
q=q/sum(q(:));

%r - pdf of glitch
gb=round(z*y'/df)-1;
gb0=gb+N+m; %make zero fdot  end at  (N-1)-th f bin (offset by m to remove any effect of convolution)
L=gb0(M); % max shift due to fdot
ind=repmat(1:2*N+m,M,1); %add N zeros at the end
n=find(ind<=gb0);
r0=zeros(size(ind));
r0(n)=1;
%convolve with q
r=convn(r0,q,'same');
r(:,1:m)=[]; % remove first  m columns
  

