function [q,r,mu]=fokker_plank_b(sigma,z,x,y,df,dfdot,gi)
%evaluate discrete backward distribution over rectangular grid
Kq=length(y);
N=length(x);
K=(Kq-1)/4;
d=3;% d std
%m=max(1,min(N,ceil(sigma*z.^(3/2)/3^(1/2)*(2*d+1)/df)));
m=max(1,min(N,ceil(sigma*z.^(1/2)*(2*d+1)/df)));
if (m/2)==round(m/2), m=m-1;end
xq=x(1:m);

%V=[z^3/3 z^2/2;z^2/2 z]*sigma^2;
V = [z 0; 0 1]*sigma^2;
s=[xq((m+1)/2) y((Kq+1)/2)];%position in the center
[xx,yy]=meshgrid(xq,y);
q=fokker_plank_pdf(xx-s(1),yy-s(2),[0 0],V);
if sum(q(:))==0
	q((Kq+1)/2,(m+1)/2)=1;
end
q=q - logsumexp(q(:));
r=[];
mu=round(z*y'/df); % mean backward displacement in bins
mu=mu(K+(1:2*K+1));
if gi
	%r - pdf of glitch
	gb=round(-z*y'/df)+1;
	gb0=gb+N; %make zero fdot  end at  (N+1)-th f bin (offset by m to remove any effect of convolution)
	ind=repmat(1:N*2+m,Kq,1); %add N zeros at the front
	n=find(ind>=gb0 & ind<gb0+1e4);
	r0=zeros(size(ind));
	r0(n)=1;
	%convolve with q
	r=convn(r0,exp(q),'same');
	r(:,N*2+m-(1:m))=[]; % remove last  m columns
	r=r(K+(1:2*K+1),:);%ignore extra fdran

	r=r/sum(r(:));%max(r(:))*max(q(:));\
    r = log(r);
end

  

