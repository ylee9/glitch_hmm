function [z,f0,fdot0,x,N]=sim_phase(f,fdot,glitch,L,par)
%%%
lambda=[1e2 1e5 1e6]; %par for exp distribution
pgap=[0.65 0.95];%pdf of mixture
sigma=1e-16;
v=1e-3;% std of the timing noise
sec=-1/f(1):1e-4:1/f(1);
f0=f(1);
fdot0=fdot(1);
k=1;
Zg=[];
for n=1:L
	b=0;
	gap=rand;
	Z=3600;
% 	if (gap<pgap(1))
% 		Z=exprnd(lambda(1),1,1)+500;
% 	elseif gap>=pgap(1) & gap<pgap(2)
% 		Z=exprnd(lambda(2),1,1);
% 	else
% 		Z=exprnd(lambda(3),1,1);
% 		if isempty(glitch)
% 		glitch=n
% 		end
% 	end
    if  n==101
		Z=1*24*3600;
	end
	t=Z+sec;	
	if any(glitch==n)
		Zg=0;
		k=k+1;
		b=par(1)*sum(exp(-(sec+1)./par(2:4)'));%transient term starts at around Z sec after previous measurement
		ff=f0(n)+fdot0(n)*t+f(k)+b;%glitch
		fdot0(n+1)=fdot(k);
	else
		if ~isempty(Zg)
			b=par(1)*sum(exp(-(sec+Zg+Z+1)./par(2:4)')-exp(-(sec+Zg+1)./par(2:4)'));%transient term
			Zg=Zg+Z;
		end
		ff=f0(n)+fdot0(n)*t+b;
		fdot0(n+1)=fdot0(n)+randn*sigma*sqrt(Z); %add noise to fdot0
	end
	phi=ff.*t-1/2*fdot0(n+1).*t.^2;
	mx=peakdet(cos(2*pi*phi),1e-3);
	if isempty(mx)
		[~,mx]=max(cos(2*pi*phi));
		mx=mx(:);
	end
	[~,m]=min(abs(sec(mx(:,1))));
	
	N(n+1)=phi(mx(m,1));
	f0(n+1)=ff(mx(m,1));
	z(n+1)=t(mx(m,1));
end
x=cumsum(z);
z(1)=[];f0(1)=[];fdot0(1)=[];x(1)=[];N(1)=[];