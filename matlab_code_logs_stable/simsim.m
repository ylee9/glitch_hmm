glitch=[0 (1:10)*1e-7];
fdot0=1.11e-12:1e-14:(1.34e-12-1e-14); % odd length
load opt_par
df=1e-8;
tau=[3 180 24*3600*3]; %transient
par=[5e-8 tau];
NN=100;
M=10;
fdot=linspace(fdot0(1),fdot0(end),M);
glitch=[1.5e-7];
for m=1:10
	clear err g E gap s
	for n=1:NN
		[z,F,dF,x,N]=sim_phase([u(1) glitch(1)],[-u(2) -fdot(m)],[43],100,par);
		z0=z+randn(1,100)*1.41e-3;
		f0=min(F)-10e-8:df:max(F)+10e-8;
		[J,BF]=seq_glitch_det(z,f0,fdot0,0);
		J(1)=[];
		[path,delta,mode,psi,phi,evidence] = viterbi_Pulse_3dcol(z0,f0,fdot0,J);
		k=find(path(:,2)~=0);
		err(n)=sqrt(mean((F(k)-f0(path(k,2))).^2));
		g(n)=length(J);
		if ~isempty(J)
			s(n,1:length(J))=J;
		end
	end
	disp(m)
	meanerr(m)=mean(err);
	pdet(m)=sum(g)/NN;
% 	
	pac(m)=sum(s(:)==43 | s(:)==42 | s(:)==44)/NN;
	uns{m}=unique(s);
end
