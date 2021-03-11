function [w,pos,mod,evidence,pmode]=colmaxf2d(d0,Q,mun,evidence,gi,R)
% Q = fokker_plank(sigma,z(t),f0(1:m(t)),fdot); t - measurement
% number,m=min(3001,ceil(sigma*z(t)^(3/2)/3^(1/2)*(2*n+1)/df)); for n=6 standard
% diviations, upper limit for 11 fdot bins is ~ 3001 ----
% mun=max(0,round((mu-f0(1))/df)+1) and mu=[f0-z(t)*fdot'];  
% L -glitch size (in freq bins)

[K,N]=size(d0);
%dynamics for Mode A
M=size(Q,2);
d=[-inf*ones(K,(M-1)/2) d0 -inf*ones(K,(M-1)/2)];
ev0=[zeros(K,(M-1)/2) evidence zeros(K,(M-1)/2)];
l=im2col(d,[K M],'sliding');
ev=im2col(ev0,[K M],'sliding');
normev=sum(evidence(:));
w=zeros(K,N); wg=zeros(K,N); 
c=zeros(K,N); cg=zeros(K,N); 
cor=zeros(K,N); corg=zeros(K,N);
mod=zeros(K,N);
 evidenceL=zeros(K,N);
 evidenceG=zeros(K,N);
e=0.0;
 A=[1-e e];
 if gi %glitch indicator
	 A=[0 1];
 end
 
%dynamics for Mode L
for k=1:K % loop over predicted fdot
	Qn=shift_q(Q,k);
	P=repmat(log(Qn(:)),1, N);
	[W,C]=max(l+P-max(log(Qn(:))));
	%indexes of maximums for each block
	cdot=modm(C,K);
	cf=ceil(C/K)-(M+1)/2+(1:N);
 	D=sum(ev.*exp(P));
	%work out the correct displacement and write it into arrays
	place=mun(cdot)+(1:N)';
    [u,nu]=unique(place);

	j=interp1(u,nu,1:N,'nearest',0);
	i=find(j>0 & j<N+1);
	w(k,i)=W(j(i));
	c(k,i)=cdot(j(i));
	cor(k,i)=cf(j(i));
	
	place=mun(k)+(1:N);
	i=find(place>0);
	evidenceL(k,place(i))=D(i);
end
pmode(:,:,1)=w+log(A(1));
pmode(:,:,2)=-inf;
if gi
	
	% for n=M+1:N
	% 	ind=max(1,n-L):n-M;
	% 	wm=max(max(d0(:,ind)));
	% 	if ~isempty(wm)
	% 		evidenceG(:,n)=mean(mean(evidence(:,ind)));
	% 		wg(:,n)=wm;
	% 		[a,b]=find(d0(:,ind)==wm);
	% 		cg(:,n)=a(end);
	% 		corg(:,n)=ind(1)+b(end)-1;
	% 	end
	% end
for n=1:N
	    Rn=R(:,(1:N)+N-n);
		if sum(Rn(:))
			dd0=d0+log(Rn)-max(log(Rn(:)));
			[wm]=max(max(dd0));
			wg(:,n)=wm;
			[cgm,corgm]=find(dd0==wm,1);
			cg(:,n)=cgm;
			corg(:,n)=corgm;
 			evidenceG(:,n)=sum(sum(evidence.*Rn));
		end
end
pmode(:,:,2)=wg+log(A(2));
%combine the scores
[w,mod]=max(pmode,[],3);
c(mod==2)=cg(mod==2);
cor(mod==2)=corg(mod==2);
end
evidence=A(1)*evidenceL+A(2)*evidenceG;


evidence=evidence/sum(evidence(:))*normev;
pos(:,:,1)=c;
pos(:,:,2)=cor;

