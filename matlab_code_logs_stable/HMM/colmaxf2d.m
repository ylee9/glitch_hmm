function [y,pos,mod,evidence,pmode]=colmaxf2d(d0,Q,mun,evidence,gi,R)
% Q = fokker_plank(sigma,z(t),f0(1:m(t)),fdot); t - measurement
% number,m=min(3001,ceil(sigma*z(t)^(3/2)/3^(1/2)*(2*n+1)/df)); for n=6 standard
% diviations, upper limit for 11 fdot bins is ~ 3001 ----
% mun=max(0,round((mu-f0(1))/df)+1) and mu=[f0-z(t)*fdot'];  
% L -glitch size (in freq bins)

[K,N]=size(d0);
%dynamics for Mode A
M=size(Q,2);
d=[-inf*ones(K,(M-1)/2) d0 -inf*ones(K,(M-1)/2)];
ev0=[-inf*ones(K,(M-1)/2) evidence -inf*ones(K,(M-1)/2)];
l=im2col(d,[K M],'sliding');
ev=im2col(ev0,[K M],'sliding');
normev=logsumexp(evidence(:));
w=zeros(K,N); wg=zeros(K,N); 
c=zeros(K,N); cg=zeros(K,N); 
cor=zeros(K,N); corg=zeros(K,N);
mod=zeros(K,N);
evidenceL=-inf*ones(K,N);
evidenceG=-inf*ones(K,N);
pmode = -inf*ones(K, N, 2);
e=0.0;
 A=[1-e e];
 if gi %glitch indicator
	 A=[e 1-e];
 end
 
%dynamics for Mode L
if(isempty(gi))
    for k=1:K % loop over predicted fdot
        Qn=shift_q(Q,k);
        P=repmat(Qn(:),1, N);
        [W,C]=max(l+P-max(Qn(:)));
        %indexes of maximums for each block
        cdot=modm(C,K);
        cf=ceil(C/K)-(M+1)/2+(1:N);
        %D=sum(ev.*exp(P));
        D = logsumexp(ev + P);
        %work out the correct displacement and write it into arrays
        place=mun(cdot)+(1:N)';
        [u,nu]=unique(place);

        j=interp1(u,nu,1:N,'nearest',0);
        i=find(j>0);
        w(k,i)=W(j(i));
        c(k,i)=cdot(j(i));
        cor(k,i)=cf(j(i));

        place=round(mun(k))+(1:N);
        i=find(place>0 & place<N+1);
        evidenceL(k,place(i))=D(i);
    end
    pmode(:,:,1)=w;
    [y,mod]=max(pmode,[],3);
    evidence=evidenceL;

    evidence=evidence - logsumexp(evidence(:)) + normev;
    pos(:,:,1)=c;
    pos(:,:,2)=cor;
else
    R= R - logsumexp(R(:));
    for n=1:N
        Rn=R(:,(1:N)+N-n);
        if ~isinf(logsumexp(Rn(:)))
    %   			Rn=Rn/sum(Rn(:));
            dd0=d0+Rn-max(Rn(:));
            [wm]=max(max(dd0));
            wg(:,n)=wm;
            [cgm,corgm]=find(dd0==wm,1);
            cg(:,n)=cgm; %cg
            corg(:,n)=corgm;% corg
            evidenceG(:,n)=(logsumexp(evidence + Rn, 2));
        end
    end
    pmode(:,:,2)=wg;
    [y,mod]=max(pmode,[],3);
    evidence=evidenceG;

    evidence=evidence - logsumexp(evidence(:)) + normev;
    pos(:,:,1)=cg;
    pos(:,:,2)=corg;
end

