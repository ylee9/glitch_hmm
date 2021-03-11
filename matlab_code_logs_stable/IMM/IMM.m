A=[0.9 0.1; 0.9 0.1];
F=[1 delta; 0 1];F(:,:,2)=[1 delta;0 1];
q=1e-4*[delta^3/3 delta^2/2;delta^2/2 delta ];
Q(:,:,1)=q;Q(:,:,2)=q+[1 0;0 0];
par=struct('models',2,'A_m',A,'F',F,'Q',Q,'H',[1 0],'R',[1 1]*1e-4);
track=struct('x',zeros(2,1),'P',zeros(2),'mu',[0 0],'mu_hat',[0 0],'xm',zeros(2),'Pm',zeros(2,2,2),'Lambda',[0 0]);
track.mu=[0.99 0.01];track.Pm(:,:,1)=eye(2);track.Pm(:,:,2)=eye(2);
track.xm=[s(1) s(1);0 0];
clear x xm1 xm2 mu zi S
for n=1:length(s)
	track = IMM_pulse(s(n),track,par);
	x(:,n)=track.x;xm1(:,n)=track.xm(:,1);xm2(:,n)=track.xm(:,2);
	mu(n,:)=track.mu;
	zi(n,:)=track.zi;
	S(n,:)=track.S;
end