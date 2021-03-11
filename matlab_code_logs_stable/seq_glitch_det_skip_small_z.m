function [J,BF,e,e0]=seq_glitch_det(z,fran,fdran, consec)
%
% N=length(z);
J=[];
m=0;
BF=[];
N=length(z);
p=1; %safeguard at the edges
% while true
	m=m+1;
	% zero hypothesis
	[~,~,alpha0,~,evidence0,lik] = HMM_Pulse_3dcol(z,fran,fdran,J);
    E0 = evidence0(end)
	e0(m)=E0;
	E=zeros(1,N);
	%look for biggest glitch
    tic
	parfor k=p:N-p+1
        evidence = [];
        if z(k) < 1e3
            k
            E(k) = evidence0(end)
        else
            if(consec)
                [~,evidence] = HMM_Pulse_3dcol_seq(z,fran,fdran,[J k k+1],alpha0,evidence0,lik,k);
            else
                [~,evidence] = HMM_Pulse_3dcol_seq(z,fran,fdran,[J k],alpha0,evidence0,lik,k);
            end
            k
            E(k) = evidence(end);
            evidence(end)
        end
	end
	toc
	e(m,:)=E;
    e(isnan(e)) = 0;
% 	if max(E/E0)<3
% 		break
% 	end
% 	plot(E/E0), drawnow
    [BF(m),k]=max(E-E0);
    J(m)=k;
% end

	
