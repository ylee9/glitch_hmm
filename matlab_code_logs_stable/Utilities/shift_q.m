function [Qn] = shift_q(Q,n);
% Q - gaussian computed in fokker_plank.m 
% n - \dot{f} bin 
[Kq,M]=size(Q); % K,M -odd
K=(Kq-1)/4;
 % remove extra bits on both sides and shift
Qn=Q(2*K+1+(1:2*K+1)-n,:);