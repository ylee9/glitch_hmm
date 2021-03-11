dname='data/';
d=dir(dname);
f0=load([dname d(3).name]);
d(1:4)=[];
for n=1:2:length(d)
	str((n+1)/2,:)=d(n).name(1:42);
end
N=size(str,1);
fdran=9.8e-13:1e-15:10.2e-13;
df=1e-9;
ff=[floor(f0(1)/df) ceil(f0(2)/df)]*df;
fran=ff(1):df:ff(2);

for n=1:N
	disp(str(n,:))
	disp(n)
% 	f0=load(['data/' str(n,:) '_f_roi.dat']);
% 	ff=round(f0(2)/df)*df; %highest value of f0 is correct starting frequency
	z=load([dname str(n,:) '_z.dat']);
% 	maxf=1e-12*sum(z); % max frequency displacement
% 	fran=(ff-10*df-maxf):df:(ff+10*df);
% 	if length(fran)>3000 continue;end
	[J,BF,E,e]=seq_glitch_det(z,fran,fdran);
	[path] = HMM_Pulse_3dcol(z,fran,fdran,J);
	fname=['results_data/' str(n,:) '_path.dat'];
	f=fran(path(:,2));
	plot(cumsum(z),f,'o-'),drawnow
	j=J
	bf=BF;
 fid=fopen(fname,'w');
    for k=1:length(j)
        if bf(k)>3
            delta_fdot(n,k)=-fdran(path(j(k),1))+fdran(path(j(k)-1,1));
            delta_f(n,k)=fran(path(j(k),2))-(fran(path(j(k)-1,2))-z(j(k))*fdran(path(j(k),1)));
            fprintf(fid,'Glitch parameters: \n');
            fprintf(fid,'Epoch = %d, Bayes Factor = %4e, Delta f = %4e, Delta_fdot =%4e\n',j(k),bf(k),delta_f(n,k),delta_fdot(n,k));
        end
    end
    if isempty(j)
        fprintf(fid,'No glitch found\n');
    end
    fprintf(fid,'Bayesian Evidence for M_0 = %10e\n',e);
    fprintf(fid,'Bayesian Evidence for M_1(j) = \n');
    fprintf(fid,'%10e\n',E);
    fprintf(fid,'\nFrequency path: \n');
    fprintf(fid,'%1.9f\n',f);
    fclose(fid);

end	
