uuid = getenv('UUID')
data_dir = getenv('DATA_DIR');
code_dir = getenv('MATLAB_CODE_DIR');
addpath(code_dir);

ps = parallel.Settings;
ps.Pool.AutoCreate = false;

cd(data_dir)

dname='data/';
%f0=load([dname 'overall_f_roi.dat']);
%fdran= 0.8e-14:1e-15:1.2e-14;
%fdran = linspace(0.9e-15, 1.1e-15, 3);
%fdran = linspace(0.154e-10, 0.156e-10, 5)

global sigma

global kappa_per_toa;

params_fid = fopen(['closed/params_' uuid '.dat'])
str = fread(params_fid, inf);
str = char(str');
fclose(params_fid)
params_lines = splitlines(str);
sigma_toa_line = split(params_lines{15}, ':');
sigma_toa = str2num(sigma_toa_line{2});

sigma_tn_line = split(params_lines{14}, ':');
sigma_tn = str2num(sigma_tn_line{2});
session_interval_line = split(params_lines{6});
session_interval = str2num(session_interval_line{2});
recovery_timescale_line = split(params_lines{4});
recovery_timescale = str2num(recovery_timescale_line{2});
sigma = max(5e-19, 1/(session_interval)*sigma_tn); % TODO sigma safety factor

freq_line = split(params_lines{1}, ':');
fdot_line = split(params_lines{16}, ':');
fdran = linspace(-2e-15, 2e-15, 11);
%fdran = linspace(0, 2e-15, 11)%*max(1, sigma_tn/1e-12);
%fran = linspace(-5e-7, 1e-7, 1e3) + str2num(freq_line{2});
fran = linspace(-1e-7, 1e-8, 1e3) + str2num(freq_line{2});
df = mean(diff(fran));
z=load([dname 'point_' uuid '_z.dat']);
z = z(2:end);
%toas = cumsum(z)/86400;
%toas = remove_toa_clusters(toas, session_interval);
%z = diff(toas)*86400;

global f_fiducial
global fd_fiducial
global fdd_fiducial
f_fiducial = str2num(freq_line{2});
%fd_fiducial = str2num(fdot_line{2})/3;
f_fiducial = 0;
fd_fiducial = 0;
fdd_fiducial = 0;
%fd_fiducial = str2num(fdot_line{2})/3;

err_f = 2*df*z;
err_fd = 2*0.5*mean(diff(fdran))*z.^2;
err_sigma_toa = 2*3*sigma_toa*(mean(fran) + f_fiducial)*ones(size(z)); %TODO Low sigma TOA
kappa_per_toa = 1./(err_f.^2 + err_fd.^2 + err_sigma_toa.^2);

cd(code_dir)
pwd
addpath('HMM')
addpath('Utilities')
[J,BF,E,e]=seq_glitch_det(z,fran,fdran);
[path,gamma,alpha,beta,evidence,s] = HMM_Pulse_3dcol(z,fran,fdran,J);

cd(data_dir)

fname=['results_data/point_' uuid '_path.dat'];
k=find(path(:,2));
f=fran(path(k,2));
%plot(cumsum(z),f,'o-'),drawnow
j=J
bf=BF;
n = 1;
fid=fopen(fname,'w');
for k=1:length(j)
    if j(k) ~= 1
        delta_fdot(n,k)=-fdran(path(j(k),1))+fdran(path(j(k)-1,1));
        delta_f(n,k)=fran(path(j(k),2))-(fran(path(j(k)-1,2))-z(j(k))*fdran(path(j(k),1)));
    else
        delta_fdot(n,k)=0;
        delta_f(n,k)=0
    end
    fprintf(fid,'Glitch parameters: \n');
    fprintf(fid,'Epoch = %d, Bayes Factor = %4e, Delta f = %4e, Delta_fdot =%4e\n',j(k),bf(k),delta_f(n,k),delta_fdot(n,k));
end
if isempty(j)
    fprintf(fid,'No glitch found\n');
end
fprintf(fid,'Bayesian Evidence for M_0 = %10e\n',e);
fprintf(fid,'Bayesian Evidence for M_1(j) = \n');
fprintf(fid,'%10e\n',E);
fprintf(fid,'\nFrequency path: \n');
fprintf(fid,'%1.15f\n',f);
fclose(fid);

exit
