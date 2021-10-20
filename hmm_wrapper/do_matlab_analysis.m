addpath(matlab_path)
addpath([matlab_path '/HMM'])
addpath([matlab_path '/Utilities'])

%addpath('/fred/oz022/ldunn_glitch/matlab_code_logs/')
%addpath('/fred/oz022/ldunn_glitch/matlab_code_logs/HMM')
%addpath('/fred/oz022/ldunn_glitch/matlab_code_logs/Utilities')


%global f_fiducial fd_fiducial fdd_fiducial sigma kappa_per_toa

freqs = load([working_prefix 'freqs.dat']);
fdots = load([working_prefix 'fdots.dat']);
zs = load([working_prefix 'zs.dat']);
kappa_per_toa = load([working_prefix 'kappas.dat']);

%f_fiducial = load('test_freq_fiducial.dat');
%fdd_fiducial = load('test_fddot_fiducial.dat');
%fd_fiducial = -load('test_fdot_fiducial.dat');

%sigma = max(1e-21, mean(diff(test_fdots))/mean(test_zs)^0.5)
%kappa_per_toa = test_kappas;

save([working_prefix 'init.mat']);

[J,BF,E,E0] = seq_glitch_det(zs, freqs', fdots');
[path,gamma,alpha,beta,evidence,s, residuals] = HMM_Pulse_3dcol(zs, freqs', fdots',[J]);

f = fopen([out_prefix 'res.dat'], 'w');
fprintf(f, "J: %d\n", J);
fprintf(f, "BF: %.16f\n", BF);
fclose(f);

f = fopen([out_prefix 'fdot_path.dat'], 'w');
fprintf(f,"%d\n", fdots(path(:,1)));
fclose(f);

f = fopen([out_prefix 'f_path.dat'], 'w');
fprintf(f,"%d\n", freqs(path(:,2)));
fclose(f);

f = fopen([out_prefix 'final_evidence.dat'], 'w');
fprintf(f, "%.16f\n", evidence);
fclose(f);

for n = 1:size(path,1) g(:, n) = logsumexp(gamma(:, :, n),1); end;
f = fopen([out_prefix 'f_posterior.dat'], 'w');
fprintf(f, "%.16f\n", g);
fclose(f);

for i = 1:(numel(J)+1)
    E(i,:) = E(i,:) - E0(i);
end
f = fopen([out_prefix 'bfs.dat'], 'w');
fprintf(f, "%.16f\n", E);
fclose(f);

save([out_prefix 'analysis.mat'])
