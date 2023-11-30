addpath(matlab_path)
addpath([matlab_path '/HMM'])
addpath([matlab_path '/Utilities'])


num_detected = 0;
for i = 0:(num_realisations-1)
    i 
    num_detected
    freqs = load([working_prefix num2str(i) '_freqs.dat']);
    fdots = load([working_prefix num2str(i) '_fdots.dat']);
    zs = load([working_prefix num2str(i) '_zs.dat']);
    kappas = load([working_prefix num2str(i) '_kappas.dat']);

    global f_fiducial fd_fiducial fdd_fiducial sigma kappa_per_toa

    %f_fiducial = load(['test_freq_fiducial_' num2str(i) '.dat']);
    %fdd_fiducial = load(['test_fdotdot_fiducial_' num2str(i) '.dat']);;
    %fd_fiducial = -load(['test_fdot_fiducial_' num2str(i) '.dat']);

    f_fiducial = f_fids(i+1);
    fd_fiducial = fd_fids(i+1);
    fdd_fiducial = fdd_fids(i+1);

    %glitch = load(['test_preglitch_toas_' num2str(i) '.dat']);

    %sigma = max(1e-21, mean(diff(fdots))/mean(zs)^0.5)
    %sigma = max(sigma, mean(diff(freqs))/mean(zs)^1.5)
    kappa_per_toa = kappas;
    glitch = glitches(i+1);
    save([working_prefix 'in_progress.mat'])
    %[J,BF,E,E0] = seq_glitch_det(zs, freqs', fdots');
    [path,gamma,alpha,beta,evidence0,s, residuals] = HMM_Pulse_3dcol_freq_noise(zs, freqs', fdots',[]);
    %[path,gamma,alpha,beta,evidence,s, residuals] = HMM_Pulse_3dcol(zs, freqs', fdots',[glitch]);
    evidence = zeros(3, length(zs)+1);
    [~,evidence(1,:)] = HMM_Pulse_3dcol_seq_freq_noise(zs,freqs',fdots',[glitch-1],alpha,evidence0,s,glitch-1);
    [~,evidence(2,:)] = HMM_Pulse_3dcol_seq_freq_noise(zs,freqs',fdots',[glitch],alpha,evidence0,s,glitch);
    [~,evidence(3,:)] = HMM_Pulse_3dcol_seq_freq_noise(zs,freqs',fdots',[glitch+1],alpha,evidence0,s,glitch+1);
    evidence(:,end) - evidence0(end)
    if sum(evidence(:,end) - evidence0(end) > log(sqrt(10))) > 0
        num_detected = num_detected + 1
    end
end
f = fopen([working_prefix 'res.dat'], 'w');
fprintf(f, "%d\n", num_detected);
fclose(f);

