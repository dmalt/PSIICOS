% Fix bootstrap data after an errorneous run of reproduce_4_comps
% AUTHOR: dmalt
% ______________________________________________________

main

time_range = [0, 1];


Upwr = ps.GetUpwrComplete(HM.gain, pwr_rnk);

for i_band = 1:length(bands)
% freq_band = beta_band;
freq_band = bands{i_band};
load(['./bootstrap_data/', band_names{i_band}, '.mat']);
clear con_inds_re;
% i_comp = 1;
    for i_resamp = 1:length(CTs)
        CT_proj = CTs{i_resamp} - Upwr * (Upwr' * CTs{i_resamp});
        % CT_proj = CT - Upwr * (Upwr' * CT);
        [u_re,~,~] = svd(real(CT_proj));

        for i_comp = 1:4
            tic;
            [CS_re, IND] = ps.PSIICOS_ScanFast(HM.gain, u_re(:,i_comp));
            toc;

            con_inds_re{i_comp}{i_resamp} = ups.threshold_connections(CS_re, threshold, IND);
        end
end
save([band_names{i_band}, '_fix.mat'], 'con_inds_re', '-v7.3');
end
return


load('./theta_band_fix.mat')
