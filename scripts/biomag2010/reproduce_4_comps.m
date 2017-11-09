main

time_range = [0, 1];

alpha_band = [8,12];
beta_band = [16,24];

freq_band = beta_band;
% freq_band = alpha_band;

ltr     = memoize(@load_trials);
calc_CT = memoize(@ups.conn.CrossSpectralTimeseries);
proj    = memoize(@ps.ProjectAwayFromPowerComplete);
msvd    = memoize(@svd);

tr_filt= ltr(data_path, time_range, freq_band, HM);
% CT = calc_CT(tr_filt(:,:,10:15), true);
% CT = ups.GetFakeCT(size(HM.gain,1), 1200);
CTs = ups.bootstrap_CT(tr_filt, 100, true, 3);
% profile on;
Upwr = ps.GetUpwrComplete(HM.gain, pwr_rnk);
for i_resamp = 1:length(CTs)
    tic
    CT_proj = CTs{i_resamp} - Upwr * (Upwr' * CTs{i_resamp});
    [u,s,v] = msvd(real(CT_proj));
    % [u,s,v] = msvd(imag(CT_proj));
    % [u,s,v] = msvd(imag(CT));
    [CS, IND] = ps.PSIICOS_ScanFast(HM.gain, u(:,i_comp));
    con_inds{i_resamp} = ups.threshold_connections(CS, threshold, IND);
    toc
end


    con = ups.Connections(con_inds, HM, CtxHHR);
    % con.CT = CT_proj;
    con.CT = CT_proj;
    con_av = con.Average();

    % figure; con.PlotViews(0.1, 2);


    % plot(v(1,:))
    % [corr, Cpvec, Upwr] = 
    % colors = ups.plt.GetColors(3);

    % figure;
    % c = zeros(4,3);
    % for i = 1:4
    %     c(i,:) = colors(i + 1).RGB;
    %     plot(v(:,i), 'LineWidth', 2, 'Color', colors(i + 1).RGB);
    %     hold on;
    % end

    % legend('comp1', 'comp2', 'comp3', 'comp4');
