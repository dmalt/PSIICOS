function [tr_pre, tr_post, CT_post, Ps, corr, Upwr, Wmat] = prepare_cond(subj_ID, condition,...
                                                        protocol_path, HM,...
                                           freq_band, time_range_pre, time_range_post,...
                                           lambda, seed_ind, SL_rnk, sig_rnk, Upwr,...
                                           cp_part, is_fast, Wmat, seed_dst_ind)
% -------------------------------------------------------------------------------- %
% Get data for prestim or poststim condition for honest_corrs script
% -------------------------------------------------------------------------------- %

    import ups.bst.LoadTrials
    import ups.conn.CrossSpectralTimeseries
    import ups.conn.DICS
    import ps.PSIICOS
    import ups.indUpperDiag2mat


    trials_pre = LoadTrials(subj_ID, condition,...
                            freq_band, time_range_pre,...
                            HM, protocol_path);

    trials_post = LoadTrials(subj_ID, condition,...
                            freq_band, time_range_post,...
                            HM, protocol_path);

    tr_pre = trials_pre.data;
    tr_post = trials_post.data;

    CT_pre = CrossSpectralTimeseries(tr_pre, true);
    CT_post = CrossSpectralTimeseries(tr_post, true);

    % CT_resh = reshape(mean(CT_pre, 2), sqrt(size(CT_pre, 1)), sqrt(size(CT_pre, 1)));
    CT_resh = reshape(mean(CT_post, 2), sqrt(size(CT_post, 1)), sqrt(size(CT_post, 1)));
    [~, Ps] = DICS(CT_resh, HM.gain, lambda, false);
    Ps = Wmat * Ps;
    % [CT_pre, Upwr] = ps.ProjectAwayFromPowerComplete(CT_pre, HM.gain, SL_rnk);
    % CT_post = ps.ProjectAwayFromPowerComplete(CT_post, HM.gain, SL_rnk);


    % CT_post_proj = ps.ProjFromCond(CT_post, CT_pre, 100);



    % CT_mean = mean(CT,2);
    % CT_mean = CT_mean / norm(CT_mean);
    %
    % Upwr = zeros(size(Upwr));
    [corr, CT_post, Upwr] = PSIICOS(CT_post, HM.gain, SL_rnk,...
                                    sig_rnk, Upwr, seed_ind,...
                                    cp_part, is_fast);

    corr.data = Wmat * corr.data;
    IND = indUpperDiag2mat(length(corr.data));

    seed_indices = IND(:,1) == seed_dst_ind | IND(:,2) == seed_dst_ind;
    IND = IND(seed_indices, :);
    corr.IND = [IND(1:seed_dst_ind, :);...
                [seed_dst_ind, seed_dst_ind];...
                IND(seed_dst_ind + 1 : end, :)];

    % corr.IND =  fix this for god's sakes!
end
