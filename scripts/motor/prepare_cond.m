function [tr, CT, Ps, corr, Upwr, Wmat] = prepare_cond(subj_ID, condition, protocol_path, HM,...
                                           freq_band, time_range,...
                                           lambda, seed_ind, SL_rnk, sig_rnk, Upwr,...
                                           cp_part, is_fast, Wmat)
% -------------------------------------------------------------------------------- %
% Get data for prestim or poststim condition for honest_corrs script
% -------------------------------------------------------------------------------- %

    import ups.LoadTrials
    import ups.CrossSpectralTimeseries
    import ups.DICS
    import ps.PSIICOS


    trials = ups.LoadTrials(subj_ID, condition,...
                            freq_band, time_range,...
                            HM, protocol_path);

    tr = trials.data;
    CT = CrossSpectralTimeseries(tr, true);

    CT_resh = reshape(mean(CT, 2), sqrt(size(CT, 1)), sqrt(size(CT, 1)));

    [~, Ps] = DICS(CT_resh, HM.gain, lambda, false);

    % CT_mean = mean(CT,2);
    % CT_mean = CT_mean / norm(CT_mean);

    [corr, ~, Upwr] = PSIICOS(CT, HM.gain, SL_rnk,...
                           sig_rnk, Upwr, seed_ind,...
                           cp_part, is_fast);
    corr.data = Wmat * corr.data;
end
