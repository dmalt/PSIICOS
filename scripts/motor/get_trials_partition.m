function [tr_cond1_part, tr_cond2_part] = get_trials_partition(tr_cond1, tr_cond2)
% ---------------------------------------------------- %
% Create random partition of trials for permutation test
% ---------------------------------------------------- %
%  FORMAT:
%    [tr_cond1_part, tr_cond2_part] = get_trials_partition(tr_cond1, tr_cond2)
%  INPUT:
%    tr_cond1        - {n_sen x n_times x n_tr} matrix of trials
%    tr_cond2        - {n_sen x n_times x n_tr} matrix of trials

%  OUTPUT:
%    tr_cond1_part
%    tr_cond2_part
% _______________________________________________________
% Dmitrii Altukhov, dm.altukhov@ya.ru

    [n_sen1, n_times1, n_tr1] = size(tr_cond1);
    [n_sen2, n_times2, n_tr2] = size(tr_cond2);

    assert(n_sen1 == n_sen2 && n_times1 == n_times2,...
           ['Trials must have the same number of sensor an time observations;',...
            ' got [', num2str(n_sen1), ',', num2str(n_times1),...
            '] and [', num2str(n_sen2), ',', num2str(n_times2), ']']);

    tr_merged = cat(3, tr_cond1, tr_cond2);
    % disp(size(tr_merged));
    ind_cond1_part = randsample(1:n_tr1+n_tr2, n_tr1);
    ind_cond2_part = setdiff(1:n_tr1+n_tr2, ind_cond1_part);
    tr_cond1_part = tr_merged(:,:,ind_cond1_part);
    tr_cond2_part = tr_merged(:,:,ind_cond2_part);
end


