main;
band_names = {'theta_band', 'alpha_band', 'beta_band', 'lowgamma_band', 'gamma_band'};
data_dir = './bootstrap_data/';
save_dir = './pics/bootstrap_pics/';
n_comps = 4;
i_band = 3;
i_comp = 1;

n_bands  = length(band_names);

for i_band = 1:n_bands
% for i_band = 4:4%length(band_names)
    s = load([data_dir, band_names{i_band}, '.mat']);
    for i_comp = 1:n_comps
        con_re = ups.Bundles(s.con_inds_re{i_comp}, HM, CtxHHR);
        cm_re = con_re.Merge();
        cm_re.conInds = {ups.CoOrientCluster(cm_re.conInds{1}, HM.GridLoc)};
        [~,~,var_re] = ups.ClustAverage(cm_re.conInds, HM.GridLoc);

        con_re = ups.Bundles(s.con_inds_re{i_comp}, HM, CtxHHR);
        cm_re = con_re.Merge();
        cm_re.conInds = {ups.CoOrientCluster(cm_re.conInds{1}, HM.GridLoc)};
        [~,~,var_re] = ups.ClustAverage(cm_re.conInds, HM.GridLoc);
        inv_var_re(i_comp,i_band) = sqrt(1 / (var_re{1} + var_re{2}));



        con_im = ups.Bundles(s.con_inds_im{i_comp}, HM, CtxHHR);
        cm_im = con_im.Merge();
        cm_im.conInds = {ups.CoOrientCluster(cm_im.conInds{1}, HM.GridLoc)};
        [~,~,var_im] = ups.ClustAverage(cm_im.conInds, HM.GridLoc);

        con_im = ups.Bundles(s.con_inds_im{i_comp}, HM, CtxHHR);
        cm_im = con_im.Merge();
        cm_im.conInds = {ups.CoOrientCluster(cm_im.conInds{1}, HM.GridLoc)};
        [~,~,var_im] = ups.ClustAverage(cm_im.conInds, HM.GridLoc);
        inv_var_im(i_comp,i_band) = sqrt(1 / (var_im{1} + var_im{2}));
        % cc = ups.Bundles(s.con_inds_re{i_comp}, HM, CtxHHR);
        % con_im = ups.Bundles(s.con_inds_im{i_comp}, HM, CtxHHR);
    end
end

c = categorical(band_names);
inv_var = zeros(5,8);
for i = 1:n_bands
    for j = 1:n_comps
        inv_var(i, 2 * j - 1 : 2 * j) =  [inv_var_re(j,i), inv_var_im(j,i)];
    end
end

c = categorical({'comp1, theta',...
                 'comp2, theta',...
                 'comp3, theta',...
                 'comp4, theta',...
                 'comp1, alpha',...
                 'comp2, alpha',...
                 'comp3, alpha',...
                 'comp4, alpha',...
                 'comp1, beta',...
                 'comp2, beta',...
                 'comp3, beta',...
                 'comp4, beta',...
                 'comp1, gamma1',...
                 'comp2, gamma1',...
                 'comp3, gamma1',...
                 'comp4, gamma1',...
                 'comp1, gamma2',...
                 'comp2, gamma2',...
                 'comp3, gamma2',...
                 'comp4, gamma2'}, 'Ordinal', false);
c = reordercats(c, {'comp1, theta',...
                    'comp2, theta',...
                    'comp3, theta',...
                    'comp4, theta',...
                    'comp1, alpha',...
                    'comp2, alpha',...
                    'comp3, alpha',...
                    'comp4, alpha',...
                    'comp1, beta',...
                    'comp2, beta',...
                    'comp3, beta',...
                    'comp4, beta',...
                    'comp1, gamma1',...
                    'comp2, gamma1',...
                    'comp3, gamma1',...
                    'comp4, gamma1',...
                    'comp1, gamma2',...
                    'comp2, gamma2',...
                    'comp3, gamma2',...
                    'comp4, gamma2'})


b = bar(c,[inv_var_re(:), inv_var_im(:)]);
ylim([0,21]);
legend('Re', 'Im');
ylabel('1 / var')
% bar(c,inv_var);