% Test three different ways to recover orientation of dipole
subj_num = 2;
trial_num = subj_num;

CT = C{subj_num}.CT;
G = ConData{trial_num}.G2dLRU;
UP = ConData{trial_num}.UP;

G = UP' * G;
N_sen_reduced = sqrt(size(CT,1));
CT = reshape(CT, N_sen_reduced, N_sen_reduced);
CT = UP'* CT * UP;
ind_loc_1 = 10;
ind_loc_2 = 51;

range1 = 2 * ind_loc_1 - 1 : 2 * ind_loc_1;
range2 = 2 * ind_loc_2 - 1 : 2 * ind_loc_2;

g1 = G(:, range1);
g2 = G(:, range2);

A = g1' * CT * g2;
[ua,sa,va] = svd(A);

B = [real(A), imag(A)];
[ub,sb,vb] = svd(B);

D = [real(A); imag(A)];
[ud,sd,vd] = svd(D);

oss_number = norm(ub(:,1)' * A * vd(:,1));
gross_number = sa(1,1);


[u,v] = FindOr(CT, g1, g2);
dmalt_number = norm(u' * A * v);

fprintf('oss_number: %s, gross_number: %s, dmalt_number: %s \n', ...
	    oss_number, gross_number, dmalt_number);
