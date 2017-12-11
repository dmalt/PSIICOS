% ---------------------------------------------- %
% This scripts is designed to test different
% strategies of finding active connections from
% cross-spectrum timeseries by the means of SVD.
% ---------------------------------------------- %
% AUTHOR: dmalt
% DATE: Wed Nov 15 14:23:16 MSK 2017


% ------------------ setup  constants ----------------- %
% paths
home = getenv('HOME');
bst_path = [home,'/Documents/MATLAB/bst/brainstorm_db'];
data_path = [home, '/Data/mentrot/MentalRotationDeLange/preprocessed'];
subj_ID = 'biomag2010';
protocol_path = [bst_path,'/mentrot'];
condition = 'raw';

% head model
isLR = true;
% isLR = false;
GainSVDTh = 0.001; % results in 45 components
ch_type = 'MEG';

% load head model
HM = ups.bst.LoadHeadModel(subj_ID, condition, protocol_path, isLR, GainSVDTh, ch_type);

% load cortical surfaces
[Ctx, CtxHR, CtxHHR] = ups.bst.GetCtx(subj_ID, protocol_path);
% ------------------------------------------------------------ %

% ---------- set network locations ---------- %
G = HM.gain;
i = 100;
j = 1221;
gi = HM.gain(:,i * 2);
gi = gi / norm(gi);
gj = HM.gain(:,j * 2);
gj = gj / norm(gj);
con = ups.Bundles([i,j], HM, CtxHHR);

% con.PlotViews(0.2);
% ------------------------------------------- %

% -------------- setup network temporal profile -------------- %
% network activity is designed to have bell-shaped envelope
% and slowly drifting in time phase shift

T = 1000; % number of time samples
t = 1:T;
mu = 500; % central sample
tseries_env = exp(-1e-10 * ((mu - t) .^ 4));
% omega = 0.5 * pi / T; % frequency
omega = 0.1 * pi / T; % frequency
delta_phi =  pi / 4 ; % network phase shift
tseries_carrier = exp(1i * (omega * t + delta_phi));
tseries = tseries_carrier .* tseries_env;

% figure;
% plot(t, real(tseries), t, imag(tseries), t, abs(tseries),...
%      t, real(tseries_carrier), t, imag(tseries_carrier),...
%      'linewidth', 2);
% legend('Re', 'Im', 'Env', 'Re(carrier)', 'Im(carrier)');
% % -------------------------------------------------------- %

% ------------- prepare cross-spectrum timeseries ------------- %
Gij = kron(gi, gj);
Gji = kron(gj, gi);
CT = (Gij + Gji) * real(tseries) + 1i * (Gij - Gji) * imag(tseries);

figure;
subplot(1,3,1);
imagesc(real(CT));
title('real');
subplot(1,3,2);
imagesc(imag(CT));
title('imag');
subplot(1,3,3);
imagesc(abs(CT));
title('abs');
% ------ try different strategies of CT stacking + svd ------ %

% 1) decompose CT as is
%    ^^^^^^^^^^^^^^^^^^
[u,s,v] = svd(CT);

con_1 = con;
[Cs_re, IND] = ps.PSIICOS_ScanFast(HM.gain, real(u(:,1)));
con_inds_re = ups.threshold_connections(Cs_re, 20, IND);
con_1.conInds{2} = con_inds_re;
[Cs_im, IND] = ps.PSIICOS_ScanFast(HM.gain, imag(u(:,1)));
con_inds_im = ups.threshold_connections(Cs_im, 20, IND);
con_1.conInds{3} = con_inds_im;
[Cs_all, IND] = ps.PSIICOS_ScanFast(HM.gain, (u(:,1)));
con_inds_all = ups.threshold_connections(Cs_all, 20, IND);
con_1.conInds{4} = con_inds_all;

% figure('Name', 'SVD of CT as is');
% con_1.PlotViews(0.1, [8;1;1;1],...
%                 [0.004; 0.001; 0.001; 0.001], [],...
%                 [1; 0.7; 0.7; 0.7]);

figure;
plot(t, s(1,1) * real(v(:,1)) / sqrt(2), '-.',...
     t, s(2,2) * real(v(:,2)) / sqrt(2), '-.',...
     t, sqrt((s(2,2) * real(v(:,2))) .^ 2 + (s(1,1) * real(v(:,1))) .^ 2) / sqrt(2), '-.',...
     t, tseries_env,...
     t, real(tseries),...
     t, imag(tseries),...
     'linewidth', 2);
title('SVD of CT as is');
legend('re_1', 're_2', 'sqrt(re_1^2 + re_2^2)', 'orig\_env', 'orig\_re', 'orig\_im');

% 2) do svd on real and imag parts separately
%    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

[u_re,s_re,v_re] = svd(real(CT));
[u_im,s_im,v_im] = svd(imag(CT));

con_2 = con;

[Cs_re_2, IND] = ps.PSIICOS_ScanFast(HM.gain, u_re(:,1));
ci_re_2 = ups.threshold_connections(Cs_re_2, 20, IND);
con_2.conInds{2} = ci_re_2;
[Cs_im_2, IND] = ps.PSIICOS_ScanFast(HM.gain, u_im(:,1));
ci_im_2 = ups.threshold_connections(Cs_im_2, 20, IND);
con_2.conInds{3} = ci_im_2;

% figure('Name', 'SVD of Re and Im separately');
% con_2.PlotViews(0.2, [8;1;1],...
%                 [0.004; 0.001; 0.001], [],...
%                 [1;0.7;0.7]);

figure;
plot(t, s_re(1,1) * v_re(:,1) / sqrt(2), '-.',...
     t, s_im(1,1) * v_im(:,1) / sqrt(2), '-.',...
     t, sqrt((s_re(1,1) * v_re(:,1)) .^ 2 + (s_im(1,1) * v_im(:,1)) .^ 2) / sqrt(2), '-.',...
     t, tseries_env,...
     t, real(tseries),...
     t, imag(tseries),...
     'linewidth', 2);
title('SVD of Re and Im separately');
legend('re', 'im', 'env', 'orig\_env', 'orig\_re', 'orig\_im');

% SUMMARY
% -------
% In the end of the day both strategies are capable of finding one single network in cross-spectrum timeseries.
% Decomposition of CT as is recovers temporal profiles of real and imaginary parts by putting activation into
% first and second right singular vectors v whereas the second strategy recovers real and imaginary parts
% separately.
% On the  down side the first strategy sometimes puts in the first component the envelope of activation and the second component in that case doesn't show meaningful activation
