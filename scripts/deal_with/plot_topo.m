ChUsed_all = PickChannels('meg');
ChUsed_1 = PickChannels('grad1');
ChUsed_2 = PickChannels('grad2');
ChUsed_grad = PickChannels('grad');
ChUsed_mag = PickChannels('mag');


hm_path = '/home/dmalt/PSIICOS_osadtchii/data/0003_pran/2/headmodel_surf_os_meg_02.mat';
HM_bst = load(hm_path);
G = HM_bst.Gain(ChUsed_all,:);
sensors_file = '/home/dmalt/ups/channel_vectorview306.mat';
MEGSensors = load(sensors_file);


ch = MEGSensors.Channel(ChUsed_all);
ChLoc = ReadChannelLocations(sensors_file, ChUsed_all);
ChLoc_grad = ReadChannelLocations(sensors_file, ChUsed_grad);

[~, G] = ReduceToTangentSpace(G, 'meg');

G_grad = G(ChUsed_grad,:);

[M,N] = size(G_grad);


% Cmm = zeros(M);
% tic
% for i = 1:N
% 	for j = 1:N
% 		Cmm = Cmm + G_grad(:,i) * G_grad(:,j)';
% 	end
% end
% toc
% save('Cmm.mat', 'Cmm', '-v7.3');

% load Cmm.mat

seed = 29;
sen_seed_xyz = ChLoc_grad(:, seed);
GridLoc =  HM_bst.GridLoc;
src_seed_ind = FindXYZonGrid(sen_seed_xyz', GridLoc);
% return;


% ---- Draw a star at seed sensor ----- %
msk = zeros(M / 2, 1);
msk(floor(seed / 2) + 1) = 1;
% ------------------------------------- %

% G_plot = G(3:3:M, src_seed_ind * 2 -1); 

% -------- Visualize topography --------------- %
% G_plot = G_grad(1:2:M, src_seed_ind * 2 ); 
% [~,~, ret_chloc] = be_viz_topo(ch(ChUsed_2), G_plot, 2, msk);
% --------------------------------------------- %

% ---------- Cross-spectrum without projection ------ %
Cmm_seed = Cmm(2:2:M, seed);
% be_viz_topo(ch(ChUsed_1), Cmm_seed, 2, msk)
% --------------------------------------------------- %
% caxis([-0.0001, 0.0001])


% ---------- Cross-spectrum without projection ------ %
% Cmm_seed = Cmm(2:2:M, seed);
% be_viz_topo(ch(ChUsed_1), Cmm_seed, 2, msk)
% % --------------------------------------------------- %
% caxis([-0.0001, 0.0001])

% ------------------- GCS -------------------------- %
G_seed_1 = G_grad(:, src_seed_ind * 2 - 1);
G_seed_2 = G_grad(:, src_seed_ind * 2);

% P_1 = eye(M) - G_seed_1 * G_seed_1';
% P_2 = eye(M) - G_seed_2 * G_seed_2';

% Cmm_p1 = P_1 * Cmm;
% Cmm_p2 = P_2 * Cmm;


% Cmm_diff_1 = Cmm - Cmm_p1;
% Cmm_diff_2 = Cmm - Cmm_p2;

% Cmm_diff_1_seed_1 = Cmm_diff_1(1:2:M, seed);
% Cmm_diff_1_seed_2 = Cmm_diff_1(2:2:M, seed);
% Cmm_diff_2_seed_1 = Cmm_diff_2(1:2:M, seed);
% Cmm_diff_2_seed_2 = Cmm_diff_2(2:2:M, seed);

% be_viz_topo(ch(ChUsed_1), Cmm_diff_2_seed_2, 2, msk)
% ------------------------------------------------------- %

% ------------- PSIICOS ------------- %
% G_seed_1 = G_grad(:, src_seed_ind * 2 - 1);
% G_seed_2 = G_grad(:, src_seed_ind * 2);

% GG_seed_1 = kron(G_seed_1, G_seed_1);
% GG_seed_2 = kron(G_seed_2, G_seed_2);

% C_P_1 = Cmm(:) - GG_seed_1 * (GG_seed_1' * Cmm(:));
% C_P_2 = Cmm(:) - GG_seed_2 * (GG_seed_2' * Cmm(:));

% Cmm_p1 = reshape(C_P_1, M, M);
% Cmm_p2 = reshape(C_P_2, M, M);

% Cmm_diff_1 = Cmm - Cmm_p1;
% Cmm_diff_2 = Cmm - Cmm_p2;

% Cmm_diff_1_seed_1 = Cmm_diff_1(1:2:M, seed);
% Cmm_diff_1_seed_2 = Cmm_diff_1(2:2:M, seed);
% Cmm_diff_2_seed_1 = Cmm_diff_2(1:2:M, seed);
% Cmm_diff_2_seed_2 = Cmm_diff_2(2:2:M, seed);

% be_viz_topo(ch(ChUsed_1), Cmm_diff_2_seed_2, 2, msk)



GG_seed = kron([G_seed_1, G_seed_2], [G_seed_1, G_seed_2]);
[u,s,v] = svd(GG_seed, 'econ');
U = u(:,1:3);
C_P = Cmm(:) - U * (U' * Cmm(:));

Cmm_p = reshape(C_P, M, M);

Cmm_diff = Cmm - Cmm_p;

Cmm_diff_seed_1 = Cmm_diff(1:2:M, seed);
Cmm_diff_seed_2 = Cmm_diff(2:2:M, seed);

be_viz_topo(ch(ChUsed_1), Cmm_diff_seed_1, 2, msk)
return

% grads = MEGSensors.Channel(ChUsed_grad);
% ReadChannelLocations(sensors_file, ChUsed_grad);
% ChLoc = ReadChannelLocations(sensors_file, ChUsed_grad);
% Loc = {grads.Loc};
% Wei = {grads.Weight};
% for i = 1:length({grads.Loc})
% 	or(:,i) = Loc{i} * Wei{i}';
% end
% quiver3(ChLoc(1,:), ChLoc(2,:), ChLoc(3,:), or(1,:), or(2,:), or(3,:))


% [or_X,or_Y] = bst_project_2d(or(1,:), or(2,:), or(3,:));
% [ch_X,ch_Y] = bst_project_2d(ChLoc(1,:), ChLoc(2,:), ChLoc(3,:));
% quiver3(ch_X, ch_Y, 0 * ch_X, or_X, or_Y, 0*or_X,'b')

% Cmm_full_1_all = Cmm_full(seed,:);
% be_viz_topo(ch, Cmm_full_1_all', 2, msk);

% 

% be_viz_topo(ch, G_f(:,src_seed_ind * 2), 2, msk);
% colormap(gca,parula(200));

% G_seed_1 = G_f(:, src_seed_ind * 2 - 1);
% G_seed_2 = G_f(:, src_seed_ind * 2 );

% P_1 = eye(M) - G_seed_1 * G_seed_1';
% P_2 = eye(M) - G_seed_2 * G_seed_2';

% Cmm_p1 = P_1 * Cmm_full;
% Cmm_p2 = P_2 * Cmm_full;


% Cmm_diff_1 = Cmm_full - Cmm_p1;
% Cmm_diff_2 = Cmm_full - Cmm_p2;

% Cmm_p1_all = Cmm_diff_1(seed,:);

% % be_viz_topo(ch, Cmm_full(:,seed), 2, msk);
% % be_viz_topo(ch, Cmm_diff_1(:,seed), 2, msk);
% % be_viz_topo(ch, Cmm_diff_2(:,seed), 2, msk);
% be_viz_topo(ch, G_seed_2(1:2:M), 3, msk);
% caxis([-0.40, 0.40])