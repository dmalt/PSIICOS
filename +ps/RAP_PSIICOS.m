function corr = RAP_PSIICOS(C, G2dU, SL_rnk, sig_rnk,...
                     Upwr, seed_ind, cp_part, is_fast, n_rap)
% --------------------------------------------------
% HAHHHAHHAHAHAHA!!!!!!!!!!!!!!!!!!
% --------------------------------------------------
% ___________________________________________________
% Dimtrii Altukhov, dm.altukhov@ya.ru
%
    import ps.PSIICOS

for i_rap = 1:n_rap
     [corr{i_rap}, C, Upwr] = PSIICOS(C, G2dU, SL_rnk,...
                                      sig_rnk, Upwr,...
                                      seed_ind, cp_part,...
                                      is_fast);

     if i_rap == 1
         Upwr_ini = Upwr;
         SL_rnk = 0;
     end
     % Find strongest connection and project C away from it
     [~, indmax] = max(corr{i_rap}.data);
     ij_max = corr{i_rap}.IND(indmax,:);

     g_i_xy = G2dU(:, ij_max(1) * 2 - 1 : ij_max(1) * 2);
     g_j_xy = G2dU(:, ij_max(2) * 2 - 1 : ij_max(2) * 2);

     G_ij = kron(g_i_xy, g_j_xy);
     G_ji = kron(g_j_xy, g_i_xy);

     GG = [G_ij, G_ji];
     GG = GG - Upwr_ini * (Upwr_ini' * GG);
     C = C - GG * (pinv(GG) * C);
     % Upwr = zeros(size(Upwr)); % don't do projection on 2-nd and subseq. iters
end
