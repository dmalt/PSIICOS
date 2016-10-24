% ---------------- Draw connectivity similarity ---------------- %
colormap viridis;
axis equal;
% --
subplot(1,3,1);

matr_full = get_con_sim_matrix(CS_full, true);
% circularGraph(matr_full);
imagesc(matr_full)
title('Full')
% colorbar;
caxis([0,0.4])
% --
subplot(1,3,2);

matr_imag = get_con_sim_matrix(CS_imag, true);

% circularGraph(matr_imag);
imagesc(matr_imag)
title('Imaginary')
% colorbar;
caxis([0,0.4])
% --
subplot(1,3,3);

matr_real = get_con_sim_matrix(CS_real, true);
% circularGraph(matr_real);
imagesc(matr_real)

title('Real')
% colorbar;
caxis([0,0.4])
% --
% -------------------------------------------------------------- %
