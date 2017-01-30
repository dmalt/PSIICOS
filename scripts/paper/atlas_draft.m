subplot(2,3,2)
con_ps.Plot(0.01, linewidth, m_radius);
title('PSIICOS');

verts_hr = CtxHHR.Vertices;
faces_hr = CtxHHR.Faces;
vx_hr = verts_hr(:,1);
vy_hr = verts_hr(:,2);
vz_hr = verts_hr(:,3);

atlas = CtxHHR.Atlas(2);
mask = logical(zeros(150,1));
mask(10) = 1;

sig_scouts = atlas.Scouts(mask);
cmap = lines(length(sig_scouts));
for i = 1:length(sig_scouts)
	verts_id_sc = sig_scouts(i).Vertices;
	faces_sc = prune_faces(faces_hr, verts_id_sc);
	hold on;
	trisurf(faces_sc, vx_hr, vy_hr, vz_hr, 'EdgeColor', 'None', 'FaceColor', cmap(i,:), 'FaceAlpha', '0.3');
end
