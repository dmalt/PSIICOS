

[rH_ind, lH_ind, isConnected(2)] = tess_hemisplit(Ctx_dst);
faces = Ctx_dst.Faces;
% faces_hr_l_orig = prune_faces(faces_hr, lH_ind);

idx = lH_ind;
mask_faces = ismember(faces(:,1), idx) | ismember(faces(:,2), idx) |  ismember(faces(:,3), idx);
faces_ = faces(mask_faces,:);
