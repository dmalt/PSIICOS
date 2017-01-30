function faces_ = prune_faces(faces, idx)
% Remove faces if one of 3 vertex indices is in idx
% ________________________________________________

	mask_faces = ismember(faces(:,1), idx) | ismember(faces(:,2), idx) |  ismember(faces(:,3), idx);
	faces_ = faces(mask_faces,:);
end