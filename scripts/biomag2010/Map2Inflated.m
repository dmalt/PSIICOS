xyz = HM.GridLoc;
xyz_hr = CtxInfl.Vertices;
xyz_mapped = zeros(size(xyz));
tic;
fprintf('Wow, look ---> ');
for i_vert_lr = 1:length(HM.GridLoc)
    ind = ups.FindXYZonGrid(xyz(i_vert_lr,:), xyz_hr);
    xyz_mapped(i_vert_lr,:) = xyz_hr(ind,:);
    if i_vert_lr > 1
        for jj = 0:log10(i_vert_lr-1)
            fprintf('\b');
        end
    end
    fprintf('%d',i_vert_lr);
end
fprintf('\n');
toc
