function vv_fl = flip_components(vv)


    cc = corrcoef(vv);
    [n_times,n_comps] = size(vv);

    for ii = 1:n_comps
        flip_sign = ((cc(:,ii) > 0) - 0.5) * 2;
        flip_matr = diag(flip_sign);
        vv_fl = vv * flip_matr;
        cc = corrcoef(vv_fl);
        mcc(ii) = mean(cc(:));
    end

    [val,key] = max(mcc);

    flip_sign = ((cc(:,key) > 0) - 0.5) * 2;
    flip_matr = diag(flip_sign);
    vv_fl = vv * flip_matr;

    % cc_fl = corrcoef(vv_fl)
    % figure;
    % subplot(1,2,1);
    % caxis([-1,1]);
    % imagesc(cc);
    % subplot(1,2,2);
    % imagesc(cc_fl);
    % caxis([-1,1]);

end
