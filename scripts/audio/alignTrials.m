function tr_aligned = alignTrials(trials_erp_data, trials_data)
	ERP = CalcERP(trials_erp_data);
	[u, s, v] = svd(ERP);
	u1 = u(:,1);

	[nSen, nTimes, nTrials] = size(trials_erp_data);

	for iTr = 1:nTrials
		tr(iTr,:) = squeeze(u1' * squeeze(trials_erp_data(:,:,iTr)));
	end

	mtr = mean(tr);
	[~, argmin] = min(mtr);

	[~, argmin_all] = min(tr,[], 2);	
	argmin_all(argmin_all < 400 | argmin_all >  500) = argmin;
	delta = argmin_all - argmin;

	starts = 250 + delta;
	% tr_aligned = zeros(nTrials, 300)
	for iTr = 1:length(delta)
		tr_aligned(:,:, iTr) = trials_data(:,starts(iTr) : starts(iTr) + 350, iTr);
	end

end