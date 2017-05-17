con_range = 100:100:1000;

load('2vs1.mat')
ind = 1;
for nConn = con_range
	for iSubj = 1:length(conInds_full)
		conInds_full_2vs1{iSubj} = conInds_full{iSubj}(1:nConn,:);
		conInds_real_2vs1{iSubj} = conInds_real{iSubj}(1:nConn,:);
		conInds_imag_2vs1{iSubj} = conInds_imag{iSubj}(1:nConn,:);
	end
	mean_full2vs1(ind) = ConnSimMetrics(conInds_full_2vs1, ChLoc);
	mean_real2vs1(ind) = ConnSimMetrics(conInds_real_2vs1, ChLoc);
	mean_imag2vs1(ind) = ConnSimMetrics(conInds_imag_2vs1, ChLoc);

	ind = ind + 1;
end

clear conInds_full, conInds_real, conInds_imag;

load('4vs1.mat')
ind = 1;
for nConn = con_range
	for iSubj = 1:length(conInds_full)
		conInds_full_4vs1{iSubj} = conInds_full{iSubj}(1:nConn,:);
		conInds_real_4vs1{iSubj} = conInds_real{iSubj}(1:nConn,:);
		conInds_imag_4vs1{iSubj} = conInds_imag{iSubj}(1:nConn,:);
	end
	mean_full4vs1(ind) = ConnSimMetrics(conInds_full_4vs1, ChLoc);
	mean_real4vs1(ind) = ConnSimMetrics(conInds_real_4vs1, ChLoc);
	mean_imag4vs1(ind) = ConnSimMetrics(conInds_imag_4vs1, ChLoc);

	ind = ind + 1;
end

figure; plot(con_range, mean_real2vs1, 'b-', con_range, mean_imag2vs1, 'r-',...
             con_range, mean_real4vs1, 'b-.', con_range, mean_imag4vs1, 'r-.');
axis([100,1000,0.1,0.5]);