NOrder  = 6;
[mean_full, std_full, CS_full] = ups.ConnSimMetrics(conInds_full, ChLoc, NOrder);
[mean_real, std_real, CS_real] = ups.ConnSimMetrics(conInds_real, ChLoc, NOrder);
[mean_imag, std_imag, CS_imag] = ups.ConnSimMetrics(conInds_imag, ChLoc, NOrder);
