subjID = 'test';
PhaseLag = pi / 20;

[HM, CT, Trials, Ctx] = ups.SimulateData(PhaseLag, 100, 0.01, 0.5, 0, true);

freqBand = [8, 12];
t_range = [0.4, 0.7];
GainSVDTh = 0.01;
isInducedOnly = true;
isLR = true;
pwr_rnk = 500;
threshold_ps = 300;
threshold_gcs = 300;
SigRnk = 0;
lambda = 100;
Upwr = [];

CT_reshape = reshape(mean(CT, 2), sqrt(size(CT,1)), sqrt(size(CT,1)));
