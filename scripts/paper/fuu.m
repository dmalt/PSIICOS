subjID = '0003_pran';
isLR = true;
freqBand = [16,25];
t_range = [-0.5, 1];
key1 = 31;
key2 = 94;

cond_main = '2';
GainSVDTh = 0.01;
isInducedOnly = true;
protocolPath = '/home/dmalt/PSIICOS_osadtchii';

% HM = ups.LoadHeadModel(subjID, cond_main, protocolPath, isLR, GainSVDTh);
% trials = ups.LoadTrials(subjID, cond_main, freqBand, t_range, GainSVDTh, protocolPath);

tr = ups.RestoreTrDim(trials.data, HM.UP);
tr1 = squeeze(tr(key1,:,:))';
tr2 = squeeze(tr(key2,:,:))';

% ntr = 2;
for ntr = 1:100

	h1 = hilbert(tr1(ntr,:));
	h2 = hilbert(tr2(ntr,:));

	times = 1:size(tr1, 2);
	phase = unwrap(angle(h1 .* conj(h2)));
	plot(times(500:700), mod(phase(500:700), pi))
	hold on;
end
	grid
	% set(gca, 'YTick', [0:pi/9:pi])
	set(gca, 'XTick', [0:50:750])