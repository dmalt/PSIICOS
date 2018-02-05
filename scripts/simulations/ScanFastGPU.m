% Run simulations to optimize scan fast
% with GPU
% AUTHOR: dmalt
% ------------------------------ %


subjID = 'test';
% PhaseLag = pi / 20;
PhaseLag = pi / 2 - pi / 20;
GainSVDTh = 0.01;

sim = memoize(@ups.SimulateData);
[HM, CT] = sim(PhaseLag, 100, GainSVDTh, 1, 0, false);

% profile on;
gg = rand(52, 30000, 'single');

% tic;
% [CS1,IND] = ps.PSIICOS_ScanFast(HM.gain, (mean(CT,2)));
% toc;
% tic;
% [CS2, IND] = PSIICOS_ScanFastGPU(HM.gain,(mean(CT,2)));
% toc;

% if all(abs(CS1 - CS2) < 1e-6)
%     disp('EQUAL')
% else
%     disp('NONEQUAL')
% end

tic;
[CS3,IND] = ps.PSIICOS_ScanFast(gg, (mean(CT,2)));
toc;
% profile on;
% profile on;
tic;
[CS4, IND] = PSIICOS_ScanFastGPU(gg, (mean(CT,2)), false, 1001);
toc;
% profile viewer
% toc
if all(abs(CS3 - CS4) < 1e-4)
    disp('EQUAL')
else
    disp('NONEQUAL')
end
% CT_sq = single(reshape(mean(CT,2),52, 52));
% G = single(HM.gain);
% gg = rand(52,30000, 'single');

  % aa        15000x15000            1800000000  single    complex
