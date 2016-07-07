function CT = ProjAwayFromCond(CT1, CT2)
% -------------------------------------------------------
% ProjAwayFromCond: project cross-spectrum CT1 from CT2
% -------------------------------------------------------
% FORMAT:
%   CT = ProjAwayFromCond(CT1, CT2) 
% INPUTS:
%   CT1        - {N_sensors_reduced ^ 2 x Ntimes} cross-spectrum
%                matrix for conditon 1
%   CT2        - {N_sensors_reduced ^ 2 x Ntimes} cross-spectrum
%                matrix for conditon 2
% OUTPUTS:
%   CT         - {N_sensors_reduced ^ 2 x Ntimes} cross-spectrum
%                matrix for condition 1 projected away from
%                condition 2
% ________________________________________
% Dmitrii Altukhov, dm.altukhov@ya.ru
	[u, s, v] = svd(CT2);
	CT = CT1 - u(:,1:6) * u(:,1:6)' * CT1;