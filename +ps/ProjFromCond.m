function CT = ProjFromCond(CT1, CT2, rnk)
% -------------------------------------------------------
% Project cross-spectrum CT1 from CT2
% -------------------------------------------------------
% FORMAT:
%   CT = ProjFromCond(CT1, CT2) 
% INPUTS:
%   CT1        - {N_sensors_reduced ^ 2 x Ntimes} cross-spectrum
%                matrix for conditon 1
%   CT2        - {N_sensors_reduced ^ 2 x Ntimes} cross-spectrum
%                matrix for conditon 2
%   rnk        - integer scalar; rank of projector from CT2
% OUTPUTS:
%   CT         - {N_sensors_reduced ^ 2 x Ntimes} cross-spectrum
%                matrix for condition 1 projected away from
%                condition 2
% ________________________________________
% Dmitrii Altukhov, dm.altukhov@ya.ru
	if nargin < 3
		rnk = 6;
	end

	[u, ~, ~] = svd(CT2);
	if rnk <= size(u, 2)
		CT = CT1 - u(:,1:rnk) * u(:,1:rnk)' * CT1;
	else
		error('InputError: ValueOutOfRange', 'rnk is bigger then size(CT2,1)');
	end
end