function CT_restored = RestoreSensorDimension(CT, UP)
% -------------------------------------------------------
% RestoreSensorDimension: transform cross-spectrum CT
% back to real sensors from artificial ones.
% -------------------------------------------------------
% FORMAT:
%   CT_restored = RestoreSensorDimension(CT, UP) 
% INPUTS:
%   CT        - {N_sensors_reduced ^ 2 x Ntimes}
%               cross-spectrum matrix on artificial sensors
%   UP        - {N_sensors x N_sensors_reduced} matrix;
%               its columns are first N_sensors_reduced
%               left singular vectors of gain matrix.
% OUTPUTS:
%   CT_restored   - {N_sensors ^ 2 x Ntimes} cross-spectrum 
%                   matrix on real sensors
% ______________________________________________________
% Dmitrii Altukhov, dm.altukhov@ya.ru
	for t = 1:size(CT,2)
		CT_reshaped = reshape(CT(:,t), sqrt(size(CT, 1)), sqrt(size(CT,1)));
	    CT_sq_restored = UP' * CT_reshaped * UP;
	    CT_restored(:,t) = CT_sq_restored(:);
	end
