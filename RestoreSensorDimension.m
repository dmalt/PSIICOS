function CT_restored = RestoreSensorDimension(CT, UP)
% -------------------------------------------------------
% RestoreSensorDimension: transform cross-spectrum CT
% back to real sensors from artificial ones.
% -------------------------------------------------------
% FORMAT:
%   CT_restored = RestoreSensorDimension(CT, UP) 
% INPUTS:
%   CT        - {N_sensors_reduced x N_sensors_reduced}
%               cross-spectrum matrix on artificial sensors
%   UP        - {N_sensors x N_sensors_reduced} matrix;
%               its columns are first N_sensors_reduced
%               left singular vectors of gain matrix.
% OUTPUTS:
%   CT_restored   - {N_sensors x N_sensors} cross-spectrum 
%                   matrix on real sensors
% ______________________________________________________
% Dmitrii Altukhov, dm.altukhov@ya.ru

	CT_reshaped = reshape(CT, sqrt(size(CT, 1)), sqrt(size(CT,1)));
    CT_restored = UP' * CT_reshaped * UP;