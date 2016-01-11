function ind = TimeAsIndex(TimeArray, Time)
% ----------------------------------------------------------
% TimeAsIndex: returnes index of instant in TimeArray
% which is the closest to specified Time
% ----------------------------------------------------------
% FORMAT:
%   ind = TimeAsIndex(TimeArray, Time) 
% INPUTS:
%   TimeArray     - {1 x Ntimes} array with descrete times
%   Time          - scalar; time instant for which we want
%                   the index of the closest element in 
%                   TimeArray  
% OUTPUTS:
%   ind           - scalar; index of an element in 
%                   TimeArray with the closest value to Time
% __________________________________________________________
% Dmitrii Altukhov, dm.altukhov@ya.ru

	if TimeArray(1) <= Time & Time <= TimeArray(end)
		[dummy, ind] = min(abs(TimeArray - Time));
	else
		fprintf('ERROR: TimeAsIndex: time %3.3f is not within [%3.3f, %3.3f] range\n',...
				 Time, TimeArray(1), TimeArray(end));
	 	ind = 0; 
	end