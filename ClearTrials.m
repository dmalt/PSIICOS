function ConData = ClearTrials(ConData)
% Removes field Trial from ConData{c} structure
	for c = 1:length(ConData)
		ConData{c} = rmfield(ConData{c}, 'Trials');
	end