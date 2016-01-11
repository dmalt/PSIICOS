function ConData = ClearTrials(ConData)
	for c = 1:length(ConData)
		ConData{c} = rmfield(ConData{c}, 'Trials');
	end