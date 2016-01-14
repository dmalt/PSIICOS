function Pairs = PlotConnections(C, ChLoc, CT_metrics_type)
	N_subjects = length(C);
	pcntg = 2 * 1e-3;
	figure
	for s = 1:N_subjects
	    %%%%%%%%%%%%%%%%% DRAWING PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	    % -------------------- Get back to real sensors from artificial ones ---------------------------- %
	    D21{s} = RestoreSensorDimension(C{s}.CT, C{s}.UP);
	    % ------------------------------------------------------------------------------------------------- %
	    if CT_metrics_type == 'real'
		    M = abs(real(D21{s}));
		elseif CT_metrics_type == 'imag'
			M = abs(imag(D21{s}))
		elseif CT_metrics_type == 'abs'
			M = abs(D21{s})
		end 	

	    %M = (ConDataBand{20+s}.wPLI-ConDataBand{s}.wPLI)-(ConDataBand{10+s}.wPLI-ConDataBand{s}.wPLI);
	    [aux, key_srt] = sort(M(:));
	    ind_max = key_srt(fix((1 - pcntg) * length(key_srt)):end);
	    % th = aux(fix((1 - pcntg) * length(key_srt)));

	    h = subplot(2, 5, s);
	    % plot3(ChLoc(1,:), ChLoc(2,:), ChLoc(3,:), '.');
	    PlotSensors(ChLoc);
	    hold on
	    Pairs{s} = [];
	    for i=1:length(ind_max)
	      [ii jj]  = ind2sub(size(D21{s}), ind_max(i));
	      Pairs{s} = [Pairs{s}; [ii jj]];
	      % plot3([ChLoc(1, ii) ChLoc(1, jj)], [ChLoc(2, ii) ChLoc(2, jj)], [ChLoc(3, ii) ChLoc(3, jj)], 'Color', 'r');
	      DrawPair(ii, jj, ChLoc);
	    end;
	    set(h, 'View', [0 90])
	    axis tight
	    axis off
	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	end