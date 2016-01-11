function ConDataBand = BandPassFilter(ConData, Band, TimeRange, Fsamp)
% -----------------------------------------------------------
% BandPassFilter: band-pass filter trial data, crop TimeRange
% interval and write the result into ConDataBand cell array.
% -----------------------------------------------------------
% FORMAT:
%   ConDataBand = BandPassFlter(ConData, Band, TimeRange, Fsamp) 
% INPUTS:
%   ConData        - {1 x N_conditions_total} cell array with each element 
%                    containing a structure for a single condition 
%                    recording from protocol.  
%   Band           - {1 x 2} array with lower and higher frequencies for
%                    band-pass filtering
%   TimeRange      - {1 x 2} array with ends of time interval we want to 
%                    crop from trial epochs
%   Fsamp          - scalar; data sampling frequency 
% OUTPUTS:
%   ConDataBand             - {1 x N_conditions_total} cell array storing filtered
%                             and cropped trial data.
%   ConDataBand{c}.Trials   - {N_sensors_reduced x Ntimes x NumTrials}
%                             cropped and filtered trials data
% ________________________________________________________________________
% Alex Ossadtchii ossadtchi@gmail.com, Dmitrii Altukhov, dm.altukhov@ya.ru
	fprintf('Filtering data...\n');
	N_conditions_total = length(ConData);
	[b,a] = butter(5, Band / (Fsamp / 2));
	 for sc = 1:N_conditions_total
	    for t = 1:size(ConData{sc}.Trials,3)
	        ind0 = TimeAsIndex(ConData{sc}.Time, TimeRange(1));
            ind1 = TimeAsIndex(ConData{sc}.Time, TimeRange(2));
	        T = ind1 - ind0 + 1; 
	        tmp = filtfilt(b, a, (ConData{sc}.Trials(:,:,t))')';
	        ConDataBand{sc}.Trials(:,:,t) = tmp(:,ind0:ind1);
	    end;
	    sc;
	end;
