function IND = UpperDiagToPairs(N)
% -------------------------------------------------------
% UpperDiagToPairs: for an {N x N} matrix make a mapping 
% index array IND which enumerates matrix elements above
% the diagonal. Pairs go rowwise, meaning that IND(1) is 
% a [1,2] pair, IND(2) is [1,3] and so on.
% Example:
% [i,j] = IND(p),
% where p is a number of a pair 
% -------------------------------------------------------
% FORMAT:
%   IND = UpperDiagToPairs(Nsites) 
% INPUTS:
%   N             - 
% OUTPUTS:
%   IND           - {N * (N - 1) / 2 x 2} matrix;
% ________________________________________
% Dmitrii Altukhov, dm.altukhov@ya.ru
	p = 1;
	for i = 1:N
		for j = i + 1: N
			IND(p,:) = [i,j];
			p = p + 1;
		end
	end