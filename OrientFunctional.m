function f = OrientFunctional(phi, A)
% -------------------------------------------------------
% OrientFunctional: compute coherence between two sources
% for a given orientation of topography vectors.
% Orientations are given by u = [cos(phi(1)); sin(phi(1))],
% v = [cos(phi(2)); sin(phi(2))]
% -------------------------------------------------------
% FORMAT:
%   f = OrientFunctional(phi, A) 
% INPUTS:
%   phi        - {2 x 1} array with angles of orientation
%                for two dipoles at two different sources 
%   A          - {2 x 2} matrix of coherence between two 
%                sources for a loose orientation of 
%                dipoles
% OUTPUTS:
%   f          - minus coherence on given orientations 
% ________________________________________
% Dmitrii Altukhov, dm.altukhov@ya.ru

	u = [cos(phi(1)); sin(phi(1))];
	v = [cos(phi(2)); sin(phi(2))];
	f = -(u' * real(A) * v) ^ 2 - (u' * imag(A) * v) ^ 2;
