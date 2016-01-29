function [C1, C2, C4] = PrepCondCT(ConData, N_subjects)
% -------------------------------------------------------
% PrepCondCT: split cross-spectral data from ConData for 
% PSIICOS protocol into 3 conditions: C1 (pseudowords),
% C2 (verbs) and C4 (nouns) 
% -------------------------------------------------------
% FORMAT:
%   [C1, C2, C4, C] = PrepCondCT(ConData, N_subjects) 
% INPUTS:
%   ConData        -
%   N_subjects     -
% OUTPUTS:
%   C1
%   C2
%   C4
%   C
% ________________________________________
% Dmitrii Altukhov, dm.altukhov@ya.ru

	for s=1:N_subjects
                C1{s}.CT = ConData{s}.CrossSpecTime;
                C1{s}.CTP = ConData{s}.CrossSpecTimeP;
                C1{s}.CT_Ind = ConData{s}.CrossSpecTimeInd;
                C1{s}.CT_IndP = ConData{s}.CrossSpecTimeIndP;
                C1{s}.UP = ConData{s}.UP;
                C1{s}.G = ConData{s}.G2dLRU;
                C1{s}.Upwr = ConData{s}.Upwr;


                C2{s}.CT = ConData{N_subjects + s}.CrossSpecTime;
                C2{s}.CTP = ConData{N_subjects + s}.CrossSpecTimeP;
                C2{s}.CT_Ind = ConData{N_subjects + s}.CrossSpecTimeInd;
                C2{s}.CT_IndP = ConData{N_subjects + s}.CrossSpecTimeIndP;
                C2{s}.UP = ConData{s}.UP;
                C2{s}.G = ConData{s}.G2dLRU;
                C2{s}.Upwr = ConData{s}.Upwr;

                C4{s}.CT = ConData{N_subjects * 2 + s}.CrossSpecTime;
                C4{s}.CTP = ConData{N_subjects * 2 + s}.CrossSpecTimeP;
                C4{s}.CT_Ind = ConData{N_subjects * 2 + s}.CrossSpecTimeInd;
                C4{s}.CT_IndP = ConData{N_subjects * 2 + s}.CrossSpecTimeIndP;
                C4{s}.UP = ConData{s}.UP;
                C4{s}.G = ConData{s}.G2dLRU;
                C4{s}.Upwr = ConData{s}.Upwr;
                % -- Cross-spectra projected away from other conditions --- % 
                C1{s}.CTfrom2 = ProjAwayFromCond(C1{s}.CT_Ind, C2{s}.CT_Ind);
                C1{s}.CTfrom4 = ProjAwayFromCond(C1{s}.CT_Ind, C4{s}.CT_Ind);

                C2{s}.CTfrom1 = ProjAwayFromCond(C2{s}.CT_Ind, C1{s}.CT_Ind);
                C2{s}.CTfrom4 = ProjAwayFromCond(C2{s}.CT_Ind, C4{s}.CT_Ind);

                C4{s}.CTfrom1 = ProjAwayFromCond(C4{s}.CT_Ind, C1{s}.CT_Ind);
                C4{s}.CTfrom2 = ProjAwayFromCond(C4{s}.CT_Ind, C2{s}.CT_Ind);
        end
