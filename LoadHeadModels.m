function ConData = LoadHeadModels(Conditions, ProtocolDir, Protocol, bUseHR)
% -------------------------------------------------------------------------
% LoadHeadModels: loads head models for each condition recording in a 
% specified protocol from brainstorm database and stores them into 
% ConData structure
% -------------------------------------------------------------------------
% FORMAT:
%   ConData = LoadHeadModels(Conditions, ProtocolDir, Protocol, bUseHR)
% INPUTS:
%   Conditions        - {1 x Nconditions} cell array storing names of 
%                       conditions
%   ProtocolDir       - folder path to the protocol which we want to 
%                       load from
%   Protocol          - brainstorm-generated structure with info about 
%                       protocol. Can be acquired via running brainstorm and 
%                       executing the command
%                       Protocol = bst_get('ProtocolStudies', ProtocolName);
%   bUseHR            - boolean flag; if True high resolution gain matrices
%                       will be loaded as well. If false loads only low res.
%                       
% OUTPUTS:
%   ConData           - {1 x N_conditions_total} cell array with each element 
%                       containing a structure for a single condition 
%                       recording from protocol. Gain matrices are stored in
%                       ConData{i}.HM_LR and ConData{i}.HM_HR 
% __________________________________________________________________________
% Alex Ossadtchi ossadtchi@gmail.com, Dmitrii Altukhov, dm.altukhov@ya.ru
    ConditionsFound = 0;
    Ncond = length(Conditions);
    sc = 1;
    for c = 1:Ncond
        for s = 1:length(Protocol.Study)
            if(strcmp(Protocol.Study(s).Name, Conditions{c}))
                fprintf('Found study condition %s \n ', Conditions{c}); 
                for hm = 1:length(Protocol.Study(s).HeadModel)

                    if(strcmp(Protocol.Study(s).HeadModel(hm).Comment,'Overlapping spheres_HR'))
                        if bUseHR
                            ConData{sc}.HM_HR = load([ProtocolDir Protocol.Study(s).HeadModel(hm).FileName]);
                        end
                    else
                        ConData{sc}.HM_LR = load([ProtocolDir Protocol.Study(s).HeadModel(hm).FileName]);
                    end
                end;
                sc = sc + 1;
            end;
        end;
    end;
