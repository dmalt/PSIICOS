function ConData = LoadTrials(ConData, Protocol, Conditions, bLoadTrials, ProtocolDir, ChUsed)
% ---------------------------------------------------------------------------------------
% LoadTrials: load trials data from brainstorm database into ConData cell array
% ---------------------------------------------------------------------------------------
% FORMAT:
%   ConData = LoadTrials(ConData, Protocol, Conditions, bLoadTrials, ProtocolDir, ChUsed)
% INPUTS:
%   ConData        - {1 x N_conditions_total} cell array with each element 
%                    containing a structure for a single condition 
%                    recording from protocol. Gain matrices are stored in
%                    ConData{i}.HM_LR and ConData{i}.HM_HR 
%   Protocol       - brainstorm-generated structure with info about 
%                    protocol. Can be acquired via running brainstorm and 
%                    executing the command
%   Conditions     - {1 x Nconditions} cell array storing names of 
%                    conditions
%   bLoadTrials    - boolean flag; if True function will actually load
%                    trials data into ConData array. Otherwise it will
%                    just go through each condition and load number of 
%                    trials for each condition
%   ProtocolDir    - folder path to the protocol which we want to 
%                    load from
%   ChUsed         - {1 x NumOfChannelsUsed} array of indeces of
%                    channels that are left for analysis; 
%                    to use gradiometers only go for
%                    ChUsed = 1:306; ChUsed(3:3:end) = [];
% OUTPUTS:
%   ConData:
%   ConData{c}.NumTials     - scalar; number of trials for condition 'c'
%   ConData{c}.Trials       - {N_sensors_reduced x Ntimes x NumTrials}
%                             trials data
%   ConData{c}.Time         - {1 x Ntimes} discrete time array with zero 
%                             locked to stimulus. Prestimulus times are 
%                             negative.
%   ConData{c}.Fsamp        - scalar; sampling frequency for condition c
% ________________________________________________________________________
% Alex Ossadtchii ossadtchi@gmail.com, Dmitrii Altukhov, dm.altukhov@ya.ru

    Ncond = length(Conditions);
    sc = 1;
    ConditionsFound = 0;
    for c = 1:Ncond
        for s = 1:length(Protocol.Study)
            if(strcmp(Protocol.Study(s).Name,Conditions{c}))
                fprintf('Found study condition %s \n : ', Conditions{c}); 
                ConData{sc}.NumTrials = length(Protocol.Study(s).Data);
                if(bLoadTrials)
                    fprintf('Loading Trials (Max %d) : ', ConData{sc}.NumTrials); 
        %            UP = ConData{fix((sc-1)/Ncond)*Ncond+1}.UP;
                    UP = ConData{sc}.UP;
                    for iTrial = 1:ConData{sc}.NumTrials
                        % In order to save memory we use aux structure on which
                        % we do PCA and store the result in structure ConData
                        aux = load([ProtocolDir Protocol.Study(s).Data(iTrial).FileName]);
                        if iTrial == 1
                             ConData{sc}.Trials = zeros(size(UP, 1), length(aux.Time));
                             ConData{sc}.Time = aux.Time;
                             ConData{sc}.Fsamp = 1. / (aux.Time(2) - aux.Time(1));
                        end;
                        %tmp = filtfilt(b,a,(UP*aux.F(ChUsed,:))')';
                        %ConData{sc}.Trials(:,:,iTrial) = tmp(:,ind0:ind1);
                        ConData{sc}.Trials(:,:,iTrial) = UP * aux.F(ChUsed,:);
                        %ConData{sc}.Trials0(:,:,iTrial) = aux.F(ChUsed,:);
                        if iTrial > 1
                            for tt=0:log10(iTrial - 1)
                                fprintf('\b'); % delete previous counter display
                            end
                        end
                        fprintf('%d', iTrial);
                    end; % trials t
                    fprintf(' -> Done\n');
                end;
                sc = sc + 1;
             end;
        end;
    end;