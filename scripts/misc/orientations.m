%
% Plot orientations of Electa MEG gradiometers
% ______________________________________________

import ups.bst.PickChannels
import ups.ReadChannelLocations

% sensors_file = '/home/dmalt/ups/data/channel_vectorview306.mat';
sensors_file = '/home/dmalt/Code/matlab/utils_psiicos/channel_vectorview306.mat';
MEGSensors = load(sensors_file);
ChUsed = PickChannels(MEGSensors.Channel,'MEG GRAD');

grads = MEGSensors.Channel(ChUsed);
ReadChannelLocations(sensors_file, ChUsed);
ChLoc = ReadChannelLocations(sensors_file, ChUsed);

Loc = {grads.Loc};
Wei = {grads.Weight};
for i = 1:length({grads.Loc})
    ori(:,i) = Loc{i} * Wei{i}';
end

quiver3(ChLoc(1,:), ChLoc(2,:), ChLoc(3,:), ori(1,:), ori(2,:), ori(3,:))
