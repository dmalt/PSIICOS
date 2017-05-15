%
% Plot orientations of Electa MEG gradiometers
% ______________________________________________

import ups.PickChannels
import ups.ReadChannelLocations

sensors_file = '/home/dmalt/ups/data/channel_vectorview306.mat';
ChUsed = PickChannels('grad');
MEGSensors = load(sensors_file);

grads = MEGSensors.Channel(ChUsed);
ReadChannelLocations(sensors_file, ChUsed);
ChLoc = ReadChannelLocations(sensors_file, ChUsed);

Loc = {grads.Loc};
Wei = {grads.Weight};
for i = 1:length({grads.Loc})
	ori(:,i) = Loc{i} * Wei{i}';
end

quiver3(ChLoc(1,:), ChLoc(2,:), ChLoc(3,:), ori(1,:), ori(2,:), ori(3,:))
