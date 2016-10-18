sensors_file = '/home/dmalt/ups/channel_vectorview306.mat';
ChUsed = PickChannels('grad');
MEGSensors = load(sensors_file);

grads = MEGSensors.Channel(ChUsed);
ReadChannelLocations(sensors_file, ChUsed);
ChLoc = ReadChannelLocations(sensors_file, ChUsed);

Loc = {grads.Loc};
Wei = {grads.Weight};
for i = 1:length({grads.Loc})
	or(:,i) = Loc{i} * Wei{i}';
end

quiver3(ChLoc(1,:), ChLoc(2,:), ChLoc(3,:), or(1,:), or(2,:), or(3,:))