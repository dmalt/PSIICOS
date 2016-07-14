sensors_file = '/home/dmalt/ups/channel_vectorview306.mat';
ChUsed = PickChannels('grad');
ChLoc = ReadChannelLocations(sensors_file, ChUsed);
