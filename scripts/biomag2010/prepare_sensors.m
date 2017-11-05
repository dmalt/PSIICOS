function meg_loc_3d = prepare_sensors(ch_path)
    ch = load(ch_path);
    ch = ch.Channel;
    ch_used = ups.bst.PickChannels(ch, 'MEG');
    ch_meg = ch(ch_used);
    n_ch = length(ch_meg);
    meg_loc_3d = zeros(n_ch, 3);
    for i_ch = 1:n_ch
        meg_loc_3d(i_ch,:) = mean(ch_meg(i_ch).Loc,2);
    end
end
