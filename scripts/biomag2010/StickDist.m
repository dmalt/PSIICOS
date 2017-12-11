function dist = StickDist(stick1, stick2)
    d1 = norm(stick1(1,:) - stick2(1,:));
    d2 = norm(stick1(2,:) - stick2(2,:));
    dist = mean([d1,d2]);
end
