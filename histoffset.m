function [Histo, offst, dt] = histoffset(dif,binsize);

dt = max(dif)-min(dif)+1;       % finding number of time bins
Histo = zeros(dt,2);            % initialising vector

[n, edges] = histcounts(dif,dt);    % creates histogram with dt time bins
Histo(:,1) = floor(edges(2:end));   % specifies time bin differences
Histo(:,2) = n;                     % counts in each time bin difference
clear n edges

range = 3/binsize;
[maxi, midi] = max(Histo(floor(dt/2-range):floor(dt/2+range),2));
midi = midi + floor(dt/2-range)-1;
offst = Histo(midi,1);