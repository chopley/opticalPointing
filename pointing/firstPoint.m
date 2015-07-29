day1_1 = read_arc('15-sep-2009:05:16:00', '15-sep-2009:06:00:06');
day1_2 = read_arc('15-sep-2009:06:16:00', '15-sep-2009:06:44:11');
day1_3 = read_arc('15-sep-2009:10:04:26', '15-sep-2009:13:30:00');

day2_1 = read_arc('16-sep-2009:05:24:41', '16-sep-2009:07:34:30');
day2_2 = read_arc('16-sep-2009:08:08:35', '16-sep-2009:09:00:00');
day2_3 = read_arc('16-sep-2009:09:20:37', '16-sep-2009:12:22:00');
day2_4 = read_arc('16-sep-2009:12:26:00', '16-sep-2009:13:14:55');


d = catstruct(1, [day1_1 day1_2 day1_3 day2_1 day2_2 day2_3 day2_4]);

keyboard;
eval(sprintf('display(''%d points on %d distinct sources'');', length(d.array.frame.utc), length(unique(d.antenna0.tracker.source))));

% 105 distinct sources.

save pointData2.mat 



% now we call the optical pointing routine.
[fm, ide, obs] = optical_pointing(d, 'cbassFirst');
