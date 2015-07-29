d1 = read_arc('25-sep-2009:03:45:00', '25-sep-2009:05:40:00');
d2 = read_arc('25-sep-2009:06:36:10', '25-sep-2009:08:02:00');
d3 = read_arc('25-sep-2009:08:04:30', '25-sep-2009:08:56:30');
d4 = read_arc('25-sep-2009:09:04:00', '25-sep-2009:09:48:30');
d5 = read_arc('25-sep-2009:09:50:58', '25-sep-2009:09:59:00');
d6 = read_arc('25-sep-2009:10:01:45', '25-sep-2009:10:56:35');
d7 = read_arc('25-sep-2009:10:58:35', '25-sep-2009:11:04:25');
d8 = read_arc('25-sep-2009:11:07:55', '25-sep-2009:12:02:55');
d8 = read_arc('25-sep-2009:12:04:55', '25-sep-2009:12:52:40');
d9 = read_arc('25-sep-2009:12:56:40', '25-sep-2009:13:00:00');
d10= read_arc('25-sep-2009:13:01:57', '25-sep-2009:13:25:00');

d = catstruct(1, [d1 d2 d3 d4 d5 d6 d7 d8 d9 d10]);

keyboard;
eval(sprintf('display(''%d points on %d distinct sources'');', length(d.array.frame.utc), length(unique(d.antenna0.tracker.source))));

% 105 distinct sources.

save pointDataSec.mat 
keyboard;

% now we call the optical pointing routine.
[fm, ide, obs] = optical_pointing(d, 'cbassFirst');
