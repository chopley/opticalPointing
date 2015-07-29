function d = checkImages(filename)

d = [];

txt= sprintf('[s w] = unix(''ls --color=never -lh /home/cbassuser/fits/%s*'');', filename);

eval(txt);


entVals = strfind(w, '.fit');

for m=1:length(entVals)
  entry{m} = w(entVals(m)-15:entVals(m)+3);
end


for m=1:length(entry)
  eval(sprintf('thisFits = fitsread2(''/home/cbassuser/fits/%s'');', entry{m}));
imagesc(thisFits);

eval(sprintf('title(''%s'');', entry{m}(10:15)));
pause(0.1);

end
keyboard;



% now we read in the data
%d = read_arc('11-feb-2010:04:43:31', '11-feb-2010:11:34:00');
%badTimes = [date2mjd(2010, 2, 11, 5, 55) date2mjd(2010, 2, 11, 6, 05);
%  date2mjd(2010, 2, 11, 6, 10) date2mjd(2010, 2, 11, 6, 17);
%  date2mjd(2010, 2, 11, 6, 18) date2mjd(2010, 2, 11, 6, 20);
%  date2mjd(2010, 2, 11, 6, 59) date2mjd(2010, 2, 11, 7, 01);
%  date2mjd(2010, 2, 11, 7, 36) date2mjd(2010, 2, 11, 7, 37);
%  date2mjd(2010, 2, 11, 8, 04) date2mjd(2010, 2, 11, 8, 08);
%  date2mjd(2010, 2, 11, 8, 36) date2mjd(2010, 2, 11, 8, 41);
%  date2mjd(2010, 2, 11, 9, 02) date2mjd(2010, 2, 11, 9, 10);
%  date2mjd(2010, 2, 11, 9, 30) date2mjd(2010, 2, 11, 9, 36, 30);
%  date2mjd(2010, 2, 11, 9, 54) date2mjd(2010, 2, 11, 9, 57);
%  date2mjd(2010, 2, 11,10, 15) date2mjd(2010, 2, 11, 10,19);
%  date2mjd(2010, 2, 11,10, 34) date2mjd(2010, 2, 11, 10, 38);
%  date2mjd(2010, 2, 11,10, 51) date2mjd(2010, 2, 11, 11, 01);
%  date2mjd(2010, 2, 11,11, 06) date2mjd(2010, 2, 11, 11, 10);
%  date2mjd(2010, 2, 11,11, 16) date2mjd(2010, 2, 11, 11, 22);
%  date2mjd(2010, 2, 11,11, 24) date2mjd(2010, 2, 11, 11, 27);
%  date2mjd(2010, 2, 11,11, 30) date2mjd(2010, 2, 11, 11, 33);
%  ];
%ind = ones(size(d.array.frame.utc));
%for m=1:size(badTimes,1)
%  f = find(d.array.frame.utc>badTimes(m,1) & ...
%      d.array.frame.utc<badTimes(m,2));
%  ind(f) = 0;
%end
%d = framecut(d, logical(ind));

    
  
  


return

% data Apr 26, 2011
d = read_arc('26-Apr-2011:04:42:54',  '26-Apr-2011:12:51:00');
% data to remove
badTimes = [date2mjd(2011, 04, 26, 10, 25) date2mjd(2011, 04, 26, 10, 29);
  date2mjd(2011, 04, 26, 10, 36) date2mjd(2011, 04, 26, 10, 40);
  date2mjd(2011, 04, 26, 11, 52) date2mjd(2011, 04, 26, 12, 31);];
ind = ones(size(d.array.frame.utc));
for m=1:size(badTimes,1)
  f = find(d.array.frame.utc>badTimes(m,1) & ...
      d.array.frame.utc<badTimes(m,2));
  ind(f) = 0;
end
d = framecut(d, logical(ind));
dc = framecut(d, d.array.frame.features>0);
save opt_point_apr26.mat dc




% data from October 14, 2010
d = read_arc('14-oct-2010:05:05:54', '14-oct-2010:13:38:00');

% data to remove
badTimes = [date2mjd(2010, 10, 14, 3, 33) date2mjd(2010, 10, 14, 3, 37);
  date2mjd(2010, 10, 14, 4, 20) date2mjd(2010, 10, 14, 4, 23);
  date2mjd(2010, 10, 14, 5, 35) date2mjd(2010, 10, 14, 5, 38);
  date2mjd(2010, 10, 14, 5, 42) date2mjd(2010, 10, 14, 5, 45);
  date2mjd(2010, 10, 14, 6, 40) date2mjd(2010, 10, 14, 6, 42);
  date2mjd(2010, 10, 14, 7, 31) date2mjd(2010, 10, 14, 7, 35);
  date2mjd(2010, 10, 14, 9, 30) date2mjd(2010, 10, 14, 9, 35);
  date2mjd(2010, 10, 14,10, 11) date2mjd(2010, 10, 14,10, 16);
  date2mjd(2010, 10, 14,11, 03) date2mjd(2010, 10, 14,11, 18);
  date2mjd(2010, 10, 14,11, 59) date2mjd(2010, 10, 14,12, 14);
  date2mjd(2010, 10, 14,12, 45) date2mjd(2010, 10, 14, 12,54);
  date2mjd(2010, 10, 14,13, 14) date2mjd(2010, 10, 14, 13,16);
  date2mjd(2010, 10, 14,13, 21) date2mjd(2010, 10, 14, 13, 24);
  ];
ind = ones(size(d.array.frame.utc));
for m=1:size(badTimes,1)
  f = find(d.array.frame.utc>badTimes(m,1) & ...
      d.array.frame.utc<badTimes(m,2));
  ind(f) = 0;
end
d = framecut(d, logical(ind));
dc = framecut(d, d.array.frame.features>0);
save point_oct14.mat dc;

%------------------------------------------------------------

d = read_arc('15-Oct-2010:02:54:02', '15-oct-2010:10:00');

badTimes = [date2mjd(2010, 10, 15,10, 50) date2mjd(2010, 10, 15,10, 53);
  date2mjd(2010, 10, 15,10, 57) date2mjd(2010, 10, 15,11, 00)];
ind = ones(size(d.array.frame.utc));
for m=1:size(badTimes,1)
  f = find(d.array.frame.utc>badTimes(m,1) & ...
      d.array.frame.utc<badTimes(m,2));
  ind(f) = 0;
end
d = framecut(d, logical(ind));
dc = framecut(d, d.array.frame.features>0);

d = read_arc('15-Oct-2010:10:00:00', '15-oct-2010:13:39');
badTimes = [date2mjd(2010, 10, 15,10, 50) date2mjd(2010, 10, 15,10, 53);
  date2mjd(2010, 10, 15,10, 57) date2mjd(2010, 10, 15,11, 00)];
ind = ones(size(d.array.frame.utc));
for m=1:size(badTimes,1)
  f = find(d.array.frame.utc>badTimes(m,1) & ...
      d.array.frame.utc<badTimes(m,2));
  ind(f) = 0;
end
d = framecut(d, logical(ind));
dc2 = framecut(d, d.array.frame.features>0);

dc = catstruct(1, [dc dc2]);

ind = diff(dc.array.frame.utc)*24*60*60>2;
ind = [1; ind];

dc = framecut(dc, ind);

save point_oct15.mat dc;


%------------------------------------------------------------
d1 = read_arc('25-Nov-2010:01:17:51', '25-nov-2010:04:27:47');
d2 = read_arc('25-Nov-2010:05:55:00', '25-nov-2010:09:12:02');
d = catstruct(1, [d1 d2]);
clear d1; 
clear d2;
badTimes = [date2mjd(2010, 11, 25,03, 30) date2mjd(2010, 11, 25,03, 35);
  date2mjd(2010, 11, 25, 7, 45) date2mjd(2010, 11, 25, 7, 52)];
ind = ones(size(d.array.frame.utc));
for m=1:size(badTimes,1)
  f = find(d.array.frame.utc>badTimes(m,1) & ...
      d.array.frame.utc<badTimes(m,2));
  ind(f) = 0;
end
d = framecut(d, logical(ind));
dc = framecut(d, d.array.frame.features>0);

ind = diff(dc.array.frame.utc)*24*60*60>2;
ind = [1; ind];

dc = framecut(dc, logical(ind));

save point_nov25.mat dc;


%------------------------------------------------------------

d = read_arc('30-May-2011:06:01:10',  '30-May-2011:12:02:00');
badTimes = [...
  date2mjd(2011, 05, 30, 7, 31) date2mjd(2011, 05, 30, 7, 40);
  date2mjd(2011, 05, 30, 7, 43) date2mjd(2011, 05, 30, 7, 48);
  date2mjd(2011, 05, 30, 7, 52) date2mjd(2011, 05, 30, 8, 06);
  date2mjd(2011, 05, 30, 8, 21) date2mjd(2011, 05, 30, 8, 45);
  date2mjd(2011, 05, 30, 8, 50) date2mjd(2011, 05, 30, 8, 59);
  date2mjd(2011, 05, 30, 9, 08) date2mjd(2011, 05, 30, 9, 15);
  date2mjd(2011, 05, 30, 9, 25) date2mjd(2011, 05, 30, 9, 34);
  date2mjd(2011, 05, 30, 9, 39) date2mjd(2011, 05, 30, 10,03);
  date2mjd(2011, 05, 30,10, 08) date2mjd(2011, 05, 30, 10,13);
  date2mjd(2011, 05, 30,10, 32) date2mjd(2011, 05, 30, 10,32) ; 
  date2mjd(2011, 05, 30,10, 50) date2mjd(2011, 05, 30, 10,57);
  date2mjd(2011, 05, 30,11, 03) date2mjd(2011, 05, 30, 11, 26);
  date2mjd(2011, 05, 30,11, 38) date2mjd(2011, 05, 30, 11, 48)];

ind = ones(size(d.array.frame.utc));
for m=1:size(badTimes,1)
  f = find(d.array.frame.utc>badTimes(m,1) & ...
      d.array.frame.utc<badTimes(m,2));
  ind(f) = 0;
end
d = framecut(d, logical(ind));
dc = framecut(d, d.array.frame.features>0);

ind = diff(dc.array.frame.utc)*24*60*60>2;
ind = [1; ind];

dc = framecut(dc, logical(ind));

ind = dc.antenna0.tracker.horiz_off(:,2)<-0.05;
dc = framecut(dc, ~ind);

save point_may30.mat dc;

%------------------------------------------------------------

d = read_arc('28-May-2011:05:11:17' , '28-May-2011:11:23:58');
dc = framecut(d, d.array.frame.features>0);

ind = diff(dc.array.frame.utc)*24*60*60>2;
ind = [1; ind];

dc = framecut(dc, logical(ind));

save point_may28_lst.mat dc;

%------------------------------------------------------------

d = read_arc('29-May-2011:05:01:30',  '29-May-2011:12:05:44');

badTimes = [...
  date2mjd(2011, 05, 29, 5, 28) date2mjd(2011, 05, 29, 5, 31);
  date2mjd(2011, 05, 29, 5, 34) date2mjd(2011, 05, 29, 5, 40);
  date2mjd(2011, 05, 29, 6, 17) date2mjd(2011, 05, 29, 6, 19);
  date2mjd(2011, 05, 29, 6, 29) date2mjd(2011, 05, 29, 6, 38);
  date2mjd(2011, 05, 29, 6, 43) date2mjd(2011, 05, 29, 6, 58);
  date2mjd(2011, 05, 29, 7, 09) date2mjd(2011, 05, 29, 7, 11);
  date2mjd(2011, 05, 29, 7, 16) date2mjd(2011, 05, 29, 7, 32);
  date2mjd(2011, 05, 29, 7, 37) date2mjd(2011, 05, 29, 7, 50);
  date2mjd(2011, 05, 29, 7, 56) date2mjd(2011, 05, 29,11, 50);]

ind = ones(size(d.array.frame.utc));
for m=1:size(badTimes,1)
  f = find(d.array.frame.utc>badTimes(m,1) & ...
      d.array.frame.utc<badTimes(m,2));
  ind(f) = 0;
end
d = framecut(d, logical(ind));
dc = framecut(d, d.array.frame.features>0);
ind = diff(dc.array.frame.utc)*24*60*60>2;
ind = [1; ind];

dc = framecut(dc, logical(ind));

save point_may29.mat dc;

%------------------------------------------------------------
regs={'array.frame.received'...
    'array.frame.utc double',...
    'array.frame.features',...
    'array.weather.utc double',...
    'array.weather.pressure double',...
    'array.weather.airTemperature double',...
    'array.weather.windSpeed[0] double',...
    'array.weather.windDirection[0] double',...
    'array.weather.status double',...
    'array.weather.relativeHumidity double',...
    'antenna0.thermal.ccHeaterCurrent double',...
    'antenna0.thermal.lsTemperatureSensors double',...
    'antenna0.thermal.ccTemperatureLoad double',...
    'antenna0.thermal.utc double',...
    'antenna0.thermal.dlpTemperatureSensors double',...
    'antenna0.tracker.lst double',...
    'antenna0.tracker.lacking double',...
    'antenna0.tracker.equat_geoc double',...
    'antenna0.tracker.horiz_topo double',...
    'antenna0.tracker.horiz_mount double',...
    'antenna0.tracker.flexure double',...
    'antenna0.tracker.horiz_off double',...
    'antenna0.tracker.tilts double',...
    'antenna0.tracker.fixedCollimation double',...
    'antenna0.tracker.encoder_off double',...
    'antenna0.tracker.sky_xy_off double',...
    'antenna0.tracker.source string',...
    'antenna0.tracker.scan_off double',...
    'antenna0.tracker.refraction double',...
    'antenna0.tracker.ut1utc double',...
    'antenna0.tracker.eqneqx double',...
    'antenna0.tracker.utc double',...
    'antenna0.tracker.time_diff double',...
    'antenna0.tracker.offSource double',...
    'antenna0.tracker.siteActual double',...
    'antenna0.tracker.siteFiducial double',...
    'antenna0.servo.utc double',...
    'antenna0.servo.fast_az_pos double',...
    'antenna0.servo.fast_el_pos double',...
    'antenna0.servo.fast_az_err double',...
    'antenna0.servo.fast_el_err double',...
    'antenna0.servo.actual_current_az1[0] double',...
    'antenna0.servo.actual_current_az2[0] double',...
    'antenna0.servo.actual_current_el1[0] double',...
    'antenna0.frame.utc double',...
    'antenna0.frame.received double',...
    'antenna0.receiver.flags double',...
    'antenna0.receiver.utc double',...
    'antenna0.receiver.data double',...
    'antenna0.receiver.diagnostics double',...
    'antenna0.receiver.drainCurrent double',...
    'antenna0.receiver.drainVoltage double',...
    'antenna0.receiver.gateVoltage double',...
    };
				
d = read_arc('28-oct-2011:02:15:23',  '28-Oct-2011:11:51:00', regs);

badTimes = [...
  date2mjd(2011, 10, 28, 4, 05) date2mjd(2011, 10, 28, 4, 08);
  date2mjd(2011, 10, 28, 7, 03) date2mjd(2011, 10, 28, 7, 07);
  date2mjd(2011, 10, 28, 7, 10) date2mjd(2011, 10, 28, 7, 14);
  date2mjd(2011, 10, 28, 8, 40) date2mjd(2011, 10, 28, 8, 49, 30);
  date2mjd(2011, 10, 28, 8, 55) date2mjd(2011, 10, 28, 9, 32);
  date2mjd(2011, 10, 28, 9, 46) date2mjd(2011, 10, 28,10, 36);
  date2mjd(2011, 10, 28,11, 08) date2mjd(2011, 10, 28,11, 11);
  date2mjd(2011, 10, 28,11, 31) date2mjd(2011, 10, 28,11, 36);
  date2mjd(2011, 10, 28,11, 39) date2mjd(2011, 10, 28,11, 48);]

ind = ones(size(d.array.frame.utc));
for m=1:size(badTimes,1)
  f = find(d.array.frame.utc>badTimes(m,1) & ...
      d.array.frame.utc<badTimes(m,2));
  ind(f) = 0;
end
d = framecut(d, logical(ind));
dc = framecut(d, d.array.frame.features>0);
ind = diff(dc.array.frame.utc)*24*60*60>2;
ind = [1; ind];

dc = framecut(dc, logical(ind));

dc1 = framecut(dc, dc.array.frame.features>0);

save point_20111028.mat dc;


d = read_arc('09-Feb-2012:02:51:45', '09-Feb-2012:12:32:00', regs);

badTimes = [...
  date2mjd(2012, 02, 09, 5, 43) date2mjd(2012, 02, 09, 5, 47);
  date2mjd(2012, 02, 09, 6, 11) date2mjd(2012, 02, 09, 6, 23);
  date2mjd(2012, 02, 09, 6, 31) date2mjd(2012, 02, 09, 6, 39);
  date2mjd(2012, 02, 09, 6, 42) date2mjd(2012, 02, 09, 6, 45, 30);
  date2mjd(2012, 02, 09, 6, 50) date2mjd(2012, 02, 09, 7, 03);
  date2mjd(2012, 02, 09, 9, 15) date2mjd(2012, 02, 09, 9, 19);
  date2mjd(2012, 02, 09, 9, 25) date2mjd(2012, 02, 09, 9, 29);
  date2mjd(2012, 02, 09, 9, 32) date2mjd(2012, 02, 09, 9, 41);
  date2mjd(2012, 02, 09,10, 32) date2mjd(2012, 02, 09,10, 37);
  date2mjd(2012, 02, 09,12, 24) date2mjd(2012, 02, 09,12, 29);]

ind = ones(size(d.array.frame.utc));
for m=1:size(badTimes,1)
  f = find(d.array.frame.utc>badTimes(m,1) & ...
      d.array.frame.utc<badTimes(m,2));
  ind(f) = 0;
end
d = framecut(d, logical(ind));
dc = framecut(d, d.array.frame.features>0);
ind = diff(dc.array.frame.utc)*24*60*60>2;
ind = [1; ind];

dc = framecut(dc, logical(ind));

dc1 = framecut(dc, dc.array.frame.features>0);

save point_20120209.mat dc1;
	  