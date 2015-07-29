function date=utc2date(utc)
% function date=utc2date(utc)
%
% give a sensible date from the utc, which is MJD

% subtract off Jan 1, 2000
utc=utc-51543;
date.year=floor(utc/365.24)+2000;
day=utc-floor(utc/365.24)*365.24;

md=[31 28 31 30 31 30 31 31 30 31 30 ...
      31];
months={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

for m=1:12
  day=day-md(m);
  if day<0
    break
  end
end
date.month=char(months(m));

date.day=floor(day)+md(m);

subday=day-floor(day);
date.hour=floor(subday*24);
subhour=24*subday-floor(24*subday);
date.min=floor(subday*60);
submin=60*subhour-floor(60*subday);
date.sec=floor(subday*60);

return

  
