function [obs ide off out] = calcOffScan(d, plotparams, field, maindir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  function [obs ide off out] = calcOffScan(d, [plotparams, field, maindir])
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==1
  plotparams.plot = 1;
  plotparams.interactive = 1;
  plotparams.save = 0;
  field = [];
  maindir = [];
end

setPlotDisplay(plotparams.plot);

% load up the current model
om(:,1)=d.antenna0.tracker.flexure(:,1);
om(:,2)=d.antenna0.tracker.flexure(:,2);
om(:,3)=d.antenna0.tracker.tilts(:,1);
om(:,4)=d.antenna0.tracker.tilts(:,2);
om(:,5)=d.antenna0.tracker.tilts(:,3);
om(:,6)=d.antenna0.tracker.fixedCollimation(:,1);
om(:,7)=d.antenna0.tracker.fixedCollimation(:,2);
om(:,8)=d.antenna0.tracker.encoder_off(:,1);
om(:,9)=d.antenna0.tracker.encoder_off(:,2);
om = mean(om);

% I. split them up into respective scans
% Ia.  find the start and stop index of each scan.
[s e] = findStartStop(d.index.radio_point_scan.slow);
betaAzAll = nan(length(s), 2, 3);
betaElAll = nan(length(s), 2, 3);


if(isempty(e) | isempty(s))
  obs.az = nan;
  obs.el = nan;
  ide.az = nan;
  ide.el = nan;
  off.az = nan;
  off.el = nan;
  goodPoint = [0 0];
  display('No crosses');
  s = [];
  return;
  
end

if(e(1)<s(1))
  obs.az = nan;
  obs.el = nan;
  ide.az = nan;
  ide.el = nan;
  off.az = nan;
  off.el = nan;
  goodPoint = [0 0];
  display('No crosses');
  s = [];
  return;
end


% determine the index of the intensity channel
switch (size(d.antenna0.receiver.data,2))
  case 6
    intIndex = [1 6];

  case 8
    intIndex = [1 8];
    
  case 10
    intIndex = [1 9];
end

plotIndex = 0;

% let's check if all the powers are negative.  if they are, we do an
% absolute value. 
if(all((all(d.antenna0.receiver.data(:,intIndex)<0))))
  display('WARNING:  Your data are negative');
  display('Taking the absolute value only for fitting offsets');
  display('This will ONLY change the sign of the data in the pointing plots');
  display('All other analysis will NOT be affected');
  d.antenna0.receiver.data = abs(d.antenna0.receiver.data);
end

for m=1:length(s)
  % cut the data for each scan
  ind = zeros(size(d.array.frame.features));
  ind(s(m):e(m)) = 1;
  ind = logical(ind);
  dc  = framecut(d, ind);

  
  
  % find the observation az/el
%  indOnSource = dc.antenna0.tracker.offSource==0;
%  indAzOff    = dc.antenna0.tracker.horiz_off(:,1) == 0 & ...
%      dc.antenna0.tracker.scan_off(:,1) ==0;
%  indElOff    = dc.antenna0.tracker.horiz_off(:,2) == 0 & ...
%      dc.antenna0.tracker.scan_off(:,2) ==0;
%  f = find(indOnSource & indAzOff & indElOff);
%  if(length(f)==0)
%    f = find(indOnSource);
%  end
%  if(length(f)>1)
%    obsvals = (dc.antenna0.tracker.horiz_mount(f(1),:)) + dc.antenna0.tracker.encoder_off(f(1),:);      
%    idevals = (dc.antenna0.tracker.horiz_topo(f(1),:));
%    
%  elseif(length(f) == 1)
%    obsvals = dc.antenna0.tracker.horiz_mount(f,:) + dc.antenna0.tracker.encoder_off(f,:);    
%    idevals = (dc.antenna0.tracker.horiz_topo(f,:));
%  end


  idevals = dc.antenna0.tracker.horiz_topo(1,:);
  [obsvals(1) obsvals(2)] = pointing_model(om, idevals(1), idevals(2));
  ide.az(m) = idevals(1);
  ide.el(m) = idevals(2);
  obs.az(m) = obsvals(1);
  obs.el(m) = obsvals(2);
  
  
  
  
  
  % next we split it into the azimuth/elevation scan
  % for the elevation scan, we fit a line and remove the slope (effectively
  % a sky dip), and fit a gaussian to the remainder
  % for the azimuth scan, we just fit a gaussian to it.

  % get the offsets from the control system.
  azApp = interp1(dc.antenna0.tracker.utc, ...
      dc.antenna0.tracker.horiz_topo(:,1), dc.antenna0.receiver.utc);
  azOffSave = azApp - (180/pi)*(unwrap(pi/180*dc.antenna0.servo.apparent(:,1)));
  
  %  [aa, ee ] =calcAzElCs(dc, 1);
%  ao = aa - dc.antenna0.servo.az;
%  azOffSave = ao;
%  azOffSave(azOffSave<-180) = azOffSave(azOffSave<-180) + 360;
%  azOffSave(azOffSave>180) = azOffSave(azOffSave>180) - 360;
   azOffSave = -azOffSave;

  elApp = interp1(dc.antenna0.tracker.utc, ...
      dc.antenna0.tracker.horiz_topo(:,2), dc.antenna0.receiver.utc);
  elOffSave = elApp - dc.antenna0.servo.apparent(:,2);
  %elOffSave = ee - dc.antenna0.servo.el;
  elOffSave = -elOffSave;
  
  % assuming the slowest we would ever scan is at 0.2 degrees per second, a
  % az scan happens when the az speed is mostly constant, increasing, and
  % greater than 0.2 degrees/s, which is 
  azInd = (deriv(dc.antenna0.servo.az))>0.1/100;
  [startaz endaz] = findStartStop(azInd);
  f = find( (endaz - startaz) == max(endaz - startaz));
  azInd = zeros(size(azInd));
  azInd(startaz(f):endaz(f)) = 1;
  azInd = logical(azInd);
  azIs  = dc.antenna0.receiver.data(azInd,intIndex);
  azIs(:,1) = smooth(azIs(:,1));
  azIs(:,2) = smooth(azIs(:,2));

  azFlags       = dc.flags.fast(azInd,[1 3]);
  azIs(azFlags) = nan;
  
  
  % if amplitude is negative, no peak in data
  
  % same goes for elevation
  elInd = deriv(dc.antenna0.servo.el)>0.1/100;
  [startel endel] = findStartStop(elInd);
  f = find( (endel-startel) == max(endel - startel));
  elInd = zeros(size(elInd));
  elInd(startel(f):endel(f)) = 1;
  elInd = logical(elInd);
  elIs  = dc.antenna0.receiver.data(elInd,intIndex);
  elIs(:,1) = smooth(elIs(:,1));
  elIs(:,2) = smooth(elIs(:,2));
  
  elFlags       = dc.flags.fast(elInd, [1 3]);
  elIs(elFlags) = nan;
  
  % if most of the scan is flagged, don't bother fitting.
  azFlagPer = length(find(azFlags))./length(azIs(:));
  elFlagPer = length(find(elFlags))./length(elIs(:));
  if(azFlagPer < 0.15 & elFlagPer < 0.15)
    doFit = 1;
  else
    doFit = 0;
  end

  aa = unique(dc.antenna0.tracker.source);
  sourceName{m} = aa{1};
  timeVal(m)    = dc.array.frame.utc(1);
  
  if(doFit)
    testVals = [120:30:800];
    for mm=1:2
      
      if(~isempty(find(azInd)))
	% fit a gaussian about the max in both cases
	azOff = azOffSave;
	azI = azIs(:,mm);
	f = find(azI == max(azI));
	azOff = azOff(azInd);
	
	stopAz = 0;
	
	azIndex = 1;
	while(stopAz==0)
	  azFit = find(azInd);
	  if(f<=testVals(azIndex) | (f+testVals(azIndex))>=length(azI))
	    badAzPoint = 1;
	    azFit = [];
	    xoff = nan;
	    azIndex = azIndex+1;
	    betaAz(1:3) = nan;
	    ampVal = nan;
	  else
	    azFit = azOff((f-testVals(azIndex):f+testVals(azIndex)));
	    azIFit = azI(f-testVals(azIndex):f+testVals(azIndex));
	    indgood= zeros(size(azI));
	    indgood(f-testVals(azIndex):f+testVals(azIndex)) = 1;
	    
	    
%	    display(sprintf('testing value: %d', testVals(azIndex)));
	    
	    OPTIONS = optimset('Display', 'off', 'MaxIter', 200000, ...
		'MaxFunEvals', 20000, 'TolX', 0.00005, 'TolFun', 0.00005);
	    LB = [0 azFit(testVals(azIndex))-0.2  -10 -inf -inf];
	    UB = [max(azIFit) azFit(testVals(azIndex))+0.2 10 inf inf];
	    x0 = [1 azFit(testVals(azIndex)-1) 1 1 1];
	    [betaAz resnorm residual exitflag] = lsqnonlin(@(x) mchisq2(azIFit, gaussfit(x, azFit)), x0, LB, UB, OPTIONS);
	    
	    amplitude = betaAz(1);
	    xoff      = betaAz(2);
	    sigma     = betaAz(3); 
	    xoffIndex = find(abs(azFit - xoff) == min(abs(azFit - xoff)));
	    ampVal    = azIFit(xoffIndex);
	    
	    % other terms remove a baseline
	    if (amplitude<0 | abs(xoff)>8 | resnorm>500)
	      badAzPoint = 1;
	    else
	      badAzPoint = 0;
	    end 
	    azFit = gaussfit(betaAz, azOff);
	    
	    if(abs(xoff)<4)
	      stopAz = 1;
	    else
	      azIndex = azIndex+1;
	    end
	  end
	  
	  if(azIndex>length(testVals))
	    stopAz = 1;
	  end
	end
	% save for weighting
	if(azIndex >= length(testVals))
	  azIndex = azIndex-1;
	end
	azWidth(m,mm) = testVals(azIndex);
	azMax(m,mm)   = ampVal;
	if(~badAzPoint)
	  azBase(m,mm)  = nanmedian(azI(~indgood));
	  azRms(m,mm)   = nanstd(azI(~indgood));
	else
	  azBase(m,mm)  = nan;
	  azRms(m,mm)   = nan;
	end
	betaAzAll(m,mm,:) = betaAz(1:3);
      
	% check one last time if it's a good point.
	sigDet  = (azMax(m,mm) - azBase(m,mm) )./azRms(m,mm);
	if(sigDet < 3 | isnan(sigDet))
	  badAzPoint = 1;
	  azBase(m,mm) = nan;
	  azRms(m,mm) = nan;
	end
      end
      
      
      if(~isempty(elInd))
	elOff = elOffSave;
	elI = elIs(:,mm);
	% first we remove the offset from slewing the atmosphere
	% from the sky dip code
	x = 1./(sind(dc.antenna0.servo.el(elInd)));
	y = abs(elI);
	x(isnan(y)) = [];
	y(isnan(y)) = [];
	[tau Tground] = linfit(x,y);
	% remove the effect
	elI = elI - Tground(1) - tau(1)./sind(dc.antenna0.servo.el(elInd));
	
	% fit a gaussian about the max 
	frange = elOff(elInd);
	indrange = abs(frange)<4;
	f = find(elI(indrange) == max(elI(indrange))) + find(indrange,1);
	elOff2 = elOff(elInd);
	
	stopEl = 0;
	elIndex = 1;
	while(stopEl==0)
	  elFit = find(elInd);
	  if(f<=testVals(elIndex) | (f+testVals(elIndex))>=length(elI))
	    badElPoint = 1;
	    elFit = [];
	    yoff = nan;
	    elIndex = elIndex+1;
	    betaEl = [nan nan nan];
	    ampVal = nan;
	  else
	    elFit = elOff2((f-testVals(elIndex):f+testVals(elIndex)));
	    elIFit = elI(f-testVals(elIndex):f+testVals(elIndex));
	    indgood= zeros(size(elI));
	    indgood(f-testVals(elIndex):f+testVals(elIndex)) = 1;
	    
%	    display(sprintf('el testing value: %d', testVals(elIndex)));
	    [betaEl] = nlinfit(elFit, elIFit, @gaussfit, [1 elFit(testVals(elIndex)) 1 1 1]);
	    amplitude = betaEl(1);
	    yoff      = betaEl(2);
	    sigma     = betaEl(3);
	    yoffIndex = find(abs(elFit - yoff) == min(abs(elFit - yoff)));
	    ampVal    = elIFit(yoffIndex);	    
	    
	    % other terms remove a baseline
	    if (amplitude<0 | abs(yoff)>8)
	      badElPoint = 1;
	    else
	      badElPoint = 0;
	    end 
	    elFit = gaussfit(betaEl, elOff2);
	    
	    if(abs(yoff)<1)
	      stopEl = 1;
	    else
	      elIndex = elIndex+1;
	    end
	  end
	  if(elIndex>length(testVals))
	    stopEl = 1;
	  end
	  
	end
	% save for weighting
	if(elIndex >= length(testVals))
	  elIndex = elIndex - 1;
	end
	elWidth(m,mm) = testVals(elIndex);
	elMax(m,mm)   = ampVal;
	if(~badElPoint)
	  elBase(m,mm)   = nanmedian(elI(~indgood));
	  elRms(m,mm)    = nanstd(elI(~indgood));
	else
	  elBase(m,mm)   = nan;
	  elRms(m,mm)    = nan;
	end
	betaElAll(m,mm,:) = betaEl(1:3);      
      end
      
%      display('here');
      
      display(sprintf('Point %d of %d', m, length(s)));
      
      
      clf
      %if(~any(isnan(betaAzAll(m,mm,:))))
      subplot(2,1,1)
      plot( (azOff)*cos(dc.antenna0.servo.el(1)*pi/180), azI');
      xlabel('Az Offset (degrees)');
      ylabel('Intensity');
      if(~isempty(azFit))
	hold on; plot(azOff*cos(dc.antenna0.servo.el(1)*pi/180), azFit, 'r'); hold off
	axis([min(azOff) max(azOff) min(azI) max(azI)]);
      end
      
      %end
      
      %if(~any(isnan(betaElAll(m,mm,:))))
      subplot(2,1,2)
      plot(elOff(elInd), elI);
      xlabel('El Offset (degrees)');
      ylabel('Intensity');
      if(~isempty(elFit))
	hold on; plot(elOff(elInd), elFit, 'r'); hold off
	axis([min(elOff) max(elOff) min(elI) max(elI)]);
      end
      %end

      eval(sprintf('gtitle(''scan on source %s:'', 0.96, 2);', ...
	  dc.antenna0.tracker.source{1}));
      
      if(plotparams.interactive)
	display(sprintf('Source Name: %s', dc.antenna0.tracker.source{1}));
	r = input('Keep This point? [y/n]: ', 's');
	
	if(strcmp(r, 'y') || strcmp(r, 'Y'))
	  badPoint(m,mm) = 0;
	  badString = 'GOOD';
	elseif(strcmp(r, 'n') ||  strcmp(r, 'N'))
	  badPoint(m,mm) = 1;
	  badString = 'BAD';
	elseif(strcmp(r, 'k'))
	  keyboard;
	end
      else
	badPoint(m,mm) = badAzPoint | badElPoint;
	if(badPoint(m,mm))
	  badString = 'BAD ';
	else
	  badString = 'GOOD';
	end
      end
      hold off;
      eval(sprintf('gtitle(''%s'', 0.96, 3);',...
	  badString));
      
      if(plotparams.save)
	plotIndex = plotIndex+1;
	if(isempty(maindir))
	  maindir = getMainDir(d, field);
	end
	dbclear if error
	set(gcf,'paperposition',[0 0 6.0 9.0])
	filename = sprintf('%s/pointing/fig%d.png', maindir, plotIndex);
	eval(sprintf('print -dpng -r70 %s', filename));
	dbstop if error
      end
    end    
    clear azIs;
    clear elIs;
    
    pause(0.1);
  else
    display('All data in this cross have been flagged');
    azWidth(m,1:2) = nan;
    azMax(m,1:2)   = nan;
    azBase(m,1:2)  = nan;
    azRms(m,1:2)   = nan;
    betaAzAll(m,1:2,1:3) = nan;
    elWidth(m,1:2) = nan;
    elMax(m,1:2)   = nan;
    elBase(m,1:2)  = nan;
    elRms(m,1:2)   = nan;
    betaElAll(m,1:2,1:3) = nan;
    badPoint(m,1:2) = 1;
  end

end


badPoint = logical(badPoint);
elOff = betaElAll(:,:,2);
azOff = betaAzAll(:,:,2);
elOff(badPoint) = nan;
azOff(badPoint) = nan;
indBad = badPoint;

% first we assign a weight with the source.
azSig = (azMax - azBase)./azRms;
azWidth = abs(betaAzAll(:,:,3)/2);
azErr = azWidth./azSig;


elSig = (elMax - elBase)./elRms;
elWidth = abs(betaElAll(:,:,3)/2);
elErr   = elWidth./elSig;


% throw out values where the error is large or negative.
elBad = elErr < 0 | elErr > 1;
azBad = azErr < 0 | azErr > 1;

elOff(elBad) = nan;
azOff(azBad) = nan;
elErr(elBad) = nan;
azErr(elBad) = nan;
indBad(elBad) = 1;
indBad(azBad) = 1;


% if hte values for both channels aren't within 10 arcmin, toss them.
dEl = abs(diff(elOff, [], 2));
dAz = abs(diff(azOff, [], 2).*cos(obs.el*pi/180)');

f = find(dEl>20/60 | dAz>20/60);
elOff(f,:) = nan;
azOff(f,:) = nan;
elErr(f,:) = nan;
azErr(f,:) = nan;
indBad(f,:) = 1;

badPoint = indBad;

indBad = sum(indBad,2)==2;

display(sprintf('Throwing out %d points', length(find(indBad))));

[off.az  off.azErr] = vwstat(azOff, azErr, 2);
[off.el  off.elErr] = vwstat(elOff, elErr, 2);

display('Done with all crosses');

% add final offests to the observed az/el
obs.az = obs.az + off.az';
obs.el = obs.el + off.el';

obsel = obs.el;

out.time   = timeVal;
out.name   = sourceName;
out.az     = obs.az;
out.el     = obs.el;

% throw out the data that is bad
ind = isnan(obs.az);
obs.az(ind) = [];
obs.el(ind) = [];
obs.name    = sourceName(~ind);
obs.timeVal = timeVal(~ind);
obs.index   = find(~ind);
ide.az(ind) = [];
ide.el(ind) = [];
off.az(ind) = [];
off.el(ind) = [];
off.azErr(ind) = [];
off.elErr(ind) = [];


%azOff(ind,:) = [];
%azErr(ind,:) = [];


out.xwidth = nan(size(azOff));
out.ywidth = nan(size(azOff));
out.xoff   = azOff./sind([obsel' obsel']);
out.yoff   = elOff;
out.errpeak= nan(size(azErr));
out.errxoff= azErr./sind([obsel' obsel']);
out.erryoff= elErr;
out.eloff  = elOff;
out.azoff  = azOff;
out.azmin  = azBase;
out.azpeak = azMax;
out.azrms  = azRms;
out.elmin  = elBase;
out.elpeak = elMax;
out.elrms  = elRms;
out.azsig  = azSig;
out.elsig  = elSig;
out.badpt  = badPoint | isnan(azOff);


return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [s e] = findStartStop(features)
%
%% start indices are when the difference in the feature value increases, stop
%% is when it decreases
%s = find(diff(double(features))>0);
%
%% end is when it decreases
%e = find(diff(double(features))<0);
%
%% last one
%lastOneBad = 0;
%if(length(e)<length(s))
%  % the dimensions don't match
%  if(s(1)<e(1))
%    % if start is less than end, that means we're usually missing hte last
%    % endpoint
%    e(length(s)) = length(features);
%    lastOneBad = 1;
%  else
%    s(1) = [];
%  end
%elseif(length(s)<length(e))
%  s(1) = [];
%end
%s = s+1;
%
%goodPoint = zeros(size(s));
%for m=1:length(s)
%  thisFeat = unique(features(s(m):e(m)));
%  if(thisFeat == 129 | thisFeat ==128)
%    goodPoint(m) = 1;
%  end
%end  
%
%if(lastOneBad==1)
%  goodPoint(length(s)) = 0;
%end
%
%
%s = s(logical(goodPoint));
%e = e(logical(goodPoint));
%
%return;


function x = mchisq2(dataVals, modelVals)

x = dataVals - modelVals;

f = find(isnan(x));
x(f) = [];
f = find(isinf(x));
x(f) = [];

if(isempty(x))
  display('damn it all to hell')
  keyboard;
end

return;

	    

