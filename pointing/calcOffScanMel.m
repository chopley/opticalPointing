function [obs ide off out] = calcOffScan(d, plotparams, field, maindir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  function [obs ide off] = calcOffScan(d, [plotparams, field, maindir])
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
  doFit = 1;
  
  % cut the data for each scan
  ind = zeros(size(d.array.frame.features));
  ind(s(m):e(m)) = 1;
  ind = logical(ind);
  dc  = framecut(d, ind);
  aa  = unique(dc.antenna0.tracker.source);
  sourceName{m} = aa{1};
  timeVal(m)    = dc.array.frame.utc(1);
  
  % values for the pointing model
  idevals = dc.antenna0.tracker.horiz_topo(1,:);
  [obsvals(1) obsvals(2)] = pointing_model(om, idevals(1), idevals(2));
  ide.az(m) = idevals(1);
  ide.el(m) = idevals(2);
  obs.az(m) = obsvals(1);
  obs.el(m) = obsvals(2);
  
  % This one fits for 2D, so we keep the az/el scans and fit together.  

  % get the offsets from the control system.
  azApp = interp1(dc.antenna0.tracker.utc, ...
      dc.antenna0.tracker.horiz_topo(:,1), dc.antenna0.receiver.utc);
  azOffSave = azApp - 180/pi*unwrap(dc.antenna0.servo.apparent(:,1)*pi/180);
  azOffSave = -azOffSave;
  
  elApp = interp1(dc.antenna0.tracker.utc, ...
      dc.antenna0.tracker.horiz_topo(:,2), dc.antenna0.receiver.utc);
  elOffSave = elApp - dc.antenna0.servo.apparent(:,2);
  elOffSave = -elOffSave;
  
  % put it in sky angles.
  azOffSave = azOffSave.*cosd(obs.el(m));
  
  % assuming the slowest we would ever scan is at 0.2 degrees per second, a
  % az scan happens when the az speed is mostly constant, increasing, and
  % greater than 0.2 degrees/s, which is 
  azInd = (deriv(dc.antenna0.servo.az))>0.2/100;
  [startaz endaz] = findStartStop(azInd);
  azInd = zeros(size(azInd));
  if(~isempty(startaz))
    f = find( (endaz - startaz) == max(endaz - startaz));
    azInd(startaz(f):endaz(f)) = 1;
  else
    doFit = 0;
  end
  azInd = logical(azInd);

  % same goes for elevation
  elInd = deriv(dc.antenna0.servo.el)>0.2/100;
  [startel endel] = findStartStop(elInd);
  elInd = zeros(size(elInd));
  if(~isempty(startel))
    % no good el scan
    f = find( (endel-startel) == max(endel - startel));
    elInd(startel(f):endel(f)) = 1;
  else
    doFit = 0;
  end
  elInd = logical(elInd);

  % combine the two that we're interested in
  fitInd   = azInd | elInd;
  fitFlags = dc.flags.fast(fitInd, [1 3]);
  azFlags  = dc.flags.fast(azInd, [1 3]);
  elFlags  = dc.flags.fast(elInd, [1 3]);
  
  % setup mel's variables
  timeAz = dc.antenna0.receiver.utc(azInd);
  if(~isempty(timeAz))
    timeAz = (timeAz - timeAz(1)).*24*60;
  end
  skyAz  = azOffSave(azInd);
  newDataAz = dc.antenna0.receiver.data(azInd, intIndex);
  
  timeEl = dc.antenna0.receiver.utc(elInd);
  if(~isempty(timeEl))
    timeEl = (timeEl - timeEl(1)).*24*60;
  end
  skyEl  = elOffSave(elInd);
  newDataEl = dc.antenna0.receiver.data(elInd, intIndex);
  
  sigmaFit = 0.73/(2*sqrt(2*log(2)));

  % we should check we actually have a scan.
  if(length(find(elInd)) < 1500  | length(find(azInd)) < 1500)
    doFit = 0;
  end

  % check how much of the data are flagged
  elflag1  = length(find(elFlags(:,1)))./length(find(elInd));
  elflag2  = length(find(elFlags(:,2)))./length(find(elInd));
  azflag1  = length(find(azFlags(:,1)))./length(find(azInd));
  azflag2  = length(find(azFlags(:,2)))./length(find(azInd));
  
  if(elflag1 > 0.7 & elflag2 > 0.7)
    doFit = 0;
  end
  if(azflag1 > 0.7 & azflag2 > 0.7)
    doFit = 0;
  end
  
  
  if(doFit)
    % Calculate the baseline
    offEl = remove_background(timeEl, newDataEl, 'StepSize', 0.1, 'WindowSize', 1);
    f1 = find(offEl(:,1) == max(offEl(:,1)));
    f2 = find(offEl(:,2) == max(offEl(:,2)));
    elOff = mean(skyEl([f1 f2]));
    straightEl = newDataEl - offEl;
    % get the slope of the background signal and other parameters
    Pel(1,:) = polyfit([skyEl(1) last(skyEl)], [straightEl(1,1) ...
	  straightEl(length(timeEl),1)], 1);
    Pel(2,:) = polyfit([skyEl(1) last(skyEl)], [straightEl(1,2) ...
	  straightEl(length(timeEl),2)], 1);
    gradEl = Pel(:,1)';
    offsetEl = min(newDataEl);
    f = find(abs(skyEl)==min(abs(skyEl)));
    peakEl = newDataEl(f,:);
  
    
    
    offAz = remove_background(timeAz, newDataAz, 'StepSize', 0.1, 'WindowSize', 1);
    f1 = find(offAz(:,1) == max(offAz(:,1)));
    f2 = find(offAz(:,2) == max(offAz(:,2)));
    az_off = mean(skyAz([f1 f2]));
    straightAz = newDataAz - offAz;
    % get the slope of the background signal
    Paz(1,:) = polyfit([skyAz(1) last(skyAz)], [straightAz(1,1) ...
	  straightAz(length(timeAz),1)], 1);
    Paz(2,:) = polyfit([skyAz(1) last(skyAz)], [straightAz(1,2) ...
	  straightAz(length(timeAz),2)], 1);
    gradAz = Paz(:,1)';
    offsetAz = min(newDataAz);
    f = find(abs(skyAz)==min(abs(skyAz)));
    peakAz = newDataAz(f,:);
    azOff  = skyAz(f);
    
    % set up the parameters for the fit
    peakVal = mean([peakAz; peakEl]);
    
    % beta has two rows, one for each channel
    Beta = [peakVal', [azOff; azOff], [sigmaFit; sigmaFit], [elOff;elOff], [sigmaFit; sigmaFit], gradAz', gradEl', ...
	  offsetAz', offsetEl'];
    Beta(isnan(Beta)) = 0;
    
    % first let's get our data ready
    skylocs = [skyAz; skyEl];
    lenAz   = length(skyAz);
    dataFit = [newDataAz; newDataEl];
    
    % we don't want sidelobes in the fit
    indLobes = abs(skylocs)>0.5 & abs(skylocs)<1.7;
    % if it's CygA, we need to check for CygX
    indCygx = zeros(size(indLobes));
    if(sourceCorrespondance(unique(dc.antenna0.tracker.source))==1)
      % check the mean on either side of the peak have similar values
      azDat = newDataAz;
      azDat(azFlags) = nan;
      f1       = find(azDat(:,1) == max(azDat(:,1)));
      if(min(skyAz)< (skyAz(f1)-1.3))
	ind_no_lobe = skyAz < (skyAz(f1) - 1.1);
      else
	ind_no_lobe = skyAz < skyAz(f1);
      end
      m1_first = sqrt(nanvar(azDat(ind_no_lobe,:)));
      if(max(skyAz)> (skyAz(f1)+ 1.3))
	ind_no_lobe = skyAz > (skyAz(f1) + 1.1);
      else
	ind_no_lobe = skyAz > skyAz(f1);
      end
      m1_sec   = sqrt(nanvar(azDat(ind_no_lobe,:)));
      m1_diff  = ((m1_first-m1_sec)./m1_first)*100;
      
      elDat = newDataEl;
      elDat(elFlags) = nan;
      f2       = find(elDat(:,1) == max(elDat(:,1)));
      if(min(skyEl)< (skyEl(f2)-1.3))
	ind_no_lobe = skyEl < (skyEl(f2) - 1.1);
      else
	ind_no_lobe = skyEl < skyEl(f2);
      end
      m2_first = sqrt(nanvar(elDat(ind_no_lobe,:)));
      if(max(skyEl) > (skyEl(f2) + 1.3))
	ind_no_lobe = skyEl > (skyEl(f2) + 1.1);
      else
	ind_no_lobe = skyEl > skyEl(f2);
      end
      m2_sec   = sqrt(nanvar(elDat(ind_no_lobe,:)));
      m2_diff  = ((m2_first-m2_sec)./m2_first)*100;      
      
      indCygaz = zeros(size(skyAz));
      indCygel = zeros(size(skyEl));
      if(any(abs(m1_diff) > 35))
	% cygnux X is present in your az scan
	if(mean(m1_diff) > 0)
	  % it's present in first half
	  % don't get rid of the lobe itself.
	  ff = find(skyAz > (skyAz(f1) - 0.5), 1); 
	  indCygaz(1:ff) = 1;
	else
	  ff = find(skyAz > (skyAz(f1) + 0.5),1);
	  indCygaz(ff:lenAz) = 1;
	end
      end
      if(any(abs(m2_diff) > 35))
	% Cygnus X is present in your el scan
	if(mean(m2_diff) > 0)
	  % present in first half
	  ff = find(skyEl > (skyEl(f2)-0.5),1);
	  indCygel(1:ff) = 1;
	else
	  ff = find(skyEl > (skyEl(f2)+0.5),1);
	  indCygel(ff:length(skyEl)) = 1;
	end
      end
      indCygx = [indCygaz; indCygel];
    else
      indCygx = zeros(size(indLobes));
    end
  end

  % Mel's Method 
  for mm=1:2
    if(doFit)
      % do not select flagged data
      indNoFit = indCygx | indLobes | fitFlags(:,mm);
      indFit   = abs(skylocs)<5 & ~indNoFit;
    else
      indFit = 0;
    end
    if(length(find(indFit)) > 1000)
      % we fit if less than half the data are flagged
      thisData = dataFit(:,mm);
%            [betanew(mm,:), resid1, J1, cov1, mse] = nlinfit([lenAz; skylocs; indFit], thisData(indFit), ...
%      	  @gauss2D_fit2, Beta(mm,:));
%            errbeta(mm,:) = sqrt(abs(diag(cov1)));
%      LB = [0 -3 0.1 -3 0.1 -inf -inf -inf -inf];
%      UB = [inf 3 2 3 2 inf inf inf inf];
%      
%      [betanew(mm,:), resnorm, resid1, exitflag, output, lambda, J1] = ...
%	  lsqcurvefit('gauss2D_fit2', Beta(mm,:), [lenAz; skylocs; indFit], ...
%	  thisData(indFit), LB, UB);
%      errbeta = ones(size(betanew));
      display('calculating scan off');
      % let's try it with matmin because it gives us errors on the fits!!!!
      [betanew(mm,:), errbeta(mm,:), gof, stat] = matmin('chisq', Beta(mm,:), [], ...
	  'gauss2D_fit2', thisData(indFit), ones(size(thisData(indFit))), ...
	  [lenAz; skylocs; indFit]);

      fitModel = gauss2D_fit2(betanew(mm,:), [lenAz; skylocs; indFit]);
      resid1   = thisData(indFit) - fitModel;
      % calculate our models
      azScanModel = gauss2D_draw(betanew(mm,:), skyAz, zeros(size(skyAz)));
      elScanModel = gauss2D_draw(betanew(mm,:), zeros(size(skyEl)), skyEl);
      % flag if it's bad
      % we'll flag on three conditions:  
      % 1.  the x/y offset is more than 1 degree
      % 2.  the residuals on the fit data are high
      % 3.  the width of the gaussian doesn't match our beam to within 20%
      thisPointBad = 0;
      badString = 'GOOD';
      if(any(abs(betanew(mm,[2 4])) >1))
	display('Flagging point');
	display('Offset too great');
	thisPointBad = 1;
      end
      residPercent = mean(abs(resid1))./min(thisData)*100;
      if(mean(abs(resid1)) > 0.3)
	display('Flagging point');
	display('Residuals too high');
	mean(abs(resid1))
	thisPointBad = 1;
      end
      beamRatio = mean(betanew(mm, [3 5]))/sigmaFit;
      if( beamRatio < 0.3 | beamRatio > 2)
	display('Flagging point');
	display('Width does not match expected beam');
	display(sprintf('beam ratio is: %3.2f', beamRatio));
%	keyboard;
	thisPointBad = 1;
      end
      
      % select the information we want to keep
      indAzFit     = indFit(1:lenAz);
      indElFit     = indFit(lenAz+1:length(indFit));
      thisDataAz   = newDataAz(:,mm);
      thisDataEl   = newDataEl(:,mm);
      
      %  we need the az and el offset (az as in real azimuth, not sky angle)
      %  need the x and y widths 
      %  need the peak value, baseline (min) value and rms of residuals
      xwidth(m,mm) = betanew(mm,3);
      ywidth(m,mm) = betanew(mm,5);
      xoff(m,mm)   = betanew(mm,2);
      yoff(m,mm)   = betanew(mm,4);
      peakheight(m,mm) = betanew(mm,1);
      errpeak(m,mm)= errbeta(mm,1);
      errxoff(m,mm)= errbeta(mm,2);
      erryoff(m,mm)= errbeta(mm,4);
      eloff(m,mm)  = betanew(mm,4);
      azoff(m,mm)  = betanew(mm,2)./cosd(obs.el(m));
      faz          = find( abs(skyAz - xoff(m,mm)) == min(abs(skyAz - xoff(m,mm))));
      azpeak(m,mm) = thisDataAz(faz);
      thisDataAz(~indAzFit) = nan;
      azmin(m,mm)  = min(thisDataAz);
      azrms(m,mm)  = rms(thisDataAz - azScanModel);
      fel          = find( abs(skyEl - yoff(m,mm)) == min(abs(skyEl - yoff(m,mm))));
      elpeak(m,mm) = thisDataEl(fel);
      thisDataEl(~indElFit) = nan;
      elmin(m,mm)  = min(thisDataEl);
      elrms(m,mm)  = rms(thisDataEl - elScanModel);
      
      % next we plot the point
      display(sprintf('Plotting Point %d of %d', m, length(s)));
      clf
      skywrap = unwrap(skyAz);
      subplot(2,1,1)
      plot(skywrap, newDataAz(:,mm));
      xlabel('Sky Az Offset (degrees)');
      ylabel('Intensity');
      if(1)
	hold on; plot(skywrap, azScanModel, 'r'); hold off
	axis([min([skywrap; skyEl]) max([skyAz; skyEl]) min(newDataAz(:,mm)) max(newDataAz(:,mm))]);
      end
      subplot(2,1,2)
      plot(skyEl, newDataEl(:,mm));
      xlabel('El Offset (degrees)');
      ylabel('Intensity');
      if(1)
	hold on; plot(skyEl, elScanModel, 'r'); hold off;
	axis([min([skyAz; skyEl]) max([skyAz; skyEl]) min(newDataEl(:,mm)) max(newDataEl(:,mm))]);
      end
      eval(sprintf('title(''Intensity Channel # %d'');', mm));
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
	badPoint(m,mm) = thisPointBad;
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
      
    else  % if we're not fitting because we don't have enough data.
      display('The majority of the data in this cross was flagged');
      xwidth(m,mm) = nan;
      ywidth(m,mm) = nan;
      xoff(m,mm)   = nan;
      yoff(m,mm)   = nan;
      peakheight(m,mm) = nan;
      errpeak(m,mm)= nan;
      errxoff(m,mm)= nan;
      erryoff(m,mm)= nan;
      eloff(m,mm)  = nan;
      azoff(m,mm)  = nan;
      azmin(m,mm)  = nan;
      azpeak(m,mm) = nan;
      azrms(m,mm)  = nan;
      elmin(m,mm)  = nan;
      elpeak(m,mm) = nan;
      elrms(m,mm)  = nan;
      badPoint(m,mm) = 1;
    end  % with check of whether we should fit
  
  end   % of loop over intensity channels
  pause(0.1);
  
end  % of loop over source scans

% first we assign a weight with the source.
azSig = (azpeak - azmin)./azrms;
azErr = 1./azSig;

elSig = (elpeak - elmin)./elrms;
elErr = 1./elSig;

% throw out values where the error is large or negative.
elBad = elErr < 0 | elErr > 1;
azBad = azErr < 0 | azErr > 1;

indBad = logical(badPoint);
indBad(elBad) = 1;
indBad(azBad) = 1;

eloff(indBad) = nan;
azoff(indBad) = nan;

% if hte values for both channels aren't within 10 arcmin, toss them.
dEl = abs(diff(eloff, [], 2));
dAz = abs(diff(azoff, [], 2).*cosd(obs.el)');
f = find(dEl>15/60 | dAz>15/60);
indBad(f,:) = 1;

% save all the data
out.time = timeVal;
out.name = sourceName;
out.az   = ide.az;
out.el   = ide.el;

badPoints = indBad;

indBad = sum(indBad,2)==2;

display(sprintf('Throwing out %d points', length(find(indBad))));
azoff(badPoints) = nan;
eloff(badPoints) = nan;

[off.az  off.azErr] = vwstat(azoff, azErr, 2);
[off.el  off.elErr] = vwstat(eloff, elErr, 2);

display('Done with all crosses');

% add final offests to the observed az/el
obs.az = obs.az + off.az';
obs.el = obs.el + off.el';

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


% all the outputs we might want.
out.peakheight = peakheight;
out.xwidth = xwidth;  
out.ywidth = ywidth; 
out.xoff   = xoff;   
out.yoff   = yoff;
out.errpeak= errpeak;
out.errxoff= errxoff;
out.erryoff= erryoff;
out.eloff  = eloff;  
out.azoff  = azoff;  
out.azmin  = azmin;  
out.azpeak = azpeak; 
out.azrms  = azrms;  
out.elmin  = elmin;  
out.elpeak = elpeak; 
out.elrms  = elrms;  
out.azsig  = azSig;
out.elsig  = elSig;
out.badpt  = badPoints;


return;

