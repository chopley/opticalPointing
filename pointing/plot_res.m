function [rmssa ind sa]=plot_res(model,ide,obs,plottitle,gof, junk)
% plot_res(model,ide,obs,plottitle,gof)
%
% Plot the model residuals.  gof is the goodness of fit for plots of
% fit models

d2r=pi/180;

% Find the model values
[mva.az,mva.el]=pointing_model(model,ide.az,ide.el);

% Find the az/el residuals
res.az=obs.az-mva.az;

ind=res.az>180; res.az(ind)=res.az(ind)-360;
ind=res.az<-180; res.az(ind)=res.az(ind)+360;
res.el=obs.el-mva.el;

% Find the space angles between observed and reference values
sa=spaceangle(obs.az,obs.el,mva.az,mva.el,'deg');
rmssa=rms(sa);
ind=sa>3*rmssa;

% Calc the cos(el) factor
ide.cosel=cos(ide.el*d2r);

setwinsize(gcf,600,600); clf


% Plot the sky positions hit
subplot(3,3,3)
polar(ide.az*d2r,90-ide.el,'+');
view(90,-90);
hold on; polar(ide.az(ind)*d2r,90-ide.el(ind),'ro'); hold off


subplot(3,3,1)
plot(ide.az,res.az,'+'); xlabel('az'); ylabel('az residual'); %ylim([-0.5 0.5]);
hold on;  plot(ide.az(ind),res.az(ind),'ro'); hold off
subplot(3,3,2);
plot(ide.el,res.az,'+'); xlabel('el'); ylabel('az residual'); %ylim([-0.5 0.5]);
hold on; plot(ide.el(ind),res.az(ind),'ro'); hold off

subplot(3,3,4)
plot(ide.az,res.az.*ide.cosel,'+'); xlabel('az'); ylabel('az residual * cos(el)'); %ylim([-0.5 0.5]);
hold on; plot(ide.az(ind),res.az(ind).*ide.cosel(ind),'ro'); hold off
subplot(3,3,5);
plot(ide.el,res.az.*ide.cosel,'+'); xlabel('el'); ylabel('az residual * cos(el)'); %ylim([-0.5 0.5]);
hold on; plot(ide.el(ind),res.az(ind).*ide.cosel(ind),'ro'); hold off

subplot(3,3,7)
plot(ide.az,res.el,'+'); xlabel('az'); ylabel('el residual'); %%ylim([-0.01,0.02]);
hold on; plot(ide.az(ind),res.el(ind),'ro'); hold off
subplot(3,3,8)
plot(ide.el,res.el,'+'); xlabel('el'); ylabel('el residual'); %%ylim([-0.01,0.02]);
hold on; plot(ide.el(ind),res.el(ind),'ro'); hold off

% Write the model parameters if available
if(exist('model'))
   subplot(3,3,6);
  modpar={'flex sin','flex cos','az tilt ha','az tilt lat',...
      'el tilt','collim x','collim y','az zero','el zero','gof','num pts'};
  disp_vals=[model,gof,sum(~isnan(obs.az))];
  for i=1:length(modpar)
    h=text(0.1,(10-i)*0.1,sprintf('%s = %f',modpar{i},disp_vals(i)));
    if(disp_vals(i)~=0)
      set(h,'Color','b');
    end
  end
  axis off
end

% Plot space angle histograms
subplot(3,3,9)
sa=sa*60;
m=max(sa)*1.000001;
%hfill(sa,20,0,m,[]);
%hfill(sa(~ind),20,0,m,[]);
xlabel('space angle error (arcmin)');
y=ylim; ylim([0,y(2)*1.1]);
axes('Position',get(gca,'Position'),'Xlim',[0,1],'Ylim',[0,1],'Visible','off');
text(0.4,0.9,sprintf('rms=%f',rms(sa(~ind))),'Color','b');
if(sum(ind)>0)
  text(0.4,0.8,sprintf('rms=%f',rms(sa)),'Color','r');
end

if(exist('plottitle'))
  gtitle(plottitle);
end

return
