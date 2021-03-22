% example pwelch spectra

clear Fm Sm
Fs = 1; % sampling rate; here: 1 day
ws = 300; % window's length; here: 300 days
win = ws*Fs; noverlap = win/2; nfft = win; % input for pwelch

[Sm,Fm] = pwelch(detrend(serie),win,noverlap,nfft,Fs); % serie: your time series

figure(1);clf % simply test period vs spectral estimates
plot(1./Fm,Sm,'r')

% to compute dof:
npts =length(serie);
nbands = 1; % if you are not doing any frequency-band averaging
ndof=round(8./3.*npts./nfft.*nbands); %% MUST BE INTEGER!
conf=chi2conf(0.95,ndof(1));

% plot variance preserving spectra (semilog) ----

pagePortrait(1,16); clf
position = [xoff he1 wi he];
axes('units','inches','position',position);

% plot spectra: ****************
semilogx(Fm,Fm.*Sm,'-','color','k','Linewidth',lw);  % Fm.*Sm % Variance preserving!!

% to plot shading for confidence interval: *****
Fm(1) = []; % to plot spectra and get rid of freq =0
Sm(1) = [];
shx=[Fm; flipud(Fm)];
shy=[Sm*conf(1).*Fm ; flipud(Sm*conf(2).*Fm)];
hf=patch(shx,shy,[0.85 0.85 0.85]);
set(hf,'edgecolor','none') 
%**********************************************

set(gca,'xscale','log')
grid on
set(gca,'box','on')
set(gca,'tickdir','out')
set(gca,'xlim',[1e-3 1e0])
set(gca,'ylim',[0 20],'ytick',[0:2:20])
ylabel('Transport Variance [Sv^{2}]') % variance preserving! Units: units of variance
xlabel('Frequency [cycles per day]') % can add period labels on top


