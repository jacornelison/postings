function cut=get_default_round1_cuts(expt)
% cut=get_default_round1_cuts(expt)
%
% get a set of vanilla cuts
% - all are listed here - although some are set to "inf" etc

if(~exist('expt','var'))
  expt=get_experiment_name;
end

cut.fb_std_p0=[0,inf];
% cut data which shows signs of FPU temp disturbance
switch expt
 case 'bicep2'
  cut.fp_std_p0_darks=[0,5];
 case 'keck'
  cut.fb_std_p0_darks=[0,10];
end
cut.fb_std_p3=[0,inf];
cut.fb_std_sd_p0=[0,inf];
cut.fb_std_sd_p3=[0,inf];
cut.fb_std_uncal=[0,inf];

% throw out any halfscan with a nan (glitch/step detected)
cut.fb_nancount=1;

% Cut row/col where fluc jump has been detected
% this is max allowed value so 0.5 turns the cut on and 1.5 turns it off
cut.is_fj_row=1.5;
cut.is_fj_col=1.5;

% don't change this - these cuts are needed
cut.syncsampnum_diff1=0;
cut.syncsampnum_diff2=0;

% Passfrac for BICEP2/Keck
% these may not be the best vals...
switch expt
  case {'bicep2','keck'}
    cut.passfrac_col=0.7;
    cut.passfrac_chan=0.7;
  case 'bicep3'
    cut.passfrac_col=0.7;
    cut.passfrac_chan=0.7;
end

return
