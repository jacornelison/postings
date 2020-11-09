function simopt = get_default_simopt(simopt)
% simopt = get_default_simopt(simopt)

% Get default simopt parameters for reduc_makesim.  Called at the top of
% reduc_makesim.m
% does NOT assign value to simopt.sernum.

% ex. simopt.sig='normal';
%     simopt.sernum='12050011';
%     simopt=get_default_simopt(simopt);
%               OR
%     simopt=get_default_simopt;

if(~exist('simopt','var'))
  simopt=[];
end

% assign default value for each sim option
if(~isfield(simopt,'noise'))
  simopt.noise='data';
end
if(~isfield(simopt,'ntags'))
  simopt.ntags=1;
end
if(~isfield(simopt,'sig'))
  simopt.sig='normal';
end
if(~isfield(simopt,'sigmaptype'))
  simopt.sigmaptype='healmap';
end
if(~isfield(simopt,'sigmapfilename'))
  simopt.sigmapfilename='input_maps/camb_wmap5yr_noB/map_cmb_n0512_rxxxx_s31p22.fits';
end
if(~isfield(simopt,'separateABmaps'))
  simopt.separateABmaps=false;
end
if(~isfield(simopt,'coord'))
  simopt.coord='C';
  warning('simopt.coord unspecified. Defaulting to ''C'' for Celestial coordinates.')
end
if(~isfield(simopt,'siginterp'))
  simopt.siginterp='taylor';
end
if(~isfield(simopt,'force_ab_int_from_common_pix'))
  simopt.force_ab_int_from_common_pix=false;
end
if(~isfield(simopt,'interpix'))
  simopt.interpix=0.0625;
end
if(~isfield(simopt,'ukpervolt'))
  simopt.ukpervolt=get_ukpervolt();
end
if(~isfield(simopt,'type'))
  simopt.type='bicep';
end
if(~isfield(simopt,'ptsrc'))
  simopt.ptsrc='none';
end
if(~strcmp(simopt.ptsrc,'none')&~isfield(simopt,'src'))
  simopt.src.name={'PKS0537-441'};
  simopt.src.ra=84.7098392; 
  simopt.src.dec=-44.0858150;
  simopt.src.int=2*[0.55;0.3]; % in units of volts not Kelvin
end
if(~isfield(simopt,'beamcen'))
  simopt.beamcen='ideal';
end
if(~isfield(simopt,'beamwid'))
  simopt.beamwid='obs'; % causes beam widths to come from beamwid.csv file
end
if(~isfield(simopt,'beammapfilename'))
  simopt.beammapfilename=[];
end
if(~isfield(simopt,'diffpoint'))
  simopt.diffpoint='ideal';
end
if(~isfield(simopt,'chi'))
  simopt.chi='ideal';
end
if(~isfield(simopt,'epsilon'))
  simopt.epsilon='ideal';
end
if(~isfield(simopt,'polofs'))
  simopt.polofs='ideal';
end
if(~isfield(simopt,'lb'))
  simopt.lb=false;
end
if(~isfield(simopt,'xtalk'))
  simopt.xtalk=false;
end
if(~isfield(simopt,'tc'))
  simopt.tc=false;
end
if(~isfield(simopt,'xtalk_relgain'))
  simopt.xtalk_relgain=true;
end
if(~isfield(simopt,'rndcen'))
  simopt.rndcen=[0,0;0,0]; % 2 lines ra/dec
end
%  if(~isfield(simopt,'rndwid'))
%    simopt.rndwid=[0,0;0,0]; % 2 lines fwhm_major, fwhm_minor
%  end
%  if(~isfield(simopt,'rndalpha'))
%    simopt.rndalpha=[0,0];
%  end
if(~isfield(simopt,'rndsig'))
  simopt.rndsig=[0,0;0,0]; % 2 lines: pair sigma, pair to pair variation; A/B differential, A/B diff variation
end
if(~isfield(simopt,'rndelc'))
  simopt.rndelc=[0,0;0,0]; % 2 lines: pair c, pair to pair variation; A/B differential, A/B diff variation
end
if(~isfield(simopt,'rndelp'))
  simopt.rndelp=[0,0;0,0]; % 2 lines: pair p, pair to pair variation; A/B differential, A/B diff variation
end
if(~isfield(simopt,'rndeps'))
  simopt.rndeps=[0,0;0,0;0,0]; % 3 lines 100s, 150s, 220s - do these 3 work for previous years?
end
if(~isfield(simopt,'rndchi'))
  simopt.rndchi=[0,0;0,0;0,0]; % 3 lines 100s, 150s, 220s
end
if(~isfield(simopt,'abgain'))
  simopt.abgain=[1,0;1,0;1,0]; % 3 lines 100s, 150s, 220s 
end
if(~isfield(simopt,'pairmeanscp'))
  simopt.pairmeanscp=0;
end
if(~isfield(simopt,'sigmapbeamfwhms'))
  simopt.sigmapbeamfwhms=[];
end
if(~isfield(simopt,'maketod'))
  simopt.maketod=0;
end
if(~isfield(simopt,'update'))
  simopt.update=0;
end
if(~isfield(simopt,'mapopt') & ~simopt.maketod)
  simopt.mapopt=[];
end
if(~isfield(simopt,'striptod'))
  simopt.striptod=0;
end
if(~isfield(simopt,'rlz'))
  simopt.rlz=0;
end
if(~isfield(simopt,'state'))
  simopt.state=[];
end

% try to make random numbers truly random
if isempty(simopt.state)
  randn('state',sum(1e6*clock));
elseif isnumeric(simopt.state) & numel(simopt.state)<=2
  randn('state',simopt.state);
end
