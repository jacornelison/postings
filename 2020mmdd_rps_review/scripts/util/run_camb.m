function run_camb_planck2013(filename,scalars,tensors, r, n_t, p_t, do_tensor_neutrinos, n_s, p_s)
% function run_camb_planck2013(filename,scalars,tensors, r, n_t, p_t)
%  
% filename  : is not created from inputs! 
%            Choose something that makes senses, i.e. camb_planck2013_r0
% scalars   : 'T' (default), 'F'
% tensors   : 'T' (default), 'F'
% r         : '0' (default), '0.1' etc.
% n_t       : tensor spectral index
%           : '0' (default, also camb default)
% p_t       : tensor pivot scale
%           : '0.002' (default, also camb default)
% 
% NOTE: on odyssey, I've installed camb. You can use it by adding:
% source /n/home01/jetolan/exportrc
% to your .bashrc
%  
% cosmological parameters from planck 2013
% taken from:
% http://arxiv.org/pdf/1303.5076v1.pdf
% Table 5, Planck+WP
% also in camb_mar2013, base_planck_lowl_lowLike.ini
%  
% can be run by driver function:
% make_spectra.m
% eg:
% run_camb('camb_planck2013_r0','T', 'F', '0')
% run_camb('camb_planck2013_r1','T', 'T', '1')

if ~exist('scalars','var')
  scalars=[];
end

if ~exist('tensors','var')
  tensors=[];
end
  
if ~exist('r','var')
  r=[];
end

if ~exist('n_t','var')
  n_t=[];
end

if ~exist('p_t','var')
  p_t=[];
end

if ~exist('n_s','var')
  n_s=[];
end

if ~exist('p_s','var')
  p_s=[];
end


if ~exist('do_tensor_neutrinos','var')
  do_tensor_neutrinos=[];
end

if(isempty(scalars))
  scalars='T';
end

if(isempty(tensors))
  tensors='T';
end

if(isempty(r))
  r='0';
end

if(isempty(n_t))
  n_t='0';
end

if(isempty(p_t))
  p_t='0.002';
end

if(isempty(n_s))
  n_s='0.9619123';
end

if(isempty(p_s))
  p_s='0.05';
end

if(isempty(do_tensor_neutrinos))
  do_tensor_neutrinos='F';
end

% camb needs to be run from within its own dir
% find camb:
[err,camb_path] = system_safe('which camb_fits');
% if camb was not found assume we are already in the right folder
if err
  camb_path='.'
else
  camb_path=fileparts(camb_path);
end
currentdir=pwd;
cd(camb_path);

% put this into a try catch statement to be able to return to the original
% path in the end even if the run fails
try 
  tempdir='./tmp_camb';
  if exist(tempdir,'dir')
    system_safe(['rm -r ' tempdir]);
    system_safe(['mkdir ' tempdir]);
  else
    system_safe(['mkdir ' tempdir]);
  end

  h=fopen([tempdir '/camb.in'],'wt');
  fprintf(h,['output_root = ' filename '\n']);
  fprintf(h,['number_of_threads = 0\n']);
  fprintf(h,['get_scalar_cls = ' scalars '\n']);
  fprintf(h,['get_vector_cls = F\n']);
  fprintf(h,['get_tensor_cls = ' tensors '\n']);
  fprintf(h,['get_transfer = F\n']);
  fprintf(h,['do_lensing = T\n']);
  fprintf(h,['do_nonlinear = 2\n']);
  fprintf(h,['l_max_scalar = 2000\n']);
  fprintf(h,['k_eta_max_scalar = 14000\n']);
  fprintf(h,['l_max_tensor = 1500\n']);
  fprintf(h,['k_eta_max_tensor = 3000\n']);
  fprintf(h,['use_physical = T\n']);
  fprintf(h,['hubble = 67.04346\n']);
  fprintf(h,['temp_cmb = 2.7255\n']);
  fprintf(h,['ombh2 = 0.0220323\n']);
  fprintf(h,['omch2 = 0.1203761\n']);
  fprintf(h,['omnuh2 = 0.0006450616\n']);
  fprintf(h,['omk = 0\n']);
  fprintf(h,['nu_mass_eigenstates = 1\n']);
  fprintf(h,['nu_mass_degeneracies = 1.01533333333\n']);
  fprintf(h,['nu_mass_fractions = 1\n']);
  fprintf(h,['omega_baryon = 0.046\n']);
  fprintf(h,['omega_cdm = 0.224\n']);
  fprintf(h,['omega_lambda = 0.6914\n']);
  fprintf(h,['omega_neutrino = 0.0\n']);
  fprintf(h,['helium_fraction = 0.2476949\n']);
  fprintf(h,['massless_neutrinos = 2.03066666667\n']);
  fprintf(h,['massive_neutrinos = 1\n']);
  fprintf(h,['w = -1.0\n']);
  fprintf(h,['cs2_lam = 1\n']);
  fprintf(h,['reionization = T\n']);
  fprintf(h,['re_use_optical_depth = T\n']);
  fprintf(h,['re_optical_depth = 0.0924518\n']);
  fprintf(h,['re_redshift = 11.52\n']);
  fprintf(h,['re_delta_redshift = 0.5\n']);
  fprintf(h,['re_ionization_frac = -1\n']);
  fprintf(h,['initial_power_num = 1\n']);
  fprintf(h,['scalar_amp(1) = 2.21545e-9\n']);
  fprintf(h,['scalar_spectral_index(1) = ',n_s,'\n']);
  fprintf(h,['scalar_nrun(1) = 0.0\n']);
  fprintf(h,['tensor_spectral_index(1) = ' n_t '\n']);
  fprintf(h,['initial_ratio(1) = ' r '\n']);
  fprintf(h,['initial_condition = 1\n']);
  fprintf(h,['initial_vector = -1 0 0 0 0\n']);
  fprintf(h,['vector_mode = 0\n']);
  fprintf(h,['COBE_normalize = F\n']);
  fprintf(h,['CMB_outputscale = 7.42835025e12\n']);
  fprintf(h,['pivot_scalar = ',p_s,'\n']);
  fprintf(h,['pivot_tensor = ' p_t '\n']);
  fprintf(h,['transfer_high_precision = F\n']);
  fprintf(h,['transfer_interp_matterpower = T\n']);
  fprintf(h,['transfer_kmax = 2\n']);
  fprintf(h,['transfer_k_per_logint = 5\n']);
  fprintf(h,['transfer_num_redshifts = 1\n']);
  fprintf(h,['transfer_redshift(1) = 0\n']);
  fprintf(h,['feedback_level = 1\n']);
  fprintf(h,['lensing_method = 1\n']);
  fprintf(h,['accurate_BB = F\n']);
  fprintf(h,['accurate_reionization = F\n']);
  fprintf(h,['do_tensor_neutrinos = ',do_tensor_neutrinos,'\n']);
  fprintf(h,['massive_nu_approx = 1\n']);
  fprintf(h,['accurate_polarization = T\n']);
  fprintf(h,['do_late_rad_truncation = T\n']);
  fprintf(h,['RECFAST_fudge = 1.14\n']);
  fprintf(h,['RECFAST_fudge_He = 0.86\n']);
  fprintf(h,['RECFAST_Heswitch = 6\n']);
  fprintf(h,['accuracy_boost = 1\n']);
  fprintf(h,['l_accuracy_boost = 1\n']);
  fprintf(h,['high_accuracy_default = T\n']);
  fprintf(h,['l_sample_boost = 1\n']);
  fprintf(h,['scalar_output_file = scalcls.dat\n']);
  fprintf(h,['vector_output_file = vectcls.dat\n']);
  fprintf(h,['tensor_output_file = tenscls.dat\n']);
  fprintf(h,['total_output_file  = totcls.dat\n']);
  fprintf(h,['lensed_output_file = lensedcls.dat\n']);
  fprintf(h,['lensed_total_output_file = lensedtotcls.dat\n']);
  fprintf(h,['lens_potential_output_file = lenspotentialcls.dat\n']);
  fprintf(h,['FITS_filename      = scalcls.fits\n']);
  fprintf(h,['transfer_filename(1)    = transfer_out_z0.dat\n']);
  fprintf(h,['transfer_matterpower(1) = matterpower_z0.dat\n']);
  fclose(h)

  %%%%%%%%%%
  %run camb
  system_safe(['camb_fits ' tempdir '/camb.in']);

  %remove temp dir
  system_safe(['rm -r ' tempdir]);
  cd(currentdir);
catch err
  cd(currentdir);
  rethrow(err);
end

return
