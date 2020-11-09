function blfile=get_bl_filename(band)
% blfile=get_bl_filename(band)
%
% Returns the path to data file(s) containing a given beam transfer
% function, B_l.
%
% INPUTS
%   band    One of the following recognized names:
%
%             WMAP 9yr bands:
%             'wmap9k','wmap9ka','wmap9q','wmap9v','wmap9w'
%
%             WMAP 7yr bands:
%             'wmap7k','wmap7ka','wmap7q','wmap7v','wmap7w'
%
%             Planck nominal release (R1.10) bands:
%             'pl100','pl143','pl217'
%
%             Planck beam files shared under the BKP MOU (DX11d):
%             'P030mou', 'P044mou', 'P070mou', 'P100mou', 'P143mou'
%             'P217mou', 'P353mou'
%
%             BK_100: 'Kuber100', 'Kuber100rev1'
%             BK_150: 'B2uberA','B2uberB','B2bbns'
%             BK_220: 'Kuber220'
%
% OUTPUTS
%   blfile    Path(s) to the requested data files. In some cases (e.g. WMAP),
%             may be a cell array of provided files (e.g. one for each
%             "detector").

  switch band
    case 'wmap9k'
      blfile = 'input_maps/wmap9/bl/wmap_ampl_bl_K1_9yr_v5p1.txt';
    case 'wmap9ka'
      blfile = 'input_maps/wmap9/bl/wmap_ampl_bl_Ka1_9yr_v5p1.txt';
    case 'wmap9q'
      blfile = {...
        'input_maps/wmap9/bl/wmap_ampl_bl_Q1_9yr_v5p1.txt', ...
        'input_maps/wmap9/bl/wmap_ampl_bl_Q2_9yr_v5p1.txt' ...
      };
    case 'wmap9v'
      blfile = {...
        'input_maps/wmap9/bl/wmap_ampl_bl_V1_9yr_v5p1.txt', ...
        'input_maps/wmap9/bl/wmap_ampl_bl_V2_9yr_v5p1.txt' ...
      };
    case 'wmap9w'
      blfile = {
        'input_maps/wmap9/bl/wmap_ampl_bl_W1_9yr_v5p1.txt', ...
        'input_maps/wmap9/bl/wmap_ampl_bl_W2_9yr_v5p1.txt', ...
        'input_maps/wmap9/bl/wmap_ampl_bl_W3_9yr_v5p1.txt', ...
        'input_maps/wmap9/bl/wmap_ampl_bl_W4_9yr_v5p1.txt' ...
      };

    case 'wmap7k'
      blfile = 'input_maps/wmap7/bl/wmap_K1_ampl_bl_7yr_v4.txt';
    case 'wmap7ka'
      blfile = 'input_maps/wmap7/bl/wmap_Ka1_ampl_bl_7yr_v4.txt';
    case 'wmap7q'
      blfile = {...
        'input_maps/wmap7/bl/wmap_Q1_ampl_bl_7yr_v4.txt', ...
        'input_maps/wmap7/bl/wmap_Q2_ampl_bl_7yr_v4.txt' ...
      };
    case 'wmap7v'
      blfile = {...
        'input_maps/wmap7/bl/wmap_V1_ampl_bl_7yr_v4.txt', ...
        'input_maps/wmap7/bl/wmap_V2_ampl_bl_7yr_v4.txt' ...
      };
    case 'wmap7w'
      blfile = {...
        'input_maps/wmap7/bl/wmap_V1_ampl_bl_7yr_v4.txt', ...
        'input_maps/wmap7/bl/wmap_V2_ampl_bl_7yr_v4.txt', ...
        'input_maps/wmap7/bl/wmap_V3_ampl_bl_7yr_v4.txt', ...
        'input_maps/wmap7/bl/wmap_V4_ampl_bl_7yr_v4.txt' ...
      };

    case 'pl100'
      blfile = 'input_maps/planck/planck_beam/b_l_100GHz.fits';
    case 'pl143'
      blfile = 'input_maps/planck/planck_beam/b_l_143GHz.fits';
    case 'pl217'
      blfile = 'input_maps/planck/planck_beam/b_l_217GHz.fits';

    case 'P030mou'
      blfile = 'input_maps/b2planck_140811/beams/Planck_B_ell_030GHz_DX11d_FULL.fits';
    case 'P044mou'
      blfile = 'input_maps/b2planck_140811/beams/Planck_B_ell_044GHz_DX11d_FULL.fits';
    case 'P070mou'
      blfile = 'input_maps/b2planck_140811/beams/Planck_B_ell_070GHz_DX11d_FULL.fits';
    case 'P100mou'
      blfile = 'input_maps/b2planck_140806/beams/Planck_B_ell_100GHz_DX11d_FULL.fits';
    case 'P143mou'
      blfile = 'input_maps/b2planck_140806/beams/Planck_B_ell_143GHz_DX11d_FULL.fits';
    case 'P217mou'
      blfile = 'input_maps/b2planck_140806/beams/Planck_B_ell_217GHz_DX11d_FULL.fits';
    case 'P353mou'
      blfile = 'input_maps/b2planck_140806/beams/Planck_B_ell_353GHz_DX11d_FULL.fits';

    case 'Kuber100'
      warning('Kuber100 is deprecated. Kuber100rev1 should be used instead.');
      blfile = 'aux_data/beams/beamfile_20141028_sum_100.fits';
    case 'Kuber100rev1'
      blfile = 'aux_data/beams/beamfile_20150321_sum_100.fits';
    case 'B2uberA'
      blfile = 'aux_data/beams/beamfile_uberchopper_A.fits';
    case 'B2uberB'
      blfile = 'aux_data/beams/beamfile_uberchopper_B.fits';
    case 'B2bbns'
      blfile = 'aux_data/beams/beamfile_20130222_sum.fits';
    case 'K24in210'
      blfile = 'aux_data/beams/beamfile_20160101_sum_210.fits';
    case 'Kuber220'
      blfile = 'aux_data/beams/beamfile_20150101_sum_220.fits';

    case 'SPTpol150'
      blfile = 'aux_data/beams/sptpol_150Ghz_beam.fits';

    case 'B24in100'
      blfile = 'aux_data/beams/beamfile_20160101_sum_100.fits';
 
    case 'B3gaussdeconv'
      blfile = 'aux_data/beams/beamfile_20171221_gausscorrfix_sum_100.fits';
      
    case 'B3polycorr'
      blfile = 'aux_data/beams/beamfile_20180206_polycorr_sum_100.fits';
    case 'K210polycorr'
      blfile = 'aux_data/beams/beamfile_20180206_polycorr_sum_210.fits';

    case 'K270prelim'
      blfile = 'aux_data/beams/beamfile_20180724_prelim_sum_270.fits';
      
    otherwise
      error('Unknown beam transfer function (B_l) name: %s', band);
  end
end

