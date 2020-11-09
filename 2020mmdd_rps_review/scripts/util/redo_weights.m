function subdatarw=redo_weights(subdata,realhs,ind,style)
%function subdatarw=redo_weights(subdata,realhs,ind,style)
%
% Takes a flavor-subset substitution ac or matrix structure and multiplies by
% the mean weight and variance across half-scans from the real pairmap.
%
% INPUTS
%   subdata    Data structure from substitution tag. The structure needs to
%              conform to the specified style parameter.
%
%   realhs     Half-scan structure from real data pairmap.
%
%   ind        Detector index structure.
%
%   style      Identifies which type of data is to be reweighted. Valid options
%              are:
%
%                  'ac'     -> Traditional pipeline ac structure (Default)
%                  'matrix' -> Matrix pipeline matrix structure
%
% OUTPUTS
%   subdatarw  Reweighted data structure.
%

  if ~exist('style','var') || isempty(style)
    style = 'ac';
  end

  % Get the mean weights and variances from the substituted map and real map
  % for pair sum and pair diff. Pair sum channel pairs are stored in ind.a and
  % pair diff channels in ind.b.

  % For weights constant over the scanset (weight3) this is exact.
  srw=nanmean(realhs.w(:,ind.a),1);
  drw=nanmean(realhs.w(:,ind.b),1);

  % The half scan variances fluctuate over the scanset (quite a lot) so there
  % is no way to recover the result which would have been obtained pre-binning
  % into pixels.
  srv=nanmean(realhs.s(:,ind.a).^2,1);
  drv=nanmean(realhs.s(:,ind.b).^2,1);

  % Since the weights and variances are taken from the real pairmaps, their
  % means across halfscans (calculated above) can still have NaNs, so change
  % those to zeros so that sparse arrays will still contain mostly zeros. This
  % saves a lot of memory and calculation time
  srw(isnan(srw))=0; drw(isnan(drw))=0;
  srv(isnan(srv))=0; drv(isnan(drv))=0;

  switch style
    case 'ac'
      subdatarw = reweight_ac(subdata, srw, drw, srv, drv);
    case 'matrix'
      subdatarw = reweight_matrix(subdata, srw, drw, srv, drv);
    otherwise
      error('redo_weights:UnknownStyle', sprintf(...
        'The style `%s` is an unknown reweighting data product.', style));
  end
end


function subacrw=reweight_ac(subac, srw, drw, srv, drv)
  % Pre-allocate structure. This saves only a small amount of time, but is
  % still a good idea.
  subacrw=subac;

  % The fields "sitime" and "ditime" do not include weights or variances, so
  % they do not appear in the loop and will remain unchanged from their values
  % in the pre-allocated structure.

  for sdir=1:size(subac,2)  % over scan direction
    for j=1:size(subac,1) % over detector pairs

      % THE FIELDS BELOW ARE CREATED IN REDUC_MAKEPAIRMAPS - AND WHEN
      % THOSE ARE CHANGED THEY NEED TO BE CHANGED HERE ALSO
      % HOWEVER - THIS WHOLE MESS CAN BE DRAMATICALLY
      % IMPROVED BY LOOPING OVER LISTS OF FIELD NAMES

      %  RE-WEIGHT PAIR SUM QUANTITIES

      % quantities needed to make T map and var-map
      subacrw(j,sdir).wsum=subac(j,sdir).wsum.*srw(j);
      subacrw(j,sdir).wz=subac(j,sdir).wz.*srw(j);
      subacrw(j,sdir).wwv=subac(j,sdir).wwv.*srw(j).^2.*srv(j);    
      % extra quantities needed to make raw and weighted
      % integration time maps
      subacrw(j,sdir).switime=subac(j,sdir).switime.*srw(j);
      subacrw(j,sdir).swmax=subac(j,sdir).swmax*srw(j);
      % poly sub and ground sub components    
      if(isfield(subac,'wz_psub'))
        subacrw(j,sdir).wz_psub=subac(j,sdir).wz_psub.*srw(j);
        subacrw(j,sdir).wz_gsub=subac(j,sdir).wz_gsub.*srw(j);
      end

      % RE-WEIGHT PAIR DIFF QUANTITIES

      % quantities needed to make Q/U maps and var-maps
      subacrw(j,sdir).w=subac(j,sdir).w.*drw(j);
      subacrw(j,sdir).wcz=subac(j,sdir).wcz*drw(j);
      subacrw(j,sdir).wsz=subac(j,sdir).wsz*drw(j);
      subacrw(j,sdir).wcc=subac(j,sdir).wcc*drw(j);
      subacrw(j,sdir).wss=subac(j,sdir).wss*drw(j);
      subacrw(j,sdir).wcs=subac(j,sdir).wcs*drw(j);
      subacrw(j,sdir).wwccv=subac(j,sdir).wwccv.*drw(j).^2.*drv(j);
      subacrw(j,sdir).wwssv=subac(j,sdir).wwssv.*drw(j).^2.*drv(j);
      subacrw(j,sdir).wwcsv=subac(j,sdir).wwcsv.*drw(j).^2.*drv(j);
      % this one is needed for pair diff signal subtraction in noise estimation
      subacrw(j,sdir).wzdiff=subac(j,sdir).wzdiff.*drw(j);
      % this one is needed to make (approx) coaddtype=3 pair diff maps
      if(isfield(subac,'wc'))
        subacrw(j,sdir).wc=subac(j,sdir).wc.*drw(j);
      end
      % extra quantities needed to make raw and weighted
      % integration time maps
      subacrw(j,sdir).dwitime=subac(j,sdir).dwitime.*drw(j);
      subacrw(j,sdir).dwmax=subac(j,sdir).dwmax.*drw(j); 
      % poly sub and ground sub components
      if(isfield(subac,'wcz_psub'))
        subacrw(j,sdir).wcz_psub=subac(j,sdir).wcz_psub.*drw(j);
        subacrw(j,sdir).wsz_psub=subac(j,sdir).wsz_psub.*drw(j);
        subacrw(j,sdir).wcz_gsub=subac(j,sdir).wcz_gsub.*drw(j);
        subacrw(j,sdir).wsz_gsub=subac(j,sdir).wsz_gsub.*drw(j);
      end
      % deprojection templates
      if(isfield(subac,'wcd')) 
        for k=1:size(subac(j,sdir).wcd,2)
          subacrw(j,sdir).wcd{k}=subac(j,sdir).wcd{k}.*drw(j);
          subacrw(j,sdir).wsd{k}=subac(j,sdir).wsd{k}.*drw(j);
        end
      end

    end % over detector pairs
  end % over scan direction
end

function submatrw=reweight_matrix(submat, srw, drw, srv, drv)
  % Pre-allocate structure. This saves only a small amount of time, but is
  % still a good idea.
  submatrw = submat;

  for sdir=1:length(submat)  % over scan direction
    for j=1:length(submat.pair) % over detector pairs
      % THE FIELDS BELOW ARE CREATED IN REDUC_MATRIX - AND WHEN THOSE ARE
      % CHANGED THEY NEED TO BE CHANGED HERE ALSO.

      %  RE-WEIGHT PAIR SUM QUANTITIES

      submatrw(sdir).pair(j).awgfha = submat(sdir).pair(j).awgfha * srw(j);
      submatrw(sdir).pair(j).awa    = submat(sdir).pair(j).awa    * srw(j);
      submatrw(sdir).pair(j).awta   = submat(sdir).pair(j).awta   * srw(j);
      submatrw(sdir).pair(j).awwva  = submat(sdir).pair(j).awwva  * srw(j)^2 * srv(j);


      %  RE-WEIGHT PAIR DIFF QUANTITIES

      submatrw(sdir).pair(j).awssa = submat(sdir).pair(j).awssa * drw(j);
      submatrw(sdir).pair(j).awcca = submat(sdir).pair(j).awcca * drw(j);
      submatrw(sdir).pair(j).awcsa = submat(sdir).pair(j).awcsa * drw(j);
      submatrw(sdir).pair(j).awpa  = submat(sdir).pair(j).awpa  * drw(j);

      submatrw(sdir).pair(j).awdiffa = submat(sdir).pair(j).awdiffa * drw(j)^2 * drv(j);
      submatrw(sdir).pair(j).awwssva = submat(sdir).pair(j).awwssva * drw(j)^2 * drv(j);
      submatrw(sdir).pair(j).awwccva = submat(sdir).pair(j).awwccva * drw(j)^2 * drv(j);
      submatrw(sdir).pair(j).awwcsva = submat(sdir).pair(j).awwcsva * drw(j)^2 * drv(j);

      submatrw(sdir).pair(j).awcgfcha_awcgfsha = submat(sdir).pair(j).awcgfcha_awcgfsha * drw(j);      
      submatrw(sdir).pair(j).awsgfcha_awsgfsha = submat(sdir).pair(j).awsgfcha_awsgfsha * drw(j);      
    end
  end
end