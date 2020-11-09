function mask = ffbm_cutting_pixels(experiment,year,rxNum)
% function mask = ffbm_cutting_pixels(experiment,rxNum,rxNum)
%
% Function to make a structure that determines which
% measurements will be used in the beam param files.
% The structure is saved, and then called in ffbm_param_manyfiles.m
% to summarize beam parameters measurements.
%
% Make sure there is a directory named mask_ffbm_year 
%
% Inputs:
%         experiment: 'bicep2' or 'keck'
%         year:       2012 or 2013 for keck
%         rxNum:      only relevant for keck
%
% CLW 16 June 2014
% CLW 18 June 2014, added code for 2014

% Find run times for all beam maps that will be used
[times dkangle sched tt1 tt2] = ffbm_findbmruntimes(experiment,'',year,rxNum);

% Make mask
mask = cutting_pixels_sub(experiment,year,rxNum,times);

% Save
save(['mask_ffbm_' num2str(year) '/mask_rx' ...
      int2str(rxNum)],'mask','rxNum')

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function det = cutting_pixels_sub(experiment,year,rxNum,times)

% Set thresholds based on experiment
switch experiment
  case 'bicep2'
    switch year
      case 2011
	[pp,ind] = get_array_info('20110201');
	pklow = 25;
	siglow = 0.10;
	sighigh = 0.30;
	elliphigh = 0.20;
      case 2012
	[pp,ind] = get_array_info('20120201');
	pklow = 25;
	siglow = 0.10;
	sighigh = 0.30;
	elliphigh = 0.20;
    end 
  
  case 'keck'
    switch year
      case 2012
	disp('2012');
      case 2013
	[pp ind] = get_array_info('20130201');
	cutind = (pp.rx == rxNum);
	pp = structcut(pp,cutind);
	ind = make_ind(pp);
	pklow = 25;
	siglow = 0.15;
	sighigh = 0.30;
	elliphigh = 0.20;
      case 2014
	[pp ind] = get_array_info('20140201');
	cutind = (pp.rx == rxNum);
	pp = structcut(pp,cutind);
	ind = make_ind(pp);
	switch pp.band(1)
	  case 100
	    pklow = 10;
	    siglow = 0.22;
	    sighigh = 0.37;
	    elliphigh = 0.20;
	  case 150
	    pklow = 25;
	    siglow = 0.15;
	    sighigh = 0.30;
	    elliphigh = 0.20;
	end
    end 
end % experiment switch

% Set up the mask
% 1 means cut completely
x = ones(1,528);
% Set rgls to 0 (usable)
x(ind.rgl) = 0;
ind.not = find(x);
clear x;

for ll = 1:length(times)
  times(ll)
  rxNum
  time = times(ll);

  % Load up the beammap for this time
  switch experiment
    case 'bicep2'
      % Get the dk angle from somewhere...
      load(['beammaps/map_' char(time) '.mat']);
      load(['beammaps/mapwin_' char(time) '.mat']);
    case 'keck'
      load(['beammaps/mapwin_' char(time) '_rx' int2str(rxNum) '.mat'])
  end
  
  % Get beam parameters, with the rotation option
  [sig e_p e_c] = ffbm_findingparam(A,pp,dk,1,experiment);
  ellip = sqrt(e_p.^2 + e_c.^2);
  siga = sig;  
  % For each of the parameters, find where the fit fails
  
  % Does the pixel exist?
  kk = 1; nanpix = [];
  for ii = 1:528;
    if isnan(A(:,ii))
      nanpix(kk) = ii;
      kk = kk + 1;
    end
  end
  
  % Peak height threshold
  kk = 1; pkcut = [];
  for ii = 1:528
    if A(1,ii) < pklow
      pkcut(kk) = ii;
      kk = kk + 1;
    end
  end
  
  % Sigma NaN
  kk = 1; nan_s = [];
  for ii = 1:528
    if isnan(sig(ii))
      nan_s(kk) = ii;
      kk = kk + 1;
    end
  end
  % Sigma thresholds
  jj = 1; badpix_sig = [];
  for ii = 1:528 
    if sig(ii) < siglow
      badpix_sig(jj) = ii;
      siga(ii) = NaN;
      jj = jj + 1;
    end
    if sig(ii) > sighigh
      badpix_sig(jj) = ii;
      siga(ii) = NaN;
      jj = jj + 1;
    end
  end
  
  % e_rgl = ellip;
  % e_rgl(ind.not) = NaN;
  
  % Ellipticity NaN
  jj = 1;
  kk = 1; nan_e = []; badpixe = [];
  for ii = 1:528
    if isnan(ellip(ii))
      nan_e(kk) = ii;
      kk = kk + 1;
    end
  end 
  % Ellipticity threshold
  for ii = 1:528
    if ellip(ii) > elliphigh
      badpixe(jj) = ii;
      % e_rgl(ii) = NaN;
      jj = jj + 1;
    end
  end
  
  % Dark SQUIDs (choose both darks and opens)
  darksquid = sort([ind.d ind.o]);
  
  if(0) % Old hard-coded version
  darksquid=[1 34 67 100 133 166 199 232 496 463 430 397 364 331 298 ...
	265] +1;
  if year == 2014
    if rxNum == 0 | rxNum == 2;
      darksquid=[...
	1 2 3     4     5     6     7     8     9    10    11    12    13 ...
	14    15    16    17    18    19    20    33    34    35    36    37    50 ...
	51    52    53    66    67    68    69    70    83    84    85    86    99 ...
	100   101   102   103   116   117   118   119   120   121   122   123   124 ...
	125   126   127   128   129   130   131   132   133   134   135   136   137 ...
	138   139   140   141   142   143   144   145   146   147   148   149   150 ...
	151   152   165   166   167   168   169   182   183   184   185   198   199 ...
	200   201   202   215   216   217   218   231   232   233   234   235   248 ...
	249   250   251   252   253   254   255   256   257   258   259   260   261 ...
	262   263   264   265   266   267   268   269   270   271   272   273   274 ...
	275   276   277   278   279   280   281   282   283   284   297   298   299 ...
	300   301   314   315   316   317   330   331   332   333   334   347   348 ...
	349   350   363   364   365   366   367   380   381   382   383   384   385 ...
	386   387   388   389   390   391   392   393   394   395   396   397   398 ...
	399   400   401   402   403   404   405   406   407   408   409   410   411 ...
	412   413   414   415   416   429   430   431   432   433   446   447   448 ...
	449   462   463   464   465   466   479   480   481   482   495   496   497 ...
	498   499   512   513   514   515   516   517   518   519   520   521   522 ...
	523   524   525   526   527   528];
    end
  end
  end % if(0)
    
  % Cut based on mirror position
  switch experiment
    case 'bicep2'
      rxpos = NaN;
      poscut = zeros(1,528);
    case 'keck'
      % Figure out where each rx is -- Pos3 is the base of drum
      rxpos = dk_rx_pos_ffbm(rxNum,dk);

      % Now cut on the far field flat
      % Two mirror positions, so use the date of the beam map to determine
      % which one to use
      % Time formatting sucks...
      mirror_change_date = '2014-02-22 02:00:00';
      bm_date = sprintf('%s-%s-%s %s:%s:%s',...
                        time{1}(1:4),time{1}(5:6),time{1}(7:8),...
		        time{1}(10:11),time{1}(12:13),time{1}(14:15));
      if datenum(bm_date) < datenum(mirror_change_date)		
	[mask100 mask150] = ffbm_maskMaker('back');
      else
	[mask100 mask150] = ffbm_maskMaker('front');
      end
      
      % Choose 100/150 GHz mask
      freq_change_date = '2014-01-01 00:00:00';
      if datenum(bm_date) < datenum(freq_change_date)
	poscut = mask150(rxpos,:);
      else
	switch pp.band(1)
	  case 100
	    poscut = mask100(rxpos,:);
	  case 150
	    poscut = mask150(rxpos,:);
	end
      end % 100/150
  end % Mirror position cut
    
  % Make a bunch of masks
  det(ll).poscut = zeros(1,528);
  det(ll).darksquid = zeros(1,528);
  det(ll).nofits = zeros(1,528);
  det(ll).badsig = zeros(1,528);
  det(ll).badellip = zeros(1,528);
  det(ll).nansigellip = zeros(1,528);
  det(ll).pkcut = zeros(1,528);
  det(ll).notrgl = zeros(1,528);
  det(ll).handcut = zeros(1,528);
  
  % Set bad pixels = 1
  % pk<25, 0.1<sig<0.3, ellip<0.15
  det(ll).rxNum = rxNum;
  det(ll).time = time;
  det(ll).dk_rx_pos = rxpos;
  det(ll).poscut = poscut;
  det(ll).darksquid(darksquid) = 1;
  det(ll).nofits(nanpix) = 1;
  det(ll).badsig(badpix_sig) = 1;
  det(ll).badellip(badpixe) = 1;
  det(ll).nansigellip(nan_s) = 1;
  det(ll).nansigellip(nan_e) = 1;
  det(ll).pkcut(pkcut) = 1;
  det(ll).notrgl(ind.not) = 1;
  
  clear pkcut nan_e nan_s badpixe badpix_sig nanpix poscut rxpos
end % Loop over beam maps

% Super fun hand cuts
% Choose the indices to kill, and set them to 1 in each beam map (ll)
switch experiment
  
  case 'bicep2'
    disp('Yay, no BICEP2 hand cuts!')
    
  case 'keck'

    switch year
      case 2012
	
	switch rxNum
	  case 0
	    handcut = [321,349];
	    det(1).handcut(handcut) = 1;
	    handcut = [84,116];
	    for ii = 1:2
	      det(ii).handcut(handcut) = 1;
	    end
	    handcut = [193,306];
	    det(3).handcut(handcut) = 1;
	    handcut = [341,342];
	    det(4).handcut(handcut) = 1;
	  case 1
	    handcut = [27,31,83,130];
	    for ii = 1:12
	      det(ii).handcut(handcut) = 1;
	    end
	    handcut = [309];
	    for ii = [3,6,10,11]
	      det(ii).handcut(handcut) = 1;
	    end
	    handcut = [321];
	    det(12).handcut(handcut) = 1;
	    handcut = [242];
	    det(5).handcut(handcut) = 1;
	    det(3).handcut(59) = 1;
	  case 2
	    handcut = [43,113,136,42,132,112,165,135,283,133];
	    det(8).handcut(handcut) = 1;
	    handcut = [51];
	    det(7).handcut(handcut) = 1;
	    handcut = [132];
	    for ii = 1:13
	      det(ii).handcut(handcut) = 1;
	    end
	    handcut = [235,311,310,461];
	    det(2).handcut(handcut) = 1;
	    det(13).handcut(232) = 1;
	    handcut = [307 408];
	    det(9).handcut(232) = 1;
	    det(3).handcut(461) = 1;
	  case 3
	    handcut = [268,350,387,383,267,283,349,366,382,284,505];
	    det(2).handcut(handcut) = 1;
	    det(3).handcut([527 516 526]) = 1;
	    det(6).handcut(18) = 1;
	    handcut = [8,57,107,222,278,282,309,325,342,391,408,...
		       7,56,106,215,221,308,324,341,357,374,390,...
		       407,440,523,506,522];
	    det(7).handcut(handcut) = 1;
	    handcut = [47,57,5];
	    det(4).handcut(handcut) = 1;
	    det(5).handcut([5,251]) = 1;
	    det(9).handcut([215,216]) = 1;
	    det(11).handcut(251) = 1;
	    det(12).handcut([99,215]) = 1;
	    det(13).handcut([99,65]) = 1;
	  case 4
	    handcut = [330];
	    for ii = 1:length(det)
	      det(ii).handcut(handcut) = 1;
	    end
	    handcut = [294];
	    det(1).handcut(handcut) = 1;
	    handcut = [320,321];
	    det(2).handcut(handcut) = 1;
	    handcut = [284,174];
	    det(3).handcut(handcut) = 1;
	    handcut = [262,307];
	    det(4).handcut(handcut) = 1;
	    handcut = [248,307];
	    det(5).handcut(handcut) = 1;
	    handcut = [175,263,294,293,174];
	    det(6).handcut(handcut) = 1;
	    handcut = [100,263];
	    det(7).handcut(handcut) = 1;
	    handcut = [303,313,302,320];
	    det(10).handcut(handcut) = 1;
	    handcut = [263];
	    det(14).handcut(handcut) = 1;
	    
	end % rx switch
	
      case 2013
	disp('No Keck 2013 hand cuts!');
	
      case 2014
	% Look at the actual beam maps
	% If there is an obviously wrong map, run
	% gcp = det2gcp(ro,co,pol,ti)
	% Then set (gcp + 1) = 1  -- Matlab indexing
	% Each rx has a different set of times, so can't go on beam map
	% number -- just choose based on actual beam map time
	% To save time, use the maskMaker to figure out which pixels have a
	% chance of making it by the position cuts and only cut those which
	% are zero
	switch rxNum
	  case 0
	    for ll = 1:length(times)
	      switch times{ll}
		case '20140218_154747'
		  handcut = [164,187,211,186,163];
		  det(ll).handcut(handcut) = 1;
		case '20140220_041958'
		  handcut = [226,207,225,238];
		  det(ll).handcut(handcut) = 1;
		case '20140222_024932'
		  handcut = [222,226,225];
		  det(ll).handcut(handcut) = 1;
		case '20140222_235122'
		  handcut = [260,444,427,425,419,494,445,428];
		  det(ll).handcut(handcut) = 1;
		case '20140223_073616'
		  handcut = [379];
		  det(ll).handcut(handcut) = 1;
		case '20140224_115627'
		  handcut = [22,470,21];
		  det(ll).handcut(handcut) = 1;
		case '20140226_183344'
		  handcut = [47,72,32,82,98,360,46,54,71,81,97,...
			     368,359];
		  det(ll).handcut(handcut) = 1;
		case '20140304_172032'
		  handcut = [494,461,451,435,420,450,442,419,358,340,...
			     375,356,357,339,374,355];
		  det(ll).handcut(handcut) = 1;
		case '20140305_010515'
		  handcut = [239,228,245,187,204,181,238,244,186,203,...
			     180,196,418,469];
		  det(ll).handcut(handcut) = 1;
		case '20140307_021320'
		  handcut = [453,224,226,222,204,214,336,360,354,223,...
			     225,221,238,203,213,345,335,359,353,337,...
			     357,341];
		  det(ll).handcut(handcut) = 1;
	      end
	    end

	  case 1
	    for ll = 1:length(times)
	      switch times{ll}
		case '20140214_100044'
		  handcut = [311,294,321,358,309,323,309,308,339];
		  det(ll).handcut(handcut) = 1;
		case '20140218_154747'
		  handcut = [30,103,100,67,350,319,311,294,272,346,...
			     319,311,272,46,54,81,102,502,396,302,...
			     318,339,380];
		  det(ll).handcut(handcut) = 1;
		case '20140218_233005'
		  handcut = [280,286,311,294,338,321,305,358,342,373,...
			     356,174,285,269,326,310,293,337,320,304,...
			     357,341,372,355];
		  det(ll).handcut(handcut) = 1;
		case '20140220_120023'
		  handcut = [308,289];
		  det(ll).handcut(handcut) = 1;
		case '20140221_092016'
		  handcut = [160,156,410,422,159,155,421];
		  det(ll).handcut(handcut) = 1;
		case '20140223_151409'
		  handcut = [263,218,381,418,443,453,490,474,455,439,215,...
			     380,425,469,489,473,454,438];
		  det(ll).handcut(handcut) = 1;
		case '20140224_194022'
		  handcut = [518,519,507];
		  det(ll).handcut(handcut) = 1;
		case '20140225_031815'
		  handcut = [525,524];
		  det(ll).handcut(handcut) = 1;
		case '20140225_110234'
		  handcut = [226,222,212,228,204,214,201,166,331,348,350,...
			     362,369,371,223,225,221,211,227,203,219,213,...
			     198,363,347,361,368,351,370,390,374];
		  det(ll).handcut(handcut) = 1;
		case '20140226_031317'
		  handcut = [499,511,501,498,510,500];
		  det(ll).handcut(handcut) = 1;
		case '20140227_193606'
		  handcut = [298,334,379,362,472,333,378,471,438];
		  det(ll).handcut(handcut) = 1;
		case '20140228_183710'
		  handcut = [131,119,100,164,387,391,397,416,402,410,130,...
			     118,85,151,182,386,390];
		  det(ll).handcut(handcut) = 1;
		case '20140301_225622'
		  handcut = [213,214];
		  det(ll).handcut(handcut) = 1;
		case '20140305_010515'
		  handcut = [22,39,49,65,53,70,67,117,38,48,64,52,69,99,...
		             383,379,352,344,327,378,351,343];
		  det(ll).handcut(handcut) = 1;
		case '20140305_141953'
		  handcut = [98,67,97,85,99,397,416,412,402,410,429,415,...
		             409];
		  det(ll).handcut(handcut) = 1;
		case '20140306_053943'
		  handcut = [459,439,469,454,438];
		  det(ll).handcut(handcut) = 1;
		case '20140307_095358'
		  handcut = [129,115,103,119,117,128,120,102,447,433,451,...
			     470,474,446,469,473];
		  det(ll).handcut(handcut) = 1;
	      end
	    end
	    
	  case 2
	    for ll = 1:length(times)
	      switch times{ll}
		case '20140215_181103'
		  handcut = [117,206,195,154,148,526,116,225,242,...
			     205,178,193,153,147,526];
		  det(ll).handcut(handcut) = 1;
		case '20140218_233005'
		  handcut = [486,457,455,439,508,485,467,473,456,...
			     454,438];
		  det(ll).handcut(handcut) = 1;
		case '20140219_071426'
		  handcut = [327,338,285,326,336,373,337];
		  det(ll).handcut(handcut) = 1;
		case '20140221_170444'
		  handcut = [362,352,360,354,375,225,210,194,...
			     508,361,345,351,376,359,370,353,374];
		  det(ll).handcut(handcut) = 1;
		case '20140222_024932'
		  handcut = [327,325,307,326,324,306];
		  det(ll).handcut(handcut) = 1;
		case '20140222_102404'
		  handcut = [181,379,362,346,329,369,352,428,180,...
			     378,361,345,328,351,335,434,458];
		  det(ll).handcut(handcut) = 1;
		case '20140223_151409'
		  handcut = [445,428,418,435,461,493,460,444,427,353,...
			     376,341];
		  det(ll).handcut(handcut) = 1;
		case '20140226_031317'
		  handcut = [245,220,214,418,443,453,219,244,213,...
			     417,442,452];
		  det(ll).handcut(handcut) = 1;
		case '20140227_021434'
		  handcut = [239,245,237,238,244,236,];
		  det(ll).handcut(handcut) = 1;
		case '20140227_193606'
		  handcut = [356,355,374,375];
		  det(ll).handcut(handcut) = 1;
		case '20140228_031217'
		  handcut = [230,220,329,313,296,311,338,321,325,356,323,...
			     219,246,473,456,328,312,302,310,320,324,322];
		  det(ll).handcut(handcut) = 1;
		case '20140307_021320'
		  handcut = [55,65,54,64,196,352,342,341,461,451,443,437,...
			     441,460,450,442,436,440];
		  det(ll).handcut(handcut) = 1;
	      end
	    end
	  
	  case 3
	    for ll = 1:length(times)
	      switch times{ll}
		case '20140219_071426'
		  handcut = [515,459,437,514,510,436,419,407];
		  det(ll).handcut(handcut) = 1;
		case '20140219_203955'
		  handcut = [43,16,42,11,110,73,62,15,281];
		  det(ll).handcut(handcut) = 1;
		case '20140221_092016'
		  handcut = [107,131,385,393,373,110,106,130,384,...
			     392,374,372];
		  det(ll).handcut(handcut) = 1;
		case '20140221_170444'
		  handcut = [140,156,206,221,173,239,123,107,90,74,...
			     57,41,24,8,133,515,499,449,433,416,...
			     400,527,511,461,445,428,412,383,367,...
			     350,334,317,301,268,395,379,362,346,...
			     329,313,280,385,368,352,319,303,286,...
			     270,7,23,40,56,73,89,106,122,139,155,...
			     172,205,221,238,254,514,498,448,432,...
			     415,399,526,510,460,444,427,411,382,...
			     366,349,333,316,300,267,394,378,361,...
			     328,295,279,384,368,351,335,318,302,...
			     285,269];
         	  det(ll).handcut(handcut) = 1;
		case '20140222_024932'
		  handcut = [133,265,379,419,403,440,423,454,438,349,...
		             361,345,378,361,384,368,392,420];
		  det(ll).handcut(handcut) = 1;
		case '20140222_102404'
		  handcut = [6,22,65,70,86,67,117,364,367,362,21,38,64,...
			     48,64,69,85,116,38,21,5,366,361,335,318,310];
		  det(ll).handcut(handcut) = 1;
		case '20140224_115627'
		  handcut = [284,371];
		  det(ll).handcut(handcut) = 1;
		case '20140225_110234'
		  handcut = [65,53,34,133,430,449,433,445,435,426,437,...
			     420,422,64,52,448,432,444,434,442,425,436,...
			     419,423,421];
		  det(ll).handcut(handcut) = 1;
		case '20140226_031317'
		  handcut = [6,32,49,37,53,34,84,395,379,369,352,393,321,...
			     305,5,31,48,36,52,413,382,394,378,368,351,...
			     343,320,304,291];
		  det(ll).handcut(handcut) = 1;
		case '20140226_105325'
		  handcut = [133];
		  det(ll).handcut(handcut) = 1;
		case '20140227_193606'
		  handcut = [281];
		  det(ll).handcut(handcut) = 1;
		case '20140228_031217'
		  handcut = [459,426,393,360,327,524,508,442,425,409,392,...
		             376,522,506,456,440,423,407,390,374,357,341,...
			     324,308,291,275];
		  det(ll).handcut(handcut) = 1;
		case '20140228_105707'
		  handcut = [418];
		  det(ll).handcut(handcut) = 1;
		case '20140304_094006'
		  handcut = [261,220,214,169,133,244,260,219,213,168];
		  det(ll).handcut(handcut) = 1;
		case '20140307_095358'
		  handcut = [164,152,147,163,151,346,345,335];
		  det(ll).handcut(handcut) = 1;
	      end
	    end
	    
	  case 4
	    for ll = 1:length(times)
	      switch times{ll}
		case '20140214_022039'
		  handcut = [243,259,258];
		  det(ll).handcut(handcut) = 1;
		case '20140215_181103'
		  handcut = [26,45,61,74,105,121,131,150,133,298,...
			     301,296,270,356,25,44,60,56,73,95,104,...
			     130,165,300,295,269,355];
		  det(ll).handcut(handcut) = 1;
		case '20140216_015133'
		  handcut = [246];
		  det(ll).handcut(handcut) = 1;
		case '20140218_080933'
		  handcut = [442,426,508,476,443,410,488,472,454,438,...
			     422,269,393,377,360,344,311,294,356,340,...
			     323,248,458,442,424,409,521,505,486,453,...
			     437,421,376,359,343,326,310,293,355,339,321];
		  det(ll).handcut(handcut) = 1;
		case '20140219_203955'
		  handcut = [160,243,255,220,159,188,254,211,219,246];
		  det(ll).handcut(handcut) = 1;
		case '20140220_041958'
		  handcut = [76,61,14,208,226,261,75,60,44,23,40,13,...
			     258,188];
		  det(ll).handcut(handcut) = 1;
		case '20140220_194118'
		  handcut = [208,243,108];
		  det(ll).handcut(handcut) = 1;
		case '20140221_170444'
		  handcut = [43,59,45,61,41,63,80,74,55,72,65,82,53,...
			     70,86,84,67,193,189,195,204,214,202,447,...
			     430,397,433,416,445,428,435,418,402,426,...
			     410,420,404,424,408,284,296,286,42,58,44,...
			     60,40,56,73,79,54,71,64,81,97,69,85,83,...
			     99,192,188,203,213,231,461,446,432,415,...
			     444,427,434,417,401,425,409,419,403,423,...
			     407,357,322,355];
		  det(ll).handcut(handcut) = 1;
		case '20140222_024932'
		  handcut = [4,18,463,480,430,449,433,445,428,418,402,...
			     387,375,358,356,340,15,17,479,462,448,432,...
			     444,427,417,401,390,374,370,357,355,339];
		  det(ll).handcut(handcut) = 1;
		case '20140222_102404'
		  handcut = [259,255,261,284,269,258];
		  det(ll).handcut(handcut) = 1;
		case '20140222_235122'
		  handcut = [228,220,214,202,216,362,352,371,375,227,215,...
			     347,361,351,370,374,455,469];
		  det(ll).handcut(handcut) = 1;
		case '20140224_115627'
		  handcut = [393,377,523,375,358,342,325,309,392,376,271,...
			     522,506,489,473,440,423,407,390,341,324,308,...
			     523,507,490,474,457,441,424,408];
		  det(ll).handcut(handcut) = 1;
		case '20140224_194022'
		  handcut = [441,422,469,438,269];
		  det(ll).handcut(handcut) = 1;
		case '20140226_031317'
		  handcut = [259,255,284,296,286,258,254,283,295,285];
		  det(ll).handcut(handcut) = 1;
		case '20140226_105325'
		  handcut = [22,5,21,15];
		  det(ll).handcut(handcut) = 1;
		case '20140226_183344'
		  handcut = [261,237,214,230,185,202,133,400,428,435,476,...
			     486,472,455,439,422,260,236,213,229,184,201,...
			     165,182,399,427,434,469,452,438,421];
		  det(ll).handcut(handcut) = 1;
		case '20140227_021434'
		  handcut = [12,8,14,11,7,13];
		  det(ll).handcut(handcut) = 1;
		case '20140228_031217'
		  handcut = [16,4,18,51,37,513,499,482,466,494,484,451,...
			     492,474,441,472,5,15,3,17,21,36,50,512,479,...
			     498,481,465,493,483,450,491,475,458,469,473,...
			     440];
		  det(ll).handcut(handcut) = 1;
		case '20140228_105707'
		  handcut = [524,502];
		  det(ll).handcut(handcut) = 1;
		case '20140228_183710'
		  % Kill tiles 1/4
		  handcut = [1:132,397:528,175,208,160,177,156,162,195,...
			     138,154,187,148,136,169,150,385,270,174,207,...
			     240,159,176,155,161,153,186,147,135,168,182];
		  det(ll).handcut(handcut) = 1;	  
		case '20140304_172032'
		  handcut = [232,264,420,404,441,424,408,472,455,439,422,...
			     458,469,452,419,403,440,423,407,454,438,421,...
			     298,315,265,350,317,301,379,346,352,393,307,...
			     315,330,314,297,349,333,316,300,378,345,368,...
			     392,322,306];
		  det(ll).handcut(handcut) = 1;
		case '20140305_141953'
		  handcut = [424,439,473,440,438,208,210,189,195,187,181,...
			     169,207,209,258,188,194,186,180,168,182,367,...
			     385,393,366,384,392];
		  det(ll).handcut(handcut) = 1;
		case '20140306_132038'
		  handcut = [169,168,182,348,334,387,347,333,368];
		  det(ll).handcut(handcut) = 1;
	      end
	    end
	
	end  % rx switch
    end % year switch
end    % experiment switch
    

return 


