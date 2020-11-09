function chopout = cleanchopref(ref)
%
% function chopout = cleanchopref(ref)
%
% TSG 2018-04-02
%
% Given a chop reference timestream (e.g. from FFBM),
% construct an "ideal" chop ref by interpolating to 
% higher precision and modulating the phase to account
% for any phase drift in the real chop ref.
%
% Input chop ref should be logical (1 for on, 0 for off).
%
% Inputs:
%      ref = logical chop ref timestream
%
% Output:
%  chopout.
%    ind_ud = FLOAT indices of expected 1 -> 0 transitions
%    ind_du = FLOAT indices of expected 0 -> 1 transitions
%    i_hires = vector containing complete hi-res FLOAT indices
%    refclean_hires = vector containing clean, hi-resolution chop ref
%    refclean = cleaned chop ref interpolated to same sampling as input
%    phasediff = phase difference (rad) btwn input and cleaned chop refs

% ideal chop ref is (this number) x real chop ref in resolution
res = 10;

if size(ref,1) == 1
  ref = ref';
end

if ~islogical(ref)
  error(['Chop reference b should be logical.']);
end

iref = (1:length(ref))';
iref_ideal = linspace(1,length(ref),length(ref)*res)';
dref = diff([ref(1); ref]);

%% Find peak frequency of chop ref (only use first chunk of data)
len = min(10000,length(ref));
maxFreq = 1 * get_maxFreq(ref,len)-9e-5;
Npercycle = 1 / maxFreq;

%% Construct ideal chop ref at high resolution.
% Start with arbitrary freq/phase, then correct the phase by
% calulcating the phase drift btwn real/ideal chop ref.
% Then after appropriate chop ref is made, recalculate phase drift
% to confirm it is flat at 0.
% First go around, use larger roll-off freq in filter to capture
% drift.  Second time, use much smaller one to flatten out phase and
% ignore outliers due to anomalous timing in chop ref.
frac_rolloff = [1.0/75 1.0/500 1.0/500];  
ph = zeros(length(iref_ideal),3);
for ii = 1:size(ph,2)
  ref_ideal = square(2*pi*iref_ideal*maxFreq - ph(:,ii));
  ref_ideal = ref_ideal > 0.5;
  dref_ideal = diff([ref_ideal(1); ref_ideal]);
  % Phase difference between real and ideal
  phase_diff = calcphasediff(iref,ref,dref,...
      iref_ideal,ref_ideal,dref_ideal,Npercycle,res);
  % In the possible scenario where the ideal ref has a transition in
  % the very beginning that doesn't exist in the real ref, this will
  % appear out of phase by 2pi.  Correct this.
  x = mean(phase_diff(1:floor(2*res/maxFreq)));
  if abs(x) > 1
    phase_diff = phase_diff - sign(x)*2*pi;
  end
  % Apply low-pass filter to phase.  Want to account for phase wander
  % and not point-to-point variations
  [bflt aflt] = butter(4,maxFreq*frac_rolloff(ii),'low');
  phase_diff_smth = filtfilt(bflt,aflt,phase_diff);
  % Apply phase correction
  ph(:,ii+1) = ph(:,ii) + (phase_diff_smth);
  % Plot for debugging
  if 0
    figure();
    plot(phase_diff,'.'); 
    %hold on; plot(phase_diff_smth,'r')
    pause
  end
end

% Check for any single hi or single low points in new ideal chop ref
if any(abs(diff(diff(ref_ideal)))==2)
  disp(['Warning: found ' num2str(sum( abs(diff(diff(ref_ideal)))>1 )) ...
	' instances of singleton chop samples in ideal chop ref (there ' ...
	' were ' num2str(sum( abs(diff(diff(ref)))>1 )) ' in the input ref).'])
end

%% Prep output
chopout.ind_du = iref_ideal(dref_ideal>0);
chopout.ind_ud = iref_ideal(dref_ideal<0);
chopout.i_hires = iref_ideal;
chopout.refclean_hires = ref_ideal;
chopout.refclean = interp1(iref_ideal,ref_ideal,iref);
chopout.refclean = chopout.refclean > 0.5;
chopout.phase_diff = phase_diff;

% Uncomment to make some plots.  Mostly for debugging for now.
%figure()
%plot(iref,ref,'b.')
%axis([0 50 -0.3 1.3])
%xlabel('Index')
%title('Original chop ref')
%figure()
%plot(iref_ideal,ref_ideal,'r-',iref,ref,'b.')
%axis([0 150 -0.3 1.3])
%xlabel('Index')
%title('Original and ideal chop ref')
%legend({'Hi-res ideal chop ref','Original chop ref'})
%figure(); hold on
%plot(phase_diff,'.')
%ylim([-2*pi 2*pi])
%xlabel('Index')
%ylabel('Phase difference')
%grid on

return

%%%%%%%%%%%%%%%%%%%%%%%%% SUB FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%
function maxFreq = get_maxFreq(ref,len)
  [chopSpec,chopFreq] = pwelch(ref(1:len),[],[],[],1);
  chopSpec(1:4) = 0;
  chopSpec((end-3):end) = 0; 
  [junk,imax] = max(chopSpec);
  maxFreq = chopFreq(imax);
return

function phase_diff = calcphasediff(iref,ref,dref,iref_ideal,ref_ideal,dref_ideal,Npercycle,res)
  % Calculate phase drift over input timestream.
  % Compute diffs for both real and ideal chop refs
  % (need to interp real to same res as ideal) and
  % calculate how far apart each 0 -> 1 transition is 
  % between real and ideal.
  ref_hires = interp1(iref,ref,iref_ideal);
  ref_hires = ref_hires >= 0.5;
  dref_hires = diff([ref_hires(1); ref_hires]); 
  % If phase drift > 2*pi, real phase drifted ahead ideal by one period;
  % if phase drift < 2*pi, real phase drifted behind ideal by one period.
  trans_real = iref_ideal(find(dref_hires==1));
  trans_ideal = iref_ideal(find(dref_ideal==1));
  minlength = min(length(trans_real),length(trans_ideal));
  pd = (trans_real(1:minlength)-trans_ideal(1:minlength))/Npercycle*2*pi;
  % If there are 2pi jumps, say due to skipped transitions in real ref,
  % the low pass filter won't perfectly capture them so manually correct
  % it here.
  %figure(); plot(pd,'.'); 
  for ii = 1:3
    jumphi = cumsum(diff(pd)>5);
    jumplo = cumsum(diff(pd)<-5);  
    pd(2:end) = pd(2:end) - 2*pi*(jumphi - jumplo);
  end
  %figure(); plot(pd,'.'); 
  % Sometimes there's a single stray point between 2*pi jumps, 
  % which otherwise escapes the algorithm.  Try to correct this
  % by diffing between points separated by two indices instead of one
  ddiff = pd(3:end) - pd(1:end-2);
  jumphi = cumsum(ddiff>5);
  jumplo = cumsum(ddiff<-5);  
  pd(3:end) = pd(3:end) - 2*pi*(jumphi - jumplo);
  % Sometimes the point-to-point scatter makes jumps double count...
  % cheap fix is to repeat the first step
  for ii = 1:3
    jumphi = cumsum(diff(pd)>5);
    jumplo = cumsum(diff(pd)<-5);  
    pd(2:end) = pd(2:end) - 2*pi*(jumphi - jumplo);
  end 
  %figure(); plot(pd,'.'); keyboard
  phase_diff = ...
      interp1(trans_real(1:minlength),pd,iref_ideal,'next','extrap');
return
