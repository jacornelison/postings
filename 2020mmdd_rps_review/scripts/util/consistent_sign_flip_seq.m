function [seq1,seq2]=consistent_sign_flip_seq(coaddopt1,coaddopt2,taglist2)
% [seq1,seq2]=consistent_sign_flip_seq(coaddopt1,coaddopt2,taglist2)
%
% Generate sign flip sequences that obey the following two constraints:
%   1) Equal total plus and minus weight. Separate (largely common)
%      sequences are made for pair sum and pair diff. Using these
%      balanced sequences should reduce signal leakage versus simple
%      random sequences although it will not be perfect - we are
%      assuming that all pairs vary together with time which is not
%      completely true.
%   2) Simultaneous tags in the two experiments get the same sign.
%
% For temporal jacks each half of the split needs to be balanced
% separately and the code below does this. One therefore needs to
% run reduc_coaddpairmaps separately for each such jack.
% 
% Inputs:
%   coaddopt1: from a map coadd in this experiment
%   coaddopt2: from a map coadd in the other experiment
%   taglist2:  tag_list.csv structure, path, or filename for the
%              other experiment
%
% Outputs:
%   seq1: sign flip sequence for this experiment
%   seq2: sign flip sequence for the other experiment

% if we've been given map file names instead of coaddopt structures, load now
if ischar(coaddopt1)
  coaddopt1=load(coaddopt1,'coaddopt');
  coaddopt1=coaddopt1.coaddopt;
end
if ischar(coaddopt2)
  coaddopt2=load(coaddopt2,'coaddopt');
  coaddopt2=coaddopt2.coaddopt;
end

% read in tag list CSV files
% For this experiment: use standard working directory layout
if ~isfield(coaddopt1,'ptag')
  coaddopt1.ptag=ParameterRead('aux_data/tag_list.csv');
end
% For the other experiment: need some specification
if ~isfield(coaddopt2,'ptag') && exist('taglist2','var') && ~isempty(taglist2)
  if isstruct(taglist2)
    coaddopt2.ptag=taglist2;
  elseif ischar(taglist2)
    if exist(taglist2,'dir')
      taglist2=fullfile(taglist2,'tag_list.csv');
    end
    coaddopt2.ptag=ParameterRead(taglist2);
  else
    error(['Bad input type for taglist2.  Should be a tag_list.csv structure, file, or path.']);
  end
end

% find tag lists for two halves of jack split
tl11=find(sum(coaddopt1.jackmask{1},2));
tl12=find(sum(coaddopt1.jackmask{2},2));
tl21=find(sum(coaddopt2.jackmask{1},2));
tl22=find(sum(coaddopt2.jackmask{2},2));

% check jacktype consistency
if coaddopt1.jacktype~=coaddopt2.jacktype
  error(['Mismatched jacktypes: 1=' coaddopt1.jackname ', 2=' coaddopt2.jackname]);
end

% call joint sign flip sequence function, either once for all, or in
% two jackknife splits.  Note that simultaneous tags in opposite
% halves of a jackknife will not get matching signs.
if(isempty(intersect(tl11,tl12)) && isempty(intersect(tl21,tl22)))
  disp(sprintf('pure temporal jack detected (%s)',coaddopt1.jackname));
  disp('making sign flip sequences for each half separately:');
  [seq1,seq2]=get_sfs_joint(coaddopt1,tl11,coaddopt2,tl21);
  [seq3,seq4]=get_sfs_joint(coaddopt1,tl21,coaddopt2,tl22);
  seq1=seq1+seq3;
  seq2=seq2+seq4;
else
  disp(sprintf('pairwise or mixed jack detected detected (%s)',coaddopt1.jackname));
  disp('making single sign flip sequence:');
  [seq1,seq2]=get_sfs_joint(coaddopt1,tl11,coaddopt2,tl21);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main function that does the work.
%   1. Ensure tag list 1 is the shorter of the two (swap if needed).
%   2. Match up tags between the two lists.
%   3. Use simple equal-weight function for shorter tag list.
%   4. Use constrained equal-weight function for longer tag list.
%   5. Verify that signs are consistent between matched tags.
function [seq1,seq2]=get_sfs_joint(coaddopt1,tl1,coaddopt2,tl2)

% 1. Let coaddopt1 be the one with the shorter tag list.  We'll choose signs for this one
% first.  If it's the other way around, swap them now.
if length(coaddopt1.tags)>length(coaddopt2.tags)
  [seq2,seq1]=get_sfs_joint(coaddopt2,tl2,coaddopt1,tl1);
  return
% If they have equal numbers of tags... make something up.  Use weights instead.
elseif length(coaddopt1.tags)==length(coaddopt2.tags)
  hstmp1=structcat(1,coaddopt1.hsmax);
  hstmp2=structcat(1,coaddopt2.hsmax);
  if nansum(hstmp1.d(:))>nansum(hstmp2.d(:))
    [seq2,seq1]=get_sfs_joint(coaddopt2,tl2,coaddopt1,tl1);
  end
end

% 2. Find simultaneous tags in the two lists
j2=tag_match_up(coaddopt2.tags,coaddopt2.ptag,coaddopt1.tags,coaddopt1.ptag);
j2=j2(tl2);

% 3. Choose sign flips first for only the expt with fewer tags.  This way it will be
% weighted as evenly as possible.
disp(['Expt. 1 (' num2str(length(tl1)) ' tags, ' num2str(sum(j2~=0)) ' in common):']);
seq1=get_sfs_simple(coaddopt1,tl1);

% 4. Call constrained sequence magic
disp(['Expt. 2 (' num2str(length(tl2)) ' tags, ' num2str(sum(j2~=0)) ' in common):']);
seq2=get_sfs_constrained(coaddopt2,tl2,j2,seq1);

% 5. Check it worked
if ~check_seq_consistency(coaddopt1,tl1,seq1,coaddopt2,tl2,seq2)
  error(['Failed to construct consistent sign flip sequence for two experiments.  Number of non-overlap tags may be too small.']);
else
  disp(['Verified sign flips are consistent for ' num2str(sum(j2~=0)) ' common tags.']);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Return true if all matching tags have consistent signs, false otherwise
function c=check_seq_consistency(coaddopt1,tl1,seq1,coaddopt2,tl2,seq2)

j1=tag_match_up(coaddopt1.tags(tl1),coaddopt1.ptag,coaddopt2.tags(tl2),coaddopt2.ptag);
seq1=seq1(:,tl1);
seq2=seq2(:,tl2);
ns=sum(seq1(1,j1>0)~=seq2(1,j1(j1>0)));
nd=sum(seq1(2,j1>0)~=seq2(2,j1(j1>0)));

c=(ns==0) & (nd==0);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a modified version of get_sfs_simple.  It uses
% the same basic approach, but instead of randomly permuting all
% the tags, it pre-sorts the ones that have already been assigned
% and randomly permutes the others.  The permutation (and therefore
% the resulting sequence) are as consistent as possible between
% pair sum and pair diff.
function seq=get_sfs_constrained(coaddopt,tl,j,seq0)

% concatenate the weights
hsmax=structcat(1,coaddopt.hsmax);

% take the sum of the weights in each tag
ind=coaddopt.ind;
ws=nansum(hsmax.w(:,ind.la)');
wd=nansum(hsmax.w(:,ind.lb)');

% setup sign sequences for sum and diff
ss=false(size(coaddopt.tags)); sd=ss;

% choose a random permutation but don't apply yet
perm=randperm(length(tl));

% handle pair sum and pair diff separately, but with same
% random permutation for unconstrained tags

% 1. pair sum

% put in order:
%   - tags constrained to be seq=1
%   - tags unconstrained, with random permutation
%   - tags constrained to be seq=0
seqc=zeros(size(tl));
seqc(j>0)=seq0(1,j(j>0));
tlc1=tl(j>0 & seqc==1);
tlc2=tl(j==0);
[tmp i]=sort(perm(j==0));
tlc2=tlc2(i);
tlc3=tl(j>0 & seqc==0);
tlc=[tlc1;tlc2;tlc3];

% take cumulative sum of weights in this order
wsc=cumsum(ws(tlc));

% find the 50% transition point
ts=find(diff((wsc<wsc(end)/2)));

% we make all the tags before this point in randomized list + and
% the ones after -
ss(tlc(1:ts))=true;

% 2. pair diff

% put in order
seqc=zeros(size(tl));
seqc(j>0)=seq0(2,j(j>0));
tlc1=tl(j>0 & seqc==1);
tlc2=tl(j==0);
[tmp i]=sort(perm(j==0));
tlc2=tlc2(i);
tlc3=tl(j>0 & seqc==0);
tlc=[tlc1;tlc2;tlc3];

% take cumulative sum of weights in this order
wdc=cumsum(wd(tlc));

% find the 50% transition point
td=find(diff((wdc<wdc(end)/2)));

% we make all the tags before this point in randomized list + and
% the ones after -
sd(tlc(1:td))=true;

% make sum/diff sequences rows of output array
seq=[ss;sd];

% check it worked
ss=ss(tl); sd=sd(tl);
ws=ws(tl); wd=wd(tl);
disp(sprintf('  residual fractional pair sum signal %.4f',(sum(ws(ss))-sum(ws(~ss)))/sum(ws)));
disp(sprintf('  residual fractional pair diff signal %.4f',(sum(wd(sd))-sum(wd(~sd)))/sum(wd)));

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is directly from get_sign_flip_seq.m.
function seq=get_sfs_simple(coaddopt,tl)

% concatenate the weights
hsmax=structcat(1,coaddopt.hsmax);

% take the sum of the weights in each tag
ind=coaddopt.ind;
ws=nansum(hsmax.w(:,ind.la)');
wd=nansum(hsmax.w(:,ind.lb)');

% setup sign sequences for sum and diff
ss=false(size(coaddopt.tags)); sd=ss;

% randomly permute tag list
tl(randperm(length(tl)))=tl;

% take cummulative sum of weights in this order
wsc=cumsum(ws(tl));
wdc=cumsum(wd(tl));

% find the 50% transition points
ts=find(diff((wsc<wsc(end)/2)));
td=find(diff((wdc<wdc(end)/2)));

% we make all the tags before this point in randomized list + and
% the ones after -
ss(tl(1:ts))=true;
sd(tl(1:td))=true;

% make sum/diff sequences rows of output array  
seq=[ss;sd];

% check it worked
ss=ss(tl); sd=sd(tl);
ws=ws(tl); wd=wd(tl);
disp(sprintf('  residual fractional pair sum signal %.4f',(sum(ws(ss))-sum(ws(~ss)))/sum(ws)));
disp(sprintf('  residual fractional pair diff signal %.4f',(sum(wd(sd))-sum(wd(~sd)))/sum(wd)));

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find matching tags between the two experiments' tag lists.
% j1 and j2 will be indices into the other tag list; if zero,
% there was no matching tag.  Don't assume tag names will
% be the same.  Instead use start times with 10 minute tolerance.
function [j1,j2]=tag_match_up(tags1,p1,tags2,p2)
tic

t1=datenum(p1.tstart,'dd-mmm-yyyy:HH:MM:SS');
t2=datenum(p2.tstart,'dd-mmm-yyyy:HH:MM:SS');

[c1 i1]=ismember(tags1,p1.tag);
[c2 i2]=ismember(tags2,p2.tag);

t1c=t1(i1);
t2c=t2(i2);

DELTA_MAX=10/60/24; % match up to +/-10 minutes

j1=interp1(t2c,1:length(t2c),t1c,'nearest','extrap');
j1(abs(t1c-t2c(j1))>DELTA_MAX)=0;

if nargout>1
  j2=interp1(t1c,1:length(t1c),t2c,'nearest','extrap');
  j2(abs(t2c-t1c(j2))>DELTA_MAX)=0;
end

s=toc;
% disp(['Spent ' num2str(s) ' seconds matching up tag lists.']);

return


