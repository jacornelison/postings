function [tagsublist,mtl]=get_tag_sublist(coaddopt,passcuts,firstorlast)
% [tagsublist,mtl]=get_tag_sublist(coaddopt)
%
% Make a tag "flavor" substitution list
%
% e.g.
% x=load('maps/1446/real_y_filtp3_weight3_gs_dp1100_jack0')
% coaddopt.tagsublist=get_tag_sublist(x.coaddopt)
% reduc_coaddpairmaps(x.coaddopt.tags,coaddopt)
%
% The passcuts option is to make all mtl pass cuts.  

if ~exist('passcuts','var') | isempty(passcuts)
  passcuts=true;
end
if ~exist('firstorlast','var') | isempty(firstorlast)
  firstorlast='last';
end


% make cuts
if passcuts
  pass=any(coaddopt.c.c2.overall,2);
end

% make a two column array of el and dk offsets
% (eloff has some values which differ by epsilon so round them)
pl=[round(coaddopt.traj.eloff*1000)/1000 coaddopt.traj.dk];
if isfield(coaddopt.traj, 'racen')
  % Add racen if it exists, rounded to the nearest multiple of 6 degrees
  % (approximately half the distance the sky drifts in a 50 minute scan set).
  pl=[pl round(coaddopt.traj.racen/6)*6];
end
if passcuts
  pl_cut = pl(pass,:);
end

% find the unique rows and index arrays -
% For B2 there are 164 entries corresponding to 8 deck angles x 20
% elevation steps (4 "early" deck angles and 4 "normal" ones) plus 4
% scansets at eloff=5 deg which occur during the G phase for the
% "normal" deck angle set.
[ur,i,j]=unique(pl,'rows',firstorlast);
if passcuts; [ur_cut,ic,jc]=unique(pl_cut,'rows',firstorlast); end

% find the corresponding minimal list of tags - this is somewhat
% "messy" - presumably it could consist only of B and C phases but
% instead it contains a funny mix because a given el/dk offset
% first occurs during some other phase
if passcuts
  tags_tmp=coaddopt.tags(pass);
  mtl=tags_tmp(ic);

  %get the ordering right.  in principle, they should be the exact same because they are ordered.  
  %However, you can imagine a case where the cuts completely remove a configuration 
  [lia,locb]=ismember(ur,ur_cut,'rows');
  if any(lia==0)
    warning('you have some strange dk/el configs that do not pass cuts.  ignoring them');
    locb(locb==0)=1;
  end
  mtl=mtl(locb);

else
  mtl=coaddopt.tags(i);
end

% make the subsitution list - each tag in first row should be
% subsituted with the corresponding one in the second row
tagsublist=[coaddopt.tags;mtl(j)];

return
