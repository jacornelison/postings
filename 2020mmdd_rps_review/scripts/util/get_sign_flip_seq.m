function sign_flip_seq=get_sign_flip_seq(coaddopt)
% sign_flip_seq=get_sign_flip_seq(coaddopt)
%
% Generate sign flip sequences which have equal total plus and
% minus weight. Separate (largely common) sequences are made for
% pair sum and pair diff. Using these balanced sequences should
% reduce signal leakage versus simple radnom sequences although it
% will not be perfect - we are assuming that all pairs vary
% together with time which is not completely true.
%
% For temporal jacks each half of the split needs to be balanced
% separately and the code below does this. One therefore needs to
% run reduc_coaddpairmaps separately for each such jack:
% 
% e.g.
% x=load('maps/1607/real_a_filtp3_weight3_gs_dp1100_jack2')
% coaddopt.sign_flip_seq=get_sign_flip_seq(x.coaddopt);
% coaddopt.jacktype='024689abc';
% reduc_coaddpairmaps(tags,coaddopt);
% for x=['1','3','5','7']
%   x=load(['maps/1607/real_a_filtp3_weight3_gs_dp1100_jack',x]);
%   coaddopt.sign_flip_seq=get_sign_flip_seq(x.coaddopt);
%   coaddopt.jacktype=x;
%   reduc_coaddpairmaps(tags,coaddopt);
% end

% find tag lists for two halves of jack split
tl1=find(sum(coaddopt.jackmask{1},2));
tl2=find(sum(coaddopt.jackmask{2},2));

if(isempty(intersect(tl1,tl2)))
  disp(sprintf('pure temporal jack detected (%s)',coaddopt.jackname));
  disp('making sign flip sequences for each half separately:');
  sign_flip_seq=get_sfs(coaddopt,tl1)+get_sfs(coaddopt,tl2);
else
  disp(sprintf('pairwise or mixed jack detected (%s)',coaddopt.jackname));
  disp('making single sign flip sequence:');
  sign_flip_seq=get_sfs(coaddopt,tl1);
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sign_flip_seq=get_sfs(coaddopt,tl)

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
sign_flip_seq=[ss;sd];

% check it worked
ss=ss(tl); sd=sd(tl);
ws=ws(tl); wd=wd(tl);
disp(sprintf('residual fractional pair sum signal %.4f',(sum(ws(ss))-sum(ws(~ss)))/sum(ws)));
disp(sprintf('residual fractional pair diff signal %.4f',(sum(wd(sd))-sum(wd(~sd)))/sum(wd)));

return
