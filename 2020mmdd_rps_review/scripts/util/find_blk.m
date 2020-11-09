function sc=find_blk(x,sampratio)
% sc=find_blk(x,sampratio)
%
% find start/end points of period x is >0

% Vectorized algorithm to find block margins
% (replaces while loop over TOD in which vectors grow)

% BEHAVIOR CHANGE:if a block ends at the end of x, the old algorithm took
% its end (incorrectly) to be length(x)-1.  The new algorithm correctly
% takes it to be length(x).

df = diff([0; (x(:)>0); 0]);
% df(i) = 1 if x(i) starts a block with x>0
sc.s = find(df==1);
% df(i+1) = -1 if x(i) ends a block with x>0
sc.e = find(df==-1) - 1;

% Introduce sf,ef - start,end fast sample
% (In BICEP 1st fast sample is contemporaneous with 1st slow sample
% according to antenna0.time.utcfast/slow)
if(exist('sampratio','var'))
  sc.sf=(sc.s-1)*sampratio+1;
  sc.ef=(sc.e-1)*sampratio+1;
end  

return
