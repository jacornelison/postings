function [s,bpwf]=calc_expvals(s,inpmod,bpwf,supfac)
% [s,bpwf]=calc_expvals(s,inpmod,bpwf,supfac)
%
% use inpmod and bpwf to calculate expectation values of bandpowers

% multiply the bpwf by the filter/beam supression factor
if(exist('supfac'))
  bpwf=bpwf_supfac(bpwf,supfac);
end

%  content of the input model:
%  TT, TE, EE, BB, TB, EB, ET, BT, BE
%  content of the bpwf:
%  TT, TP, EE>EE, BB>BB, EE>BB, BB>EE, (PT)

% find the expectation values
for i=1:length(s)

  nbin=size(bpwf(i).Cs_l,2);
  
  % don't assume inpmod and bpwf start from same ell, shift as needed
  % inpmod starts with ell=0, but bpwf can start at higher ell
  nl=size(bpwf(i).Cs_l,1);
  indx=1:nl;
  indx=indx+find(bpwf(i).l(1)==inpmod.l)-1;

  % TT
  x=bpwf(i).Cs_l(:,:,1).*repmat(inpmod.Cs_l(indx,1),[1,nbin]);
  s(i).expv(:,1)=sum(x)';
  
  % TE
  x=bpwf(i).Cs_l(:,:,2).*repmat(inpmod.Cs_l(indx,2),[1,nbin]);
  s(i).expv(:,2)=sum(x)';

  % EE
  x1=bpwf(i).Cs_l(:,:,3).*repmat(inpmod.Cs_l(indx,3),[1,nbin]); % EE>EE * EE
  x2=bpwf(i).Cs_l(:,:,6).*repmat(inpmod.Cs_l(indx,4),[1,nbin]); % BB>EE * BB
  s(i).expv(:,3)=sum([x1;x2])';

  % BB
  x1=bpwf(i).Cs_l(:,:,4).*repmat(inpmod.Cs_l(indx,4),[1,nbin]); % BB>BB * BB
  x2=bpwf(i).Cs_l(:,:,5).*repmat(inpmod.Cs_l(indx,3),[1,nbin]); % EE>BB * EE
  s(i).expv(:,4)=sum([x1;x2])';
  
  % TB bpwf not explicitly calculated - use TE
  x=bpwf(i).Cs_l(:,:,2).*repmat(inpmod.Cs_l(indx,5),[1,nbin]);
  s(i).expv(:,5)=sum(x)';
  
  % EB bpwf not explicitly calculated - use EE
  x=bpwf(i).Cs_l(:,:,3).*repmat(inpmod.Cs_l(indx,6),[1,nbin]);
  s(i).expv(:,6)=sum(x)';
  
  % if freq/rx cross include exp vals of alternate crosses
  % field are getting rename in the reduc_final output, catch this here
  if isfield(s(1),'Cs_l')
    nspec = size(s(i).Cs_l,2);
  else
    nspec = size(s(i).real,2);
  end
  if(nspec==9)
    s(i).expv(:,7)=s(i).expv(:,2);
    s(i).expv(:,8)=s(i).expv(:,5);
    s(i).expv(:,9)=s(i).expv(:,6);    
  end
  
  % hmmm - there is an ET BPWF which could be used for ET
  % expval... why not?
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bpwf=bpwf_supfac(bpwf,supfac)

for i=1:length(bpwf)
  % interpolate the supfac values to every ell
  l=[0;supfac(i).l]; % prepend a zero to stop it going negative
  sftt=interp1(l,[0;supfac(i).tt],cvec(bpwf(i).l),'linear','extrap');
  sfee=interp1(l,[0;supfac(i).ee],cvec(bpwf(i).l),'linear','extrap');
  sfbb=interp1(l,[0;supfac(i).bb],cvec(bpwf(i).l),'linear','extrap');
  sfte=interp1(l,[0;supfac(i).te],cvec(bpwf(i).l),'linear','extrap');
  sftb=interp1(l,[0;supfac(i).tb],cvec(bpwf(i).l),'linear','extrap');
  sfeb=interp1(l,[0;supfac(i).eb],cvec(bpwf(i).l),'linear','extrap');
  
  % expand out to match the 3rd dim of bpwf
  %sf=[sftt,sfte,sfee,sfbb,sftb,sfeb];
  
  % the above is wrong elements 5 and 6 of bpwf are ee->bb and
  % bb->ee
  sf=[sftt,sfte,sfee,sfbb,sfee,sfbb];
  
  if(size(bpwf(i).Cs_l,3)>6)
    % reduc_bpwf may cause bpwf.Cs_l to go to 7 but the supfac.et is not
    % computed for jackf because it is an autospectrum
    if ~isfield(supfac,'et')
      % prevent crash
      warning('removing ET bpwf')
      bpwf(i).Cs_l=bpwf(i).Cs_l(:,:,1:6);
    else
      sfet=interp1(l,[0;supfac(i).et],cvec(bpwf(i).l),'linear','extrap');
      sf=[sf,sfet]; % ET
    end
  end
  
  % expand to match 2nd dim of bpwf
  n=size(sf);
  sf=reshape(sf,[n(1),1,n(2)]);
  
  % mult into bpwf
  n=size(bpwf(i).Cs_l);
  bpwf(i).Cs_l=bpwf(i).Cs_l.*repmat(sf,[1,n(2),1]);
end

% (re)normalize bpwf to unity sum using the same function used to do
% it at the end of reduc_bpwf.m
bpwf=norm_bpwf(bpwf);

return
