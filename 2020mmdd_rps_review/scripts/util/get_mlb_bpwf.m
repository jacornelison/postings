function bpwf=get_mlb_bpwf(f)
% bpwf=get_mlb_bpwf(f)
%
% Get MLB bandpower window function from file and convert to same
% style as mine

lab={'100GHz','150GHz','cross'};

if(strfind(f,'070305'))
  for i=1:length(lab)
    x=load(sprintf('/data/mlb/s11_simple_spectra/windows/%s_tt.dat',lab{i}));
    x=reshape(x,[2698,30,3]); bpwf(i).wl(:,:,1)=x(:,:,3);
    x=load(sprintf('/data/mlb/s11_simple_spectra/windows/%s_te.dat',lab{i}));
    x=reshape(x,[2698,30,3]); bpwf(i).wl(:,:,2)=x(:,:,3);
    x=load(sprintf('/data/mlb/s11_simple_spectra/windows/%s_ee.dat',lab{i}));
    x=reshape(x,[2698,30,4]); bpwf(i).wl(:,:,3)=x(:,:,3); bpwf(i).wl(:,:,5)=x(:,:,4);
    x=load(sprintf('/data/mlb/s11_simple_spectra/windows/%s_bb.dat',lab{i}));
    x=reshape(x,[2698,30,4]); bpwf(i).wl(:,:,4)=x(:,:,3); bpwf(i).wl(:,:,6)=x(:,:,4);
    bpwf(i).l=x(:,1,2);
  end
else
  for i=1:length(lab)
    x=load(sprintf('mlb_spec/bpwf_%s_tt.dat',lab{i}));
    x=reshape(x,[2700,30,3]); bpwf(i).wl(:,:,1)=x(:,:,3);
    x=load(sprintf('mlb_spec/bpwf_%s_te.dat',lab{i}));
    x=reshape(x,[2700,30,3]); bpwf(i).wl(:,:,2)=x(:,:,3);
    x=load(sprintf('mlb_spec/bpwf_%s_ee.dat',lab{i}));
    x=reshape(x,[2700,30,4]); bpwf(i).wl(:,:,3)=x(:,:,3); bpwf(i).wl(:,:,5)=x(:,:,4);
    x=load(sprintf('mlb_spec/bpwf_%s_bb.dat',lab{i}));
    x=reshape(x,[2700,30,4]); bpwf(i).wl(:,:,4)=x(:,:,3); bpwf(i).wl(:,:,6)=x(:,:,4);
    bpwf(i).l=x(:,1,2);
  end
end

for i=1:length(bpwf)
  % convert from stupid W_l to quantity we actually mult theory
  % spectrum by
  l=repmat(bpwf(i).l,[1,30,6]);
  bpwf(i).Cs_l=bpwf(i).wl.*2*pi./(l.^2.*(l+1));
end  

return
