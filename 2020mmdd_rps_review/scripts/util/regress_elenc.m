function d=regress_elenc(d,ind,fs)
% Regress the el encoder val from each bolo scan
% Doesn't seem to work - it tries to remove atmos var instead

  for i=1:length(fs.s)
  s=fs.sf(i); e=fs.ef(i);

  v=double(d.lockin.adcData(s:e,ind));


  % contruct the regressor
  clear X
  eloff=d.eloff(s:e);
  eloff=eloff-mean(eloff);
  X=[ones(size(eloff)),eloff];

  % loop over all channels
  for j=1:size(v,2)

    % regress
    b=regress(v(:,j),X);
    % subtract the correlated component
    c=v(:,j)-b(2)*eloff;

    if(1)
      [i j]
      b
      clf
      subplot(3,1,1)
      plot(v(:,j));
      subplot(3,1,2)
      plot(X*b,'r');
      subplot(3,1,3)
      plot(c,'g');
      pause
    end

    v(:,j)=c;

    bb(i,j)=b(2);
  end

  d.lockin.adcData(s:e,ind)=v;
  end

  keyboard

  return

