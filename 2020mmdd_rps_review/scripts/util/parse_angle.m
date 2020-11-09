% V = parse_angle(S)
%
% Parse an angle which may be in formats such as
% deg:min:sec
function v = parse_angle (angstr)
  s = 0;
  sgn = +1;
  pval = 1;
  angstr0 = angstr;
  while ~isempty(angstr)
    angstr = deblank(angstr);
    [tok angstr] = strtok(angstr, ':');
    tmp = str2num(tok);
    if isempty(tmp)
      error(['Could not parse position ' angstr0]);
    end
    if any(tok=='-')
      if pval==1
        sgn = -1;
      else
        error(['Could not parse position ' angstr0]);
      end
    end
    tmp = abs(tmp);
    s = s + tmp*pval;
    pval = pval / 60;
  end
  v = s * sgn;

  return

