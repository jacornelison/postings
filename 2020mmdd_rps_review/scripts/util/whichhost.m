function host=whichhost()
% host=whichhost
% returns string 'spud', 'odyssey' or 'itasca'
% if it cannot determine which host, assume it
% is odyssey.

  persistent pHost;
  if isempty(pHost)
    h = hostname();
    % omega0 is just like spud
    if(~isempty(strfind(h,'spud')) | ~isempty(strfind(h,'omega0')))
      host='spud';
    elseif(~isempty(strfind(h,'node')))
      host='itasca';
    elseif(~isempty(strfind(h,'harvard')))
      host='odyssey';
    else
      host='odyssey';
    end
    pHost = host;
  end
  host = pHost;
return
