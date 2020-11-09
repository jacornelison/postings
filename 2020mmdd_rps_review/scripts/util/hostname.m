function host=hostname()
% host=hostname()
%
% Returns the hostname of the machine.

  persistent host_p;

  if isempty(host_p)
    host = getenv('HOSTNAME');
    if isempty(host) || isempty(strfind(host, '.'))
      [r,host] = system_safe('hostname -f');
      host = strtrim(host);
    end
    host_p = host;
  end

  host = host_p;
end

