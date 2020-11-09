function user=whoami()
% user=whoami()
%
% Returns the user name of the current user

  persistent user_p;

  if isempty(user_p)
    user = getenv('USER');
    if isempty(user)
      [r,user] = system_safe('whoami');
      user = strtrim(user);
    end
    user_p = user;
  end

  user = user_p;
end
