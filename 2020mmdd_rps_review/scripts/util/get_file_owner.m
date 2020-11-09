function owner = get_file_owner(filename)
  owner='';
  try
    [status,result]=system_safe(['ls -l ',filename]);
  catch ex
    status = 1;
  end  
  
  if ~status
    f = strfind(result,' ');
    owner = result(f(2)+1:f(3)-1);
  end 
return
