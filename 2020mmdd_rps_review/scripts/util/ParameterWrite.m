function ParameterWrite(filename,p,k)
% ParameterWrite(filename,p,k)
%
% Matlab version of Ganga's IDL ParameterWrite.pro
%                   Kiwon's IDL aux_data_write.pro
%
% e.g.: ParameterWrite('fp_pars_20050108_with_notes.csv',p,k)
%
%  ['integer', 'long', 'float', 'double', 'string', 'long64']
%
% original: 3/21/07 JMK
% modified:


file_id = fopen(filename,'wt');

% begin header portion of file -------------------

% write out comment lines
for index = 1:length(k.comments)
  if ~any(k.comments{index}(1:2)=='#')  % prepend '# ' if not there already
    k.comments{index} = ['# ' k.comments{index}];
  end
  fprintf(file_id,[k.comments{index} '\n']);
end

% write out keywords and values
keywords = fields(rmfield(k,{'comments','fields','units','formats'}));
max_keyword_width = size(char(keywords),2);
for index = 1:length(keywords)
  %eval(['line = [blanks(max_keyword_width+2) k.' keywords{index} '];'])
  eval(['val = k.' keywords{index} ';']);
  if ~isstr(val)
    val = num2str(val);
  end
  line = [blanks(max_keyword_width+2) val];
  line(1:length(keywords{index})+1) = [upper(keywords{index}) ','];  % make uppercase
  fprintf(file_id,'%s \n',line);
end

% begin tabular portion of file -------------------
fprintf(file_id,'BODY \n');

% assemble the column for each field as a cell array,
%   then convert to a uniform-width char array and append it
fieldnames = k.fields;
ncols = length(fieldnames);
eval(['nrows = length(p.' fieldnames{1} ');'])
table = [];
for col_index = 1:ncols
  if col_index==ncols
    col_sep = '';
  else
    col_sep = ', ';
  end
  colhead{1}=[upper(k.fields{col_index}) col_sep];
  colhead{2}=[k.units{col_index} col_sep];
  colhead{3}=[k.formats{col_index} col_sep];
  if strcmp(k.formats{col_index},'string')
    % putting col_sep in braces makes char() preserve whitespace
    colval = char(strcat(p.(fieldnames{col_index}), {col_sep}));
  else
    % roundabout syntax gets left-alignment, preserves whitespace
    colval = strtrim(cellstr(num2str(reshape(p.(fieldnames{col_index}),[],1),'%.15g')));
    colval = char(strcat(colval, {col_sep}));
  end
  table = [table char(colhead{:},colval)];  % append column to table 
end

% write out the table
for index = 1:size(table,1)
  fprintf(file_id,'%s \n',table(index,:));
end

fclose(file_id);

