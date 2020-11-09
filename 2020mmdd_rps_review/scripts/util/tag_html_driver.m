function tag_html_driver(tag_csv_file,stats_csv_file,tags,html_dir,reducdir_url,experiment)

%written as a driver for tag_html so you can farm smaller portions of tag browser html
%  Process tags in pt.tags min_ii:max_ii.  must input min_ii and max_ii
%
% tag_csv_file:    file with tag_list
% stats_csv_file:  file with tag_stats
% html_dir:        directory of the browser (reducplots/browser)
% reducdir_url:    (reducplots)
% experiment:      keck or BICEP2

if ~exist('tag_csv_file','var') || isempty(tag_csv_file)
  tag_csv_file=fullfile('aux_data','tag_list.csv');
end
if ~exist('stat_csv_file','var') || isempty(stats_csv_file)
  stats_csv_file=fullfile('aux_data','tag_stats.csv');
end
if ~exist('html_dir','var') || isempty(html_dir)
  html_dir=fullfile('reducplots','browser');
end
if ~exist(html_dir,'dir')
  mkdir(html_dir);
end
if ~exist('reducdir_url','var') || isempty(reducdir_url)
  reducdir_url=fullfile('..','.');
end
if ~exist('experiment','var')
  experiment=get_experiment_name;
end 

%load in the parameters
[pt,kt]=ParameterRead(tag_csv_file);
[ps,ks]=ParameterRead(stats_csv_file);

for ii=1:length(tags)

  disp(tags{ii})
  %write the html
  write_tag_html(pt,ps,ks,tags{ii},html_dir,reducdir_url,experiment);

end
