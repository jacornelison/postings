function my_mask=get_mask_from_dead_cfg(file_name)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Read in a file from dead_lists and make a 33*16 mask from it.  You will
% need to point to the file yourself because these are not normally linked
% from the pipeline.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid=fopen(file_name);

% Skip through lines until the main part begins.
header_line=fgetl(fid);
while ~strcmp(header_line,'mask = [')
    header_line = fgetl(fid);
end

% Make 32 cell arrays that are 41 rows long.
cfg_data=cell(32,1);
for jj=1:32
    cfg_data(jj,:)=textscan(fid,'%u8',41,'delimiter',',','CommentStyle',{'/*','*/'});
end

fclose(fid);

my_mask=false(41,32);
for jj=1:32
    my_mask(:,jj)=logical(cfg_data{jj});
end

% We only want the first 33 rows and the first 16 columns.
my_mask(34:41,:)=[];
my_mask(:,17:32)=[];