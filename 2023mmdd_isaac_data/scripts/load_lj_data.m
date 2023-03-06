function data = load_lj_data(fname)

% Labjack takes the form of 
%%
f = fopen(fname,'r');
header = fgetl(f);

format = repmat('%f',1,length(strsplit(header,',')));
%format = [format '%s'];
headnames = strsplit(header,',');

data_raw = textscan(f,format,'delimiter',',');%,'HeaderLines',1)
%keyboard()
data = [];
for hi = 2:length(headnames)
data.(headnames{hi}) = data_raw{hi};
end




fclose(f);
