function fitswrite(data,filename,m,n)
% fitswrite(data,filename)
%
% In the beginning there was fitswrite function written by the guy
% below. JET modified this to make fitswrite_table and
% fitswrite_bintable which are in the healpix subdir. They are
% driven by the badly named write_fits_map and write_fits_alms
% which write fits table files suitable for input into anafast and
% synfast respectively.
%
% Here I will lightly modify the original function to allow to write
% our ra/dec rectangular grid maps. Recent versions of Matlab have
% their own fitswrite function but this does not seem to easily
% allow the addition of user keywords making it pretty
% useless. (Looks like this could be gotten around by using the
% undocumented low level functions...)
%
% If data is a string will read file of that name
% If data is structure will make a stack of images
% T/Q/U/Tvar/Qvar/Uvar/QUcovar
% If data is array of structures will write element n to file
%
% m is map defn structure such that axes can be specified using
% standard fits keywords
%
% e.g:
% fitswrite(cat(3,map.T,map.Q,map.U),'test.fits',m)
% x=fitsread('test.fits');
% imagesc(x(:,:,1));

%FITSWRITE saves a Matlab matrix as a FITS image
%
%Usage: fitswrite(data, filename)
%
%Known deficiencies
%
%1. Doesn't support arrays of imaginary numbers.
%2. Only handles simple 2 dimensional arrays of data.
%
%Author: R. G. Abraham, Institute of Astronomy, Cambridge University
%        abraham@ast.cam.ac.uk

% if data is character string read from file
if(ischar(data))
  fn=data;
  load(data);
  if(exist('ac','var'))
    map=make_map(ac,m);
    if(~isfield(coaddopt,'ukpv_applied'))
      map=cal_coadd_maps(map,get_ukpervolt);
    end
  end
  if(~exist('n','var'))
    n=1;
  end
  data=map(n);
end

% if data is structure make stack of maps
if(isstruct(data))
  mapstack=true;
  data=cat(3,data.T,data.Q,data.U,data.Tvar,data.Qvar,data.Uvar,data.QUcovar);
else
  mapstack=false;
end

% seems to be needed to get correct shape when read back with
% standard fitsread function
data=permute(data,[2,1,3]);

nrow=size(data,1);
ncol=size(data,2);
nlayer=size(data,3);

header_cards = make_card('SIMPLE','T');

header_cards = [header_cards; make_card('BITPIX',-32)];     
  
if(nlayer==1)
  header_cards = [header_cards; make_card('NAXIS',2)];
else
  header_cards = [header_cards; make_card('NAXIS',3)];
end

header_cards = [header_cards; make_card('NAXIS1',nrow)];
header_cards = [header_cards; make_card('NAXIS2',ncol)];
if(nlayer>1)
  header_cards = [header_cards; make_card('NAXIS3',nlayer)];
end

% if map defn struct m provided add extra keywords to define the
% axes
if(exist('m','var'))
  header_cards = [header_cards; make_card('CRPIX1',0)];
  header_cards = [header_cards; make_card('CRPIX2',0)];
  header_cards = [header_cards; make_card('CRVAL1',m.x_tic(1))];
  header_cards = [header_cards; make_card('CRVAL2',m.y_tic(1))];
  header_cards = [header_cards; make_card('CDELT1',m.x_tic(2)-m.x_tic(1))];
  header_cards = [header_cards; make_card('CDELT2',m.y_tic(2)-m.y_tic(1))];
  header_cards = [header_cards; make_card('CTYPE1','RA')];
  header_cards = [header_cards; make_card('CTYPE2','DEC')];
  header_cards = [header_cards; make_card('POLCCONV','IAU')];
end

header_cards = [header_cards; make_card('BSCALE',1.0)];
header_cards = [header_cards; make_card('BZERO',0.0)];
header_cards = [header_cards; make_card('CREATOR','Matlab')];
header_cards = [header_cards; make_card('DATE',date)];

if(exist('fn','var'))
  [pathstr,fname]=fileparts(fn);
  header_cards = [header_cards; make_card('COMMENT',...
                 sprintf('Original filename %s',fname))];
end

% If these seem to be B2 maps add some comments
if(mapstack)
  header_cards = [header_cards; ...
                 make_card('COMMENT','BICEP2 maps in RA/Dec pixelization')];
  header_cards = [header_cards; ...
                 make_card('COMMENT','NB: IAU polarization convention!')];
  header_cards = [header_cards; ...
                 make_card('COMMENT','Stack of maps T/Q/U/Tvar/Qvar/Uvar/QUcovar')];
end
  
header_cards = [header_cards; make_card('END')];

header_record = make_header_record(header_cards);
%[ncards,dummy]=size(header_cards);
%fprintf(header_record(1,:));

fid=fopen(filename,'w');
fwrite(fid,header_record','char');

% try to figure out if we need to swap bytes. This is
% imperfect as I don't know the endian-ness of each
% architecture, so I'm only looking for ones I know for 
% sure are big-endian.
friend = computer;
if strmatch(friend,'PCWIN')
   bswap = 'b';
elseif strmatch(friend,'LNX86')
   bswap = 'b';   
elseif strmatch(friend,'ALPHA')
   bswap = 'b';
else
   bswap = 'b';
end

% write the data
fwrite(fid,data,'float',bswap);

% pad to integer * 2880
sz=numel(data)*4;
addn=2880-rem(sz,2880);
pad=ones(addn/4,1);
fwrite(fid,pad, 'float', bswap);

fclose(fid);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function card=make_card(keyword,value)
%MAKE_CARD turns a set of strings into a valid FITS card

%Make keyword field 8 bytes long
lk=length(keyword);
if (lk > 8) & (nargin>1)
	error('Keyword must be less than or equal to 8 characters!')
elseif (lk < 8 )
	keyword=[keyword,setstr(ones(1,8-lk)*32)];
end;

%Deal with both straight keyword and keyword/value pair
if (nargin==1) | strcmp(keyword,'END     ')
	%Keyword without a value
	card=keyword;	
else
	%Key/value pair has an equal sign and space at bytes 9 and 10
	card=[keyword,'= '];

	%Now output the value. The FITS standard wants things to start 
	%in different columns depending on what type of data the
	%value holds, according to the following rules:
	%
	%  Logical: T or F in column 30
	%
	%  Character string: A beginning quote in column 11 and an
	%  ending quote between columns 20 and 80.
	%
	%  Real part of an integer or floating point number: right 
	%  justified, ending in column 30.
	%
	%  Imaginary part: right justified, ending in
	%  column 50, starting after column 30 (NB. I won't bother 
	%  to support an imaginary part in this M-file, and will 
	%  let some radio astronomer who needs it add it if they want).

	if isstr(value)
  	    %Test for logical. If logical it goes in column 30 
		if (length(value)==1) & (strmatch(upper(value),'T') | strmatch(upper(value),'F'))
 			card=[card,setstr(ones(1,19)*32),value];	
		else	
			%Value must be a character string. Pad if less than 8
			%characters long.
			lv=length(value);
		    if (lv > 70)
		   error('Value must be less than 70 characters long!')
		    elseif (lv < 10 )
 	 	   value=[value,setstr(ones(1,8-lv)*32)];
		    end;
			card=[card,'''',value,''''];
		end;	
	else
		%Value must be a number. Convert to a string. Maximum
		%precision is set to 10 digits
		value=num2str(value,10);
		lv=length(value);
	
		%Write this out so it will end on column 30
		card=[card,setstr(ones(1,20-lv)*32),value];	
	end;
end;

%Now pad the output to make it exactly 80 bytes long
card=[card,setstr(ones(1,80-length(card))*32)];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hrec=make_header_record(card_matrix)

[nrow,ncol] = size(card_matrix);
n_blanks = 36 - rem(nrow,36);
blank_line = setstr(ones(1,80)*32);
hrec = [card_matrix; repmat(blank_line,n_blanks,1)];
