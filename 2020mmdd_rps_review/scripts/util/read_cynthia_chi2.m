function ra=read_cynthia_chi2(i,rf,ra,old)
% function ra=read_cynthia_chi2(i,rf,ra,old)
  
%function to read cynthia's table of chi2^2 and pte.
% for new format, the name should be something like cynthia_chi2_poly3_5bins
%
%specifying the old keyword  switches to reading in the old format
%original format  of chi2 tables such as  in cynthia_jack_chi2_080519.txt
%or  cynthia_jackknife_spectra_0_080520_II.dat  is obsolete now.

  switch old
   case 0
    % read new format order is TT,EE,BB, TE, TB, EB
    filejack=strcat(rf{i},'_jack.txt');
    filelcdm=strcat(rf{i},'_lcdm.txt');
    fid=fopen(sprintf('final/chi2/%s',filelcdm));
    C0=textscan(fid, '%s%n%n%n %n%n%n %n%n%n','commentStyle','#');
    fclose(fid);
    fid=fopen(sprintf('final/chi2/%s',filejack));
    C=textscan(fid, '%s%n%n%n %n%n%n %n%n%n','commentStyle','#');
    fclose(fid);
    sim100100=reshape(([C0{5}',C{5}']),6,6)';
    sim150150=reshape(([C0{6}',C{6}']),6,6)';
    sim100150=reshape(([C0{7}',C{7}']),6,6)';
    teo100100=reshape(([C0{8}',C{8}']),6,6)';
    teo150150=reshape(([C0{9}',C{9}']),6,6)';
    teo100150=reshape(([C0{10}',C{10}']),6,6)';
              
    ra{1,i}(1).ptes=sim100100(1,:);
    ra{1,i}(2).ptes=sim150150(1,:);
    ra{1,i}(3).ptes=sim100150(1,:);
    ra{1,i}(1).ptet=teo100100(1,:);
    ra{1,i}(2).ptet=teo150150(1,:);
    ra{1,i}(3).ptet=teo100150(1,:);
    
    ra{2,i}(1).ptes=sim100100(4,:);
    ra{2,i}(2).ptes=sim150150(4,:);
    ra{2,i}(3).ptes=sim100150(4,:);
    ra{2,i}(1).ptet=teo100100(4,:);
    ra{2,i}(2).ptet=teo150150(4,:);
    ra{2,i}(3).ptet=teo100150(4,:);
    
    ra{3,i}(1).ptes=sim100100(2,:);
    ra{3,i}(2).ptes=sim150150(2,:);
    ra{3,i}(3).ptes=sim100150(2,:);
    ra{3,i}(1).ptet=teo100100(2,:);
    ra{3,i}(2).ptet=teo150150(2,:);
    ra{3,i}(3).ptet=teo100150(2,:);
    
    ra{4,i}(1).ptes=sim100100(6,:);
    ra{4,i}(2).ptes=sim150150(6,:);
    ra{4,i}(3).ptes=sim100150(6,:);
    ra{4,i}(1).ptet=teo100100(6,:);
    ra{4,i}(2).ptet=teo150150(6,:);
    ra{4,i}(3).ptet=teo100150(6,:);
    
    ra{5,i}(1).ptes=sim100100(5,:);
    ra{5,i}(2).ptes=sim150150(5,:);
    ra{5,i}(3).ptes=sim100150(5,:);
    ra{5,i}(1).ptet=teo100100(5,:);
    ra{5,i}(2).ptet=teo150150(5,:);
    ra{5,i}(3).ptet=teo100150(5,:);
    
    ra{6,i}(1).ptes=sim100100(3,:);
    ra{6,i}(2).ptes=sim150150(3,:);
    ra{6,i}(3).ptes=sim100150(3,:);
    ra{6,i}(1).ptet=teo100100(3,:);
    ra{6,i}(2).ptet=teo150150(3,:);
    ra{6,i}(3).ptet=teo100150(3,:);
    
    ra{7,i}(1).ptes=NaN(1,6);
    ra{7,i}(1).ptet=NaN(1,6);
    
    % re-order to my TT,TE,EE,BB,TB,EB
    if(1)
      for j=1:size(ra,1)
        for f=1:length(ra{j,i})
          ra{j,i}(f).ptes=ra{j,i}(f).ptes([1,4,2,3,5,6]);
          ra{j,i}(f).ptet=ra{j,i}(f).ptet([1,4,2,3,5,6]);
        end
      end
    end
    
   case 1
    % read old format
    % Cynthia's PTE order is TT,TE,EE,BB
    fid=fopen(sprintf('final/chi2/%s',rf{i}));
    C=textscan(fid, '%d8%n%n%n%n%n%n%n%n','commentStyle','#');
    fclose(fid);
    c=[C{6} C{7} C{8} C{9} NaN(15,2)];
    
    ra{1,i}(1).ptes=NaN(1,6);
    ra{1,i}(2).ptes=NaN(1,6);
    ra{1,i}(3).ptes=NaN(1,6);
    
    ra{2,i}(1).ptes=c(3,:);
    ra{2,i}(2).ptes=c(3+5,:);
    ra{2,i}(3).ptes=c(3+10,:);
    
    ra{3,i}(1).ptes=c(1,:);
    ra{3,i}(2).ptes=c(1+5,:);
    ra{3,i}(3).ptes=c(1+10,:);
    
    ra{4,i}(1).ptes=c(5,:);
    ra{4,i}(2).ptes=c(5+5,:);
    ra{4,i}(3).ptes=c(5+10,:);

    ra{5,i}(1).ptes=c(4,:);
    ra{5,i}(2).ptes=c(4+5,:);
    ra{5,i}(3).ptes=c(4+10,:);
    
    ra{6,i}(1).ptes=c(2,:);
    ra{6,i}(2).ptes=c(2+5,:);
    ra{6,i}(3).ptes=c(2+10,:);
        
    ra{7,i}(1).ptes=NaN(1,6);
    
    
  end