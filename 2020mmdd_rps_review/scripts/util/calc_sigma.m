% function to calculate the +1 sigma and -1sigma away from the peak of a 
% normalized likelihood function
% returns index of +sig and -sig as [is1, is2]
% unlike calc_sigma, this is done by taking contours of constant likelihood (starting from the peak and iteratively
% finding the contour of constand likelyhood that  contains around the peak 68.27% of the likelhood.
function chi=calc_sigma(chi)
  
  for k=1:size(chi,1)
    for j=1:size(chi,2)
      
      if size(chi(k,j).dof,2) == 0
        continue
      end
      
      % renormalize the likelihood so the whole sum sums to 1.
      chi(k,j).rlik=chi(k,j).rlik/sum(chi(k,j).rlik,2);
      
      x=chi(k,j).r;
      % get index of the max of the likelihood
      [m,im]=max(chi(k,j).rlik);
      
      %define 2sigma=68.27%
      sigma=0.6827;
      
      % find index of contour of constant likelihood and interate the PDF under that curve
      new_contour = m;
      dr = 0.00001
      c=0;
      % repeat until we found countour that encompasses 68% of likelhood curve
      while  c <= sigma
        new_contour =  new_contour - 0.001;
        q=find(abs(chi(k,j).rlik - new_contour) < dr);
        
        %change the size of the interval untilw e get only 2 points
        while (size(q,2) ~=2 && size(q,2) ~=3 && size(q,2) ~=4)
          if size(q,2) < 2
            inc = +0.00001
          elseif size(q,2) > 4 
          inc = -0.00001
        end
          dr = dr + inc;
          q=find(abs(chi(k,j).rlik - new_contour) < dr);
        end
        
        is1=q(1);is2=q(end); 
        c=sum(chi(k,j).rlik(is1:is2)) ;
      end
      
      sprintf('peak at %1.3f. +%1.3f -%1.3f',x(im), -x(is1)+x(im),-x(im)+x(is2))
      %save m, +1sig, -1sig
      chi(k,j).irmax=im;
      chi(k,j).isigp=is1
      chi(k,j).isigm=is2;
      chi(k,j).rmax=x(im);
      chi(k,j).sigp=abs(x(is1)-x(im));
      chi(k,j).sigm=abs(x(im)-x(is2));

      %repeat for all sims
      for i=1:size(chi(k,j).slik,2)
        
        % renormalize the likelihood so the whole sum sums to 1.
        chi(k,j).slik(:,i)=chi(k,j).slik(:,i)./sum(chi(k,j).slik(:,i),1);
        
        [m,im]=max(chi(k,j).slik(:,i));
        
        % find index of contour of constant likelihood and interate the PDF under that curve
        new_contour = m;
        c=0;
        dr=0.0008;
      
        % repeat until we found countour that encompasses 68% of likelhood curve
        
        while  c < sigma
          new_contour =  new_contour - 0.001;
          q=find(abs(chi(k,j).slik(:,i) - new_contour) < dr);
          
          %change the size of the interval until we get only 2 points
          while (size(q,1) ~=2 && size(q,1) ~=3 && size(q,1) ~=4)
            if size(q,1) < 2
              inc = +0.00001;
            elseif size(q,1) > 4
              inc = -0.00001;
            end
            dr = dr + inc;
            q=find(abs(chi(k,j).slik(:,i) - new_contour) < dr);
            
          end
        
          is1=q(1);is2=q(end); 
          c=sum(chi(k,j).slik(is1:is2,i)) ;
        end
        
        %save m, +1sig, -1sig
        chi(k,j).srmax(i)=x(im);
        chi(k,j).ssigp(i)=x(is1)-x(im);
        chi(k,j).ssigm(i)=x(im)-x(is2);
      end
      
    end
  end
  
return