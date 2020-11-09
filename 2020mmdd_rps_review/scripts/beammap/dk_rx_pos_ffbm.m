function pos = dk_rx_pos_ffbm(rxNum,dkangle)
% function pos = dk_rx_pos_ffbm(rxNum,dkangle)
%
% Finds the receiver position for our standard far field beam maps
% given the dk angle or schedule letter
% Position 3 is the position at base of drum, centered
%
% map  Dk        pos1 pos2 pos3 pos4 pos5 pos6 pos7 pos8 pos9 pos10
% a    -122      rx1       rx0       rx4       rx3       rx2
% b    -50       rx2       rx1       rx0       rx4       rx3
% c    +22       rx3       rx2       rx1       rx0       rx4
% d    +94       rx4       rx3       rx2       rx1       rx0
% e    +166      rx0       rx4       rx3       rx2       rx1
% f    -86            rx1       rx0       rx4       rx3       rx2
% g    -14            rx2       rx1       rx0       rx4       rx3
% h    +58            rx3       rx2       rx1       rx0       rx4 
% i    +130           rx4       rx3       rx2       rx1       rx0 
% j    +202/-158      rx0       rx4       rx3       rx2       rx1 
% 
% mask = ffbm_maskMaker() gives position cuts for 10 positions

for ii = 1:length(dkangle)
  
  % Allow schedule position input also (a-j)
  if ischar(dkangle(ii))
    switch dkangle(ii)
      case 'a'
	dk = -122;
      case 'b'
	dk = -50;
      case 'c'
	dk = 22;
      case 'd'
	dk = 94;
      case 'e'
	dk = 166;
      case 'f'
	dk = -86;
      case 'g'
	dk = -14;
      case 'h'
	dk = 58;
      case 'i'
	dk = 130;
      case 'j'
	dk = -158;
    end
  else
    dk = dkangle(ii);
    dk = round(dk);
  end

  if dk < -180
    dk = dk + 360;
  elseif dk > 180
    dk = 360 - dk;
  end 

  rxlist = [4 3 2 1 0];
  rxposlist_A = [1 3 5 7 9];
  rxposlist_B = [2 4 6 8 10];

  switch dk
    case -122 
      rxlist = circshift(rxlist,[0,-3]);
      rxpos =  rxposlist_A(find(rxlist==rxNum));
    case -50
      rxlist = circshift(rxlist,[0,-2]);
      rxpos = rxposlist_A(find(rxlist==rxNum));
    case 22
      rxlist = circshift(rxlist,[0,-1]);
      rxpos = rxposlist_A(find(rxlist==rxNum));
    case 94
      rxlist=circshift(rxlist,[0,0]);
      rxpos=rxposlist_A(find(rxlist==rxNum));
    case 166
      rxlist=circshift(rxlist,[0,1]);
      rxpos=rxposlist_A(find(rxlist==rxNum));
    case -86
      rxlist=circshift(rxlist,[0,-3]);
      rxpos=rxposlist_B(find(rxlist==rxNum));
    case -14
      rxlist=circshift(rxlist,[0,-2]);
      rxpos=rxposlist_B(find(rxlist==rxNum));
    case 58
      rxlist=circshift(rxlist,[0,-1]);
      rxpos=rxposlist_B(find(rxlist==rxNum));
    case 130
      rxlist=circshift(rxlist,[0,0]);
      rxpos=rxposlist_B(find(rxlist==rxNum));
    case {-158,158} % dk 202: 360 - 202 = 158
      rxlist=circshift(rxlist,[0,1]);
      rxpos=rxposlist_B(find(rxlist==rxNum));
  end
  
  pos(ii)=rxpos;
end
      
return
  
