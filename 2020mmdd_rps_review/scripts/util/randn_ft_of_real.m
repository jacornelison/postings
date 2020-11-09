function rnd=randn_ft_of_real(N)
% rnd=randn_ft_of_real(N)
%
% Make an array of normally distributed unit variance numbers
% which have a purely real Fourier transform

% Generate arrays of complex random numbers with
% unit power spectrum ie. whose mean(abs^2)=1
m=N-1; n=N/2-1;
a=complex(randn(m,n),randn(m,n))/sqrt(2); % main block
b=complex(randn(n,1),randn(n,1))/sqrt(2); % half col
c=complex(randn(n,1),randn(n,1))/sqrt(2); % half col
d=complex(randn(1,n),randn(1,n))/sqrt(2); % half row

% Fill N*N array appropriately
% Paper and pencil required to figure out this is right...
rnd=zeros(N);
rnd(2:N,2:N/2)=a;
rnd(2:N,N/2+2:N)=rot90(conj(a),2);

rnd(2:N/2,1)=b;
rnd(N/2+2:N,1)=flipud(conj(b));

rnd(2:N/2,N/2+1)=c;
rnd(N/2+2:N,N/2+1)=flipud(conj(c));

rnd(1,2:N/2)=d;
rnd(1,N/2+2:N)=fliplr(conj(d));

rnd(1,1)=randn; rnd(N/2+1,1)=randn;
rnd(1,N/2+1)=randn; rnd(N/2+1,N/2+1)=randn;
return
