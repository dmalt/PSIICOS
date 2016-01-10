function [ A, key,Coh ] = CrossSpectralTimeseries(X,bInducedOnly)

if(nargin==1)
    bInducedOnly = false;
end;

[Nch, Ns  Ntr] = size(X);

if(bInducedOnly)
    ERP = mean(X,3);
    X = X-bsxfun(@minus, X,ERP);
end;
    
Xfft = fft(X,[],2);
h  = zeros(1,Ns); % nx1 for nonempty. 0x0 for empty.
if Ns > 0 && 2*fix(Ns/2) == Ns
  % even and nonempty
  h([1 Ns/2+1]) = 1;
  h(2:Ns/2) = 2;
elseif Ns>0
  % odd and nonempty
  h(1) = 1;
  h(2:(Ns+1)/2) = 2;
end
HF = repmat(h,[Nch,1,Ntr]);
XH = ifft(Xfft.*HF,[],2);
Xph = XH; %./(abs(XH)+0.0001*mean(abs(XH(:))));

clear XH;
A = zeros(Nch*(Nch-1)/2,Ns);
Nch_2 = Nch/2;


XphConj = conj(Xph);
Ntr = size(Xph,3);
range = 1:Nch-1;
trs = 1:Ntr;
k = 1;
KEY = reshape(1:Nch*Nch,Nch,Nch);
% we will take the diagonal as well
A = zeros(Nch*(Nch+1)/2,Ns);
fprintf('Calculating vectorised form of the cross spectral matrix upper triangle ... \n');
fprintf('Reference sensor (max %d): ',Nch); 
Coh = 1;

for i=1:Nch
    mn = (mean( bsxfun(@times,(Xph(1:Nch,:,:)),XphConj(i,:,:)),3));
    A(k:k+Nch-1,:) = mn;
    key(k:k+Nch-1) = KEY(1:Nch,i);
    k = k+Nch;
    if i>1
      for j=0:log10(i-1)
          fprintf('\b'); % delete previous counter display
      end
     end
     fprintf('%d', i);
end;
fprintf('\n');






