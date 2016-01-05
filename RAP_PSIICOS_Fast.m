function [ indep_topo, c_ss_hat, PVU, SubC,INDrap, Cp, Upwr] = RAP_PSIICOS_Fast( C, G2dU,RAPIts, Rnk, Upwr)
%Input:
% Power and Shift Independent Imaging of Coherent Sources (PSIICOS)
% in MEG data 
% G2dU - forward model matrix such that each source is served by two
% columns of this matrix corresponding to the topographies of dipoles in
% the tangential plane
% C - is a sensor space cross-spectral matrix
% Rnk - limits the dimension of Volume Conduction subspace
%%Output:
% Cs - the upper diagonal of source space cross-spectrum
% IND - vector of pair of dipole pairs corresponding to the elements of Cs
% Cp - projected away from the VC subspace sensor space cross-soectral
% matrix
% Upwr - VC subspace basis matrix. Columns of Upwr span the VC subspace

%% Preparatory steps
if(nargin<4)
    Rnk  = 350;
end;

if(nargin<5)
    Upwr  = [];
end;

if(isempty(Rnk))
    if(isempty(Upwr))
        Rnk = 350;
    else
        Rnk = size(Upwr,2);
    end;
end;

Nsrc = size(G2dU,2)/2; % two topography columns per each source of the grid
Nch = size(G2dU,1);

%% perform projection of the coherence matrix away from the power only
if(isempty(Upwr))
    [Cpvec, Upwr] = ProjectAwayFromPowerFixedOr(C(:), G2dU,Rnk);
else % use the existing matrix if profided
    c = Upwr'*C(:);
    Cpvec  = C(:)-Upwr*c;
end;

 k =1;
 for rnk = 1:10:size(Upwr,2)
     Cimp  = imag(C(:))-Upwr(:,1:rnk)*(Upwr(:,1:rnk)'*imag(C(:)));
     nrmim(k) = norm(Cimp(:))/norm(imag(C(:)));
     k = k+1;
end;

Cp = reshape(Cpvec,size(C,1),size(C,2));

%% normalize forward matrix
% 
% for i = 1:Nsrc
%     range_i = i*2-1:i*2;
%     G2dU(:,range_i(1)) = G2dU(:,range_i(1))/norm(G2dU(:,range_i(1)));
%     G2dU(:,range_i(2)) = G2dU(:,range_i(2))/norm(G2dU(:,range_i(2)));
% end;

%% scan all pairs with efficient vectorised implementation
Cprap = Cp;
Urap = [];
range2 = 1:2;
indep_topo = zeros(prod(size(Cprap)),RAPIts*2);
c_ss_hat = zeros(1,RAPIts*2);
for rap = 1:RAPIts
    [Cs(rap,:), IND, Cs0] =PSIICOS_ScanFast(G2dU, Cprap);
    [val_max ind_max] = max(Cs(rap,:));
    pair_max = IND(ind_max,:);
    i = IND(ind_max,1); 
    j = IND(ind_max,2); 
    range_i = i*2-1:i*2;
    range_j = j*2-1:j*2;
    ai = G2dU(:,i*2-1:i*2);
    aj = G2dU(:,j*2-1:j*2);
    cs = ai'*Cprap*aj;
    [u s v] = svd(cs);
    % u and v are complex but the orientation vectors of dipoles are physical and
    % therefore can not be defined over the field of complex numbers. 
    % Do this trick to force SVD to real 2x1 vectors(first left, ten right) to find 
    % the real orientations
    csr = real(cs);
    csi = imag(cs);
    [uL sL vL] =svd([ csr csi]);
    [uR sR vR] =svd([ csr;csi]);
    ai_or = ai*uL(:,1);
    aj_or = aj*vR(:,1);
    % not that norm(ai_or'*Cprap*aj_or) = s(1,1)= s_max and therefore the
    % fast scan implemented is valid
    qij = ai_or*aj_or';
    qji = aj_or*ai_or';
    qijp = qij(:) - Upwr*(Upwr'*qij(:));
    qjip = qji(:) - Upwr*(Upwr'*qji(:));
    qp = [qijp,qjip];
    c_ss_hat(1,range2) = (pinv(qp)*Cprap(:))';
    SubC(rap) = subcorr(Cprap(:),qp);
    Cprap_vec = Cprap(:)-qp*(pinv(qp)*Cprap(:));
    Cprap = reshape(Cprap_vec,size(Cprap));
    indep_topo(:,range2) = qp;
    PVU(rap) = norm(Cp(:) - indep_topo*c_ss_hat')/norm(Cp(:));
    INDrap(rap,1) = i; 
    INDrap(rap,2) = j; 
    range2 = range2+2;
  end

    

