function [Co,X,Y,Uc,Vc] = subcorr(A,B,NOSVD,Sa,Va,Sb,Vb)
%SUBCORR Calculate the subspace correlation between two matrices
% function Co = subcorr(A,B);
%  or
% function [Co,X,Y,Uc,Vc] = subcorr(A,B,NOSVD,Sa,Va,Sb,Vb);
%  or
% function [Co,X,Y] = subcorr(A,B,C);
% 
% return the correlations between the two subspaces defined as the 
% ranges of A and B, A is m x p, B is m x q
%
% Co is the vector of correlations, 0 <= s(min(p,q)) <= ... <= s(1)
% X and Y are the associated translators, such that
% A*X(:,1) is most correlated with B*Y(:,1), and both
%  resultant vectors are unit norm.  The matrix A*X is orthogonal, as
%  is B*Y, and the corresponding pairs are the principal vectors.
%
% If (optional) NOSVD is set to 1, then A and B are assumed to already be
% orthogonal matrices, A'*A = I and B'*B = I.  Default is 0.
% In case of rank deficient matrices, the full rank subspace is used.
% If X and Y are desired, the user must provide the rest of the decompositions
%  of A and B in the inputs Sa,Va,Sb,Vb.  Sa and Sb may be vectors.
% 
% If (optional) C is a matrix with the same number of rows as A and B,
% then it is used as a constraint between the two spaces.  The two spaces
% are projected away from the constraint.
%
% SEE SUBANGLE, SUBSPACE
% from Golub and Van Loan, page 584

% Copyright (c) 1990-1995, The Regents of the University of California.
% This software was produced under a U.S. Government contract
% (W-7405-ENG-36) by Los Alamos National Laboratory, which is operated
% by the University of California for the U.S. Department of Energy,
% and was funded in part by NIH grant R01-MH53213 through the University
% of Southern California to Los Alamos National Laboratory, 
% and was funded in part by NIH grant R01-EY08610 to Los Alamos
% National Laboratory.
% The U.S. Government is licensed to use, reproduce, and distribute this
% software.  Permission is granted to the public to copy and use this
% software without charge, provided that this Notice and any statement
% of authorship are reproduced on all copies.  Neither the Government
% nor the University makes any warranty, express or implied, or assumes
% any liability or responsibility for the use of this software.
%
% Author: John C. Mosher, Ph.D.
% Los Alamos National Laboratory
% Group ESA-MT, MS J580
% Los Alamos, NM 87545
% email: mosher@LANL.Gov

%  John C. Mosher  8/7/90
% JCM modified 2/3/92 to do svd(.,0) for faster computations
% JCM modified 7/19/93 for NOSVD flag.
% JCM 7/12/95 handle zero conditions
% JCM 1/30/96 Constraint option
% JCM 7/19/96 Derived SUBCORR from SUBANGLE to return X and Y, not u and v
%  and return correlations, not subspace angles.


[mA,p] = size(A);
[mB,q] = size(B);
if(mA~=mB),error('A and B must have same number of rows'),end

if(exist('Sa', 'var') == 1),		% user provided decompositions
  if(min(size(Sa)) > 1),	% user gave diagonal matrix
    Sa = diag(Sa);
  end
  % assume  Sb must exist as well
  if(min(size(Sb)) > 1),	% user gave diagonal matrix
    Sb = diag(Sb);
  end
end

if(exist('NOSVD', 'var') ~= 1),	% user gave nothing
  NOSVD = [];			% make empty
  NOSVD_FLAG = 0; 		% default, calculate SVD of A and B
end

if(isempty(NOSVD)),		% user gave blank constraint
  NOSVD_FLAG = 0;		% need to svd
  Constraint = [];
end

if(size(NOSVD,1) == size(A,1)),
  Constraint = NOSVD;			% constraint
  NOSVD_FLAG = 0;		% need SVDs
  inv_Constraint = pinv(Constraint);		% pseudo-inverse
end

if(~NOSVD_FLAG),  % need to calculate the svd
  if(~isempty(Constraint)),		% user gave a constraint
    % we have determined the rank of both matrices.  The constraint sets up
    % the following question: what are the best subspace angles of the pair
    % [A C] and B, where C is already presumed to span a subspace within B.
    % If C were completely contained within B, the numerical problems would
    % be easier; however, C may only "mostly" lie in B.  As such, different
    % linear combinations of A may combine with C to better perfect the
    % primary angles between C and B.  If we project away the space spanned
    % by C from both sides, we can focus on the contribution that some
    % linear combination of A brings to the model between C and B.
    % Note that the althernative is to simply call this routine as
    % subcorr([A C],B), which will answer all the ways that A may combine.
    Sa = svd(A);		% just the singular values
    Sb = svd(B);
    tolA = mA * Sa(1) * eps;	% numerical tolerance for rank
    tolB = mB * Sb(1) * eps;
    rA = sum(Sa > tolA); 	% rank of matrix A
    rB = sum(Sb > tolB); 	% rank of matrix B

    A = A - Constraint*(inv_Constraint*A); 	% project away the constraint
    B = B - Constraint*(inv_Constraint*B);
    % full Svd these projected data sets
    if(mA >= p),
      [Ua,Sa,Va] = svd(A,0); 
    else
      [Va,Sa,Ua] = svd(A',0);
    end
    if(length(Sa) >1),Sa = diag(Sa);end
    if(mB >= q),
      [Ub,Sb,Vb] = svd(B,0); 
    else
      [Vb,Sb,Ub] = svd(B',0);
    end
    if(length(Sb) >1),Sb = diag(Sb);end
    rB = rB - size(Constraint,2);	% truncate the projected away space
    % inverters
    if(nargout > 1), 		% user wants translators
      VaSai = Va(:,1:rA)*diag(1../Sa(1:rA));
      VbSbi = Vb(:,1:rB)*diag(1../Sb(1:rB));
    end
    
  else				% user gave no constraint

    if(mA >= p),
      [Ua,Sa,Va] = svd(A,0); 
    else
      [Va,Sa,Ua] = svd(A',0);
    end
    if(length(Sa) >1),Sa = diag(Sa);end
    if(mB >= q),
      [Ub,Sb,Vb] = svd(B,0); 
    else
      [Vb,Sb,Ub] = svd(B',0);
    end
    if(length(Sb) >1),Sb = diag(Sb);end
    tolA = mA * Sa(1) * eps;
    tolB = mB * Sb(1) * eps;
    rA = sum(Sa > tolA); 	% rank of matrix A
    rB = sum(Sb > tolB); 	% rank of matrix B
    % inverters
    if(nargout > 1), 		% user wants translators
      VaSai = Va(:,1:rA)*diag(1../Sa(1:rA));
      VbSbi = Vb(:,1:rB)*diag(1../Sb(1:rB));
    end

  end				% constraint or no constraint

else 				% user has already calculated the SVDs

  Ua = A;
  Ub = B;
  rA = p;  % width of A, which is assumed to be Ua
  rB = q;
  if(nargout > 1),		% user wants translators
    VaSai = Va(:,1:rA)*diag(1../Sa(1:rA));
    VbSbi = Vb(:,1:rB)*diag(1../Sb(1:rB));
  end

end

% bogus condition, probably a null matrix given
if(rA == 0 || rB == 0), 		% one of these has no rank
  Co = NaN; 			% not defined
  return
end

C = Ua(:,1:rA)'*Ub(:,1:rB);  % cross product of subspaces

if(nargout < 2)  % want only angles
  Co = svd(C);
  too_big = find(Co > 1);	% numerical error
  Co(too_big) = ones(length(too_big),1);
  too_small = find(Co < 0);
  Co(too_small) = zeros(length(too_small),1);
  return
else % nargout is 2 or more
  [Uc,Co,Vc] = svd(C);
  if(min(size(Co))==1),Co = Co(1);else Co = diag(Co);end
  too_big = find(Co > 1);	% numerical error
  Co(too_big) = ones(length(too_big),1);
  too_small = find(Co < 0);	% numerical error
  Co(too_small) = zeros(length(too_small),1);
  X = VaSai*Uc;
  Y = VbSbi*Vc;
  return
end
