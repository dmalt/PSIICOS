% Test for different methods of finding optimal orientations
% for dipoles in tangential plane. Gross in his article
% "Dynamic imaging of coherent sources: Studying neural
% interactions in the human brain" takes the biggest singular
% value and in our case it seems wrong because orientations 
% should be real vectors and the biggest singular value 
% for orientations is delievered only by first left and right 
% singular vectors which in general are complex. 
% The latter should hold since first left and right singular vector in
% the real case correspond to max power orientations of dipoles
% (see Golub, Van Loan)
% So here we test three different strategies for calculating 
% this coherence between to sources estimating the real 
% orientations.

import ups.OrientFunctional
 
s = 2;

a = rand(s) - ones(s) * 0.5 + 1i * (rand(s) - ones(s) * 0.5);
b = [real(a), imag(a)]; % horcat
c = [real(a); imag(a)]; % vertcat

[ua,sa,va] = svd(a);
[ub,sb,vb] = svd(b);
[uc,sc,vc] = svd(c);

sb(1,1);
sc(1,1);

gross_number = sa(1,1);
fprintf('gross_number = %f\n', gross_number); 


oss_number = norm(ub(:,1)' * a * vc(:,1));
fprintf('oss_number = %f\n', oss_number); 


OPTIONS.TolX = 1e-7;
OPTIONS.TolFun = 1e-7;

[phi_max,f] = fminsearch(@(phi) OrientFunctional(phi, a), [0,0], OPTIONS);
u = [cos(phi_max(1)); sin(phi_max(1))];
v = [cos(phi_max(2)); sin(phi_max(2))];

dmalt_number = norm(u' * a * v);
fprintf('dmalt_number = %f\n', dmalt_number); 
anglea = angle(a);
uv = [u, v];
uv_os = [ub(:,1), vc(:,1)];
