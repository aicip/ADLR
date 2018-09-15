function [Aest,sest]=do_nmfdecomp(mixed,c,M,N)

% remove noise
[UU, SS, VV] = svds(mixed,c);
Lowmixed = UU'*mixed;
mixed = UU*Lowmixed;

% HySime algorithm 
verbose = 'on';
[A_vca, EndIdx] = vca(mixed,'Endmembers', c,'verbose','on');

% FCLS
warning off;
AA = [1e-5*A_vca;ones(1,length(A_vca(1,:)))];
s_fcls = zeros(length(A_vca(1,:)),M*N);
for j=1:M*N
    r = [1e-5*mixed(:,j); 1];
    % s_fcls(:,j) = nnls(AA,r);
    s_fcls(:,j) = lsqnonneg(AA,r);
end

% use vca to initiate
Ainit = A_vca;
sinit = s_fcls;

% PCA
% [PrinComp, meanData] = pca(mixed', 0);     
[PrinComp, pca_score] = princomp(mixed',0);
% meanData = mean(mixed(:,1:D));
meanData = mean(mixed'); 

% test mvcnmf
tol = 1e-6;
maxiter = 10;%100;
T = 0.015;
showflag = 0;

% use conjugate gradient to find A can speed up the learning
[Aest, sest] = mvcnmf_new_me(mixed,Ainit,sinit,UU,PrinComp,meanData,T,tol,maxiter,showflag,2,1);

% show results
do_nmfresultdisp(c,Aest,sest,showflag,UU,Lowmixed,M,N);

end
