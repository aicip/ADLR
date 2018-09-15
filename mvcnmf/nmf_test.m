clear all;
clear; clc

load puremixed;
mixed = mixed';
M = 64;
N = 64;
% mixed = reshape(mixed,[M,N]);

% % read data
% load('A');
% load BANDS;
% 
% type = 5;  %TODO what's type?
c = 4;
SNR = 50;
% A = A(BANDS,[1:c]); 
% [mixed, abf] = getSynData(A, 7, 0);
% [M,N,D] = size(mixed);
% mixed = reshape(mixed,M*N,D);
% 
% % add noise
% SNR = 15; %dB
% variance = sum(mixed(:).^2)/10^(SNR/10)/M/N/D;
% n = sqrt(variance)*randn([D M*N]);
% mixed = mixed' + n;
% clear n;

% remove noise
[UU, SS, VV] = svds(mixed,c);
Lowmixed = UU'*mixed;
mixed = UU*Lowmixed;
%EM = UU'*A;

% HySime algorithm 
% estimate noise
verbose = 'on';
%noise_type = 'poisson';
%[w Rn] = estNoise(mixed,noise_type,verbose);
% Ek: eigen vector 
%[c Ek]=hysime(mixed,w,Rn,verbose);    % estimate the p
%A_vca = Ek;
%vca algorithm
c=4;
[A_vca, EndIdx] = vca(mixed,'Endmembers', c,'SNR', SNR,'verbose','on');

%sisal algorithm
% [Y,Up,my,sing_val] = dataProj(mixed,c,'proj_type','affine');
% tau=10;
% [A_est] = sisal(Y,c, 'spherize', 'yes','MM_ITERS',40, 'TAU',tau, 'verbose',2);
% A_vca = A_est;

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

% % random initialization
% idx = ceil(rand(1,c)*(M*N-1));
% Ainit = mixed(:,idx);
% sinit = zeros(c,M*N);

% PCA
% [PrinComp, meanData] = pca(mixed', 0);     
[PrinComp, pca_score] = princomp(mixed',0);
% meanData = mean(mixed(:,1:D));
meanData = mean(mixed');

% test mvcnmf
tol = 1e-6;
maxiter = 100;
T = 0.015;
showflag = 0;

% use conjugate gradient to find A can speed up the learning
[Aest, sest] = mvcnmf(mixed,Ainit,sinit,UU,PrinComp,meanData,T,tol,maxiter,showflag,2,1);

figure, 
for i=1:c
   subplot(c,2,2*i-1),
   plot(Aest(:,i),'g'); 
   if i==1 title('Estimated end-members'); end
   subplot(c,2,2*i),
   imagesc(reshape(sest(i,:,:),M,N));
   if i==1 title('Estimated abundance'); end
end

% visualize endmembers in scatterplots
showflag = 0;

d = 4;
if showflag,
	Anmf = UU'*Aest;
	figure,
	for i=1:d-1
        for j=i+1:d-1
           subplot(d-2,d-2,(i-1)*(d-2)+j-i),
           plot(Lowmixed(i,1:6:end),Lowmixed(j,1:6:end),'rx');
           hold on, plot(EM(i,:), EM(j,:),'go','markerfacecolor','g');
           plot(Anmf(i,:),Anmf(j,:),'bo','markerfacecolor','b');
        end
	end
end

        
% permute results
CRD = corrcoef([A Aest]);
DD = abs(CRD(c+1:2*c,1:c));  
perm_mtx = zeros(c,c);
aux=zeros(c,1);
for i=1:c
    [ld cd]=find(max(DD(:))==DD); 
    ld=ld(1);cd=cd(1); % in the case of more than one maximum
    perm_mtx(ld,cd)=1; 
    DD(:,cd)=aux; DD(ld,:)=aux';
end
Aest = Aest*perm_mtx;
sest = sest'*perm_mtx;
Sest = reshape(sest,[M,N,c]);
sest = sest';

% show the estimations
showflag = 1;
if showflag,
	figure, 
	for i=1:c
       subplot(c,4,4*i-3),
       plot(A(:,i),'r'); axis([0 300 0 1])
       if i==1 title('True end-members'); end
       subplot(c,4,4*i-2),
       plot(Aest(:,i),'g'); axis([0 300 0 1])
       if i==1 title('Estimated end-members'); end
       subplot(c,4,4*i-1),
       imagesc(reshape(abf(i,:),M,N));
       if i==1 title('True abundance'); end
       subplot(c,4,4*i),
       imagesc(Sest(:,:,i));
       if i==1 title('Estimated abundance'); end
	end
end


% quantitative evaluation of spectral signature and abundance

% rmse error of abundances
E_rmse = sqrt(sum(sum(((abf-sest).*(abf-sest)).^2))/(M*N*c))

% the angle between abundances
nabf = diag(abf*abf'); 
nsest = diag(sest*sest');
ang_beta = 180/pi*acos( diag(abf*sest')./sqrt(nabf.*nsest));
E_aad = mean(ang_beta.^2)^.5

% cross entropy between abundance
E_entropy = sum(abf.*log((abf+1e-9)./(sest+1e-9))) + sum(sest.*log((sest+1e-9)./(abf+1e-9)));
E_aid = mean(E_entropy.^2)^.5

% the angle between material signatures
nA = diag(A'*A);
nAest = diag(Aest'*Aest);
ang_theta = 180/pi*acos( diag(A'*Aest)./sqrt(nA.*nAest) );
E_sad = mean(ang_theta.^2)^.5

% the spectral information divergence
pA = A./(repmat(sum(A),[length(A(:,1)) 1]));
qA = Aest./(repmat(sum(Aest),[length(A(:,1)) 1])); 
qA = abs(qA); 
SID = sum(pA.*log((pA+1e-9)./(qA+1e-9))) + sum(qA.*log((qA+1e-9)./(pA+1e-9)));
E_sid = mean(SID.^2)^.5

 