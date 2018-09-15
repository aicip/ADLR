function [] = do_nmfresultdisp(c,Aest,sest,showflag,UU,Lowmixed,M,N)

% figure, 
% for i=1:c
%    subplot(c,2,2*i-1),
%    plot(Aest(:,i),'g'); 
%    if i==1 title('Estimated end-members'); end
%    subplot(c,2,2*i),
%    imagesc(reshape(sest(i,:,:),M,N));
%    if i==1 title('Estimated abundance'); end
% end

figure, 
for i=1:round(c^0.5)
    for j=1:round(c^0.5)
        idx = (i-1)*round(c^0.5)+j;
        if(idx<=c)
            subplot(round(c^0.5),round(c^0.5),idx);
            imagesc(reshape(sest(idx,:,:),M,N));
        end
    end
end

% visualize endmembers in scatterplots
d = 4;
if showflag,
	Anmf = UU'*Aest;
	figure,
	for i=1:d-1
        for j=i+1:d-1
           subplot(d-2,d-2,(i-1)*(d-2)+j-i),
           plot(Lowmixed(i,1:6:end),Lowmixed(j,1:6:end),'rx');
           %hold on, plot(EM(i,:), EM(j,:),'go','markerfacecolor','g');
           plot(Anmf(i,:),Anmf(j,:),'bo','markerfacecolor','b');
        end
	end
end

        
% % permute results
% CRD = corrcoef([A Aest]);
% DD = abs(CRD(c+1:2*c,1:c));  
% perm_mtx = zeros(c,c);
% aux=zeros(c,1);
% for i=1:c
%     [ld cd]=find(max(DD(:))==DD); 
%     ld=ld(1);cd=cd(1); % in the case of more than one maximum
%     perm_mtx(ld,cd)=1; 
%     DD(:,cd)=aux; DD(ld,:)=aux';
% end
% Aest = Aest*perm_mtx;
% sest = sest'*perm_mtx;
% Sest = reshape(sest,[M,N,c]);
% sest = sest';

% % show the estimations
% showflag = 1;
% if showflag,
% 	figure, 
% 	for i=1:c
%        subplot(c,4,4*i-3),
%        plot(A(:,i),'r'); axis([0 300 0 1])
%        if i==1 title('True end-members'); end
%        subplot(c,4,4*i-2),
%        plot(Aest(:,i),'g'); axis([0 300 0 1])
%        if i==1 title('Estimated end-members'); end
%        subplot(c,4,4*i-1),
%        imagesc(reshape(abf(i,:),M,N));
%        if i==1 title('True abundance'); end
%        subplot(c,4,4*i),
%        imagesc(Sest(:,:,i));
%        if i==1 title('Estimated abundance'); end
% 	end
% end
% 
% 
% % quantitative evaluation of spectral signature and abundance
% 
% % rmse error of abundances
% E_rmse = sqrt(sum(sum(((abf-sest).*(abf-sest)).^2))/(M*N*c));
% 
% % the angle between abundances
% nabf = diag(abf*abf'); 
% nsest = diag(sest*sest');
% ang_beta = 180/pi*acos( diag(abf*sest')./sqrt(nabf.*nsest));
% E_aad = mean(ang_beta.^2)^.5;
% 
% % cross entropy between abundance
% E_entropy = sum(abf.*log((abf+1e-9)./(sest+1e-9))) + sum(sest.*log((sest+1e-9)./(abf+1e-9)));
% E_aid = mean(E_entropy.^2)^.5;
% 
% % the angle between material signatures
% nA = diag(A'*A);
% nAest = diag(Aest'*Aest);
% ang_theta = 180/pi*acos( diag(A'*Aest)./sqrt(nA.*nAest) );
% E_sad = mean(ang_theta.^2)^.5;
% 
% % the spectral information divergence
% pA = A./(repmat(sum(A),[length(A(:,1)) 1]));
% qA = Aest./(repmat(sum(Aest),[length(A(:,1)) 1])); 
% qA = abs(qA); 
% SID = sum(pA.*log((pA+1e-9)./(qA+1e-9))) + sum(qA.*log((qA+1e-9)./(pA+1e-9)));
% E_sid = mean(SID.^2)^.5;

end
