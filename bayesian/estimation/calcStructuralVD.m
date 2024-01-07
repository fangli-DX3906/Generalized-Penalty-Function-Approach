function VD = calcStructuralVD(n, H, IRF)


IRsq = reshape(IRF.^2, n, n, H+1);
IRsqsumH = cumsum(IRsq,3);
totFEV = squeeze(sum(IRsqsumH,2));
for ih = 1:size(totFEV,2)
    vd(:,:,ih) = IRsqsumH(:,:,ih)./totFEV(:,ih);
end
VD = reshape(vd, n^2, H+1);

% n   = size(sigma,1);
% J  = [eye(n,n) zeros(n,n*(p-1))];
% TH1 = J*A^0*J'; 
% TH = TH1*chol(sigma)' * rotation; 
% TH = TH'; 
% TH2 = (TH.*TH); 
% TH3(1,:,:) = TH2;
% for i = 2:H
%     TH = J*A^(i-1)*J'*chol(sigma)' * rotation; 
%     TH = TH'; 
%     TH2  = (TH.*TH); 
%     TH3(i,:,:) = squeeze(TH3(i-1,:,:)) + TH2;
% end
% TH4 = sum(TH3,2);
% % % 1: h, 2: shocks, 3: variables
% VC  = zeros(H,n,n);
% for j = 1:n
%     VC(:,j,:) = TH3(:,j,:)./TH4;
% end
% 
% VD = zeros(n^2, H);
% for v = 1:n
%     for s = 1:n
%         id = n*(v-1)+s;
%         VD(id, :) = VC(:, s, v);
%     end
% end

end

