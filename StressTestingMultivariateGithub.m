% code and data for the paper: 
% Stress testing and systemic risk measures using elliptical conditional multivariate probabilities
% https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3575512
% Tomaso Aste 03/05/2021
clear
close all

mod180 = @(x) x.*(x<=90)+(180-x).*(x>90);

load DataGithub.mat

O = cov(X0); %X0 are the log-returns

Angle = NaN(length(SectorNames));
DXY   = NaN(length(SectorNames));
TvXY  = NaN(length(SectorNames));
for s = 1:length(SectorNames)
    iX = find(strcmp(SectorNames(s),Sectors(:,3)));
    VaRX = quantile(X0(:,iX),1-q);
    OXX = O(iX,iX);
    for r = s:length(SectorNames)
        iY = find(strcmp(SectorNames(r),Sectors(:,3)));
        VaRY = quantile(X0(:,iY),1-q);
        OXY = O(iX,iY);
        OYX = O(iY,iX);
        OYY = O(iY,iY);
        % columns | impact of 
        vLYX(r,s) = mean(OYX/OXX*(VaRX-mean(X0(:,iX)))'); %Y<-X s->r %the diagonal is the VaR (OYX/OXX=I)
        vLYX(s,r) = mean(OXY/OYY*(VaRY-mean(X0(:,iY)))'); %X<-Y r->s
        uLYX(r,s) = -mean(OYX/OXX*ones(size(iX,1),1)); %Y<-X s->r %the diagonal is the VaR (OYX/OXX=I)
        uLYX(s,r) = -mean(OXY/OYY*ones(size(iY,1),1)); %X<-Y r->s
        mLYX(r,s) = min(OYX/OXX*(VaRX-mean(X0(:,iX)))'); %Y<-X s->r %the diagonal is the VaR (OYX/OXX=I)
        mLYX(s,r) = min(OXY/OYY*(VaRY-mean(X0(:,iY)))'); %X<-Y r->s
        I(s,r) = 0.5*(sum(log2(eig(OXX)))+sum(log2(eig(OYY)))-sum(log2(eig(O))));
        % R(s,r) = corr(mean(X0(:,iX),2),mean(X0(:,iY),2)) ; % to validate I 
        Ampl(r,s) = vLYX(r,s)/mean((VaRX-mean(X0(:,iX))));
        Ampl(s,r) = vLYX(s,r)/mean((VaRY-mean(X0(:,iY))));
        Reaction(r,s) = vLYX(r,s)/mean((VaRY-mean(X0(:,iY))));
        Reaction(s,r) = vLYX(s,r)/mean((VaRX-mean(X0(:,iX))));
        %% angle 
        OYY_X = OYY - OYX/OXX*OXY ;
        if min(eig(OYY_X))<=0,fprintf('OYY_X not positive defined!\n')
        elseif ~isreal(eig(OYY_X)),fprintf('OYY_X not real !\n')
        else
        [u,e]=eig(OYY);
        [m,i]=max(max(e));
        uYY = u(:,i); 
        sYY = (m); %var(uYY'*X0(:,iY)'/sum(uYY));
        TvYY  = sum(log(diag(e)));
        [u,e]=eig(OYY_X);
        [m,i]=max(max(e));
        uYY_X = u(:,i); %/sum(u(:,i));
        sYY_X = (m); %var(uYY_X'*X0(:,iY)'/sum(uYY_X));
        TvYY_X  = sum(log(diag(e)));
        Angle(r,s) = acos(abs(uYY'*uYY_X))/pi*180; %Y <- X 
        DXY(r,s) = -(sYY_X-sYY)/sYY;
        TvXY(r,s)= (TvYY_X-TvYY);
        end
        OXX_Y = OXX - OXY/OYY*OYX ;
        if min(eig(OXX_Y))<=0,fprintf('OXX_Y not positive defined!\n')
        elseif ~isreal(eig(OXX_Y)),fprintf('OXX_Y not real !\n')
        else
        [u,e]=eig(OXX);
        [m,i]=max(max(e));
        uXX = u(:,i); %/sum(u(:,i));;
        sXX = (m); %var(uXX'*X0(:,iX)'/sum(uXX));
        TvXX  = sum(log(diag(e)));
        [u,e]=eig(OXX_Y);
        [m,i]=max(max(e));
        uXX_Y = u(:,i); %/sum(u(:,i));
        sXX_Y = (m); %var(uXX_Y'*X0(:,iX)'/sum(uXX_Y));
        TvXX_Y  = sum(log(diag(e)));
        Angle(s,r) = acos(abs(uXX'*uXX_Y))/pi*180; %X <- Y  
        DXY(s,r) = -(sXX_Y-sXX)/sXX;
        TvXY(s,r)= (TvXX_Y-TvXX);        
        end
    end
    %% Mahalanobis factor
    dX(s) = (quantile(X0(:,iX),1-q)-mean(X0(:,iX)))/OXX*(quantile(X0(:,iX),1-q)-mean(X0(:,iX)))';
    pX(s) = length(iX);
    kY = setdiff([1:N],iX)';
    OYY = O(kY,kY);
    OXY = O(iX,kY);
    OXX = O(iX,iX);
    OYX = O(kY,iX);
    vImpacted(s)  =  mean(OXY/OYY*(quantile(X0(:,kY),1-q)-mean(X0(:,kY)))'); % impact on s from all others
    vImpacting(s) =  mean(OYX/OXX*(quantile(X0(:,iX),1-q)-mean(X0(:,iX)))'); % impact of s on all others
    mImpacted(s)  =  min(OXY/OYY*(quantile(X0(:,kY),1-q) -mean(X0(:,kY)))'); % impact on s from all others
    mImpacting(s) =  min(OYX/OXX*(quantile(X0(:,iX),1-q) -mean(X0(:,iX)))'); % impact of s on all others
    uImpacted(s)  =  mean(OXY/OYY*ones(size(kY,1),1)); % impact on s from all others
    uImpacting(s) =  mean(OYX/OXX*ones(size(iX,1),1)); % impact of s on all others
end
%% Average losses 
% small NEGATIVE correlation with I ~ 16%  
M=-(vLYX-diag(diag(nan(10))))'; corr(M(triu(ones(size(vLYX,1)),1)==1),I(triu(ones(size(vLYX,1)),1)==1),'type','Pearson') 
% small POSITIVE correlation with I ~ 17%  
M1=-(vLYX-diag(diag(nan(10))));  corr(M1(triu(ones(size(vLYX,1)),1)==1),I(triu(ones(size(vLYX,1)),1)==1),'type','Pearson')
corr(M(triu(ones(size(vLYX,1)),1)==1),M1(triu(ones(size(vLYX,1)),1)==1),'type','Pearson')
% small NEGATIVE correlation with I ~ 10%  
M2=-(uLYX-diag(diag(nan(10))))'; corr(M2(triu(ones(size(vLYX,1)),1)==1),I(triu(ones(size(vLYX,1)),1)==1),'type','Pearson') 
% small POSITIVE correlation with I ~ 30%  
M3=-(uLYX-diag(diag(nan(10))));  corr(M3(triu(ones(size(vLYX,1)),1)==1),I(triu(ones(size(vLYX,1)),1)==1),'type','Pearson')
% small NEGATIVE between the two dependencies ~ 20%
corr(M2(triu(ones(size(vLYX,1)),1)==1),M3(triu(ones(size(vLYX,1)),1)==1),'type','Pearson')
% large POSITIVE correlation between the two measures ~ 80%
corr(M(triu(ones(size(vLYX,1)),1)==1),M2(triu(ones(size(vLYX,1)),1)==1),'type','Pearson')
% large POSITIVE correlation between the two measures ~ 40%
corr(M1(triu(ones(size(vLYX,1)),1)==1),M3(triu(ones(size(vLYX,1)),1)==1),'type','Pearson')
% reading: 
% vLYX(r,c)  Impact of c -> on r :: columns -> rows
% BUT THE FIGURE IS TRANSPOSED!!!!! therefore y-axis -> x-axis in the figure 
figure
M=-(vLYX-diag(diag(nan(10))))';
im=imagesc(M)
set(im,'AlphaData',~isnan(M))
colormap('hot')
colorbar
hold on
[i,j]=find(isnan(M));
plot(i,j,'xb','markersize',30)
colormap('hot')
colorbar
% h = text(5,2,'x impacting y','fontname','times','fontsize',32)
% set(h,'Rotation',-45);
% h = text(1.5,4.5,'x impacted by y','fontname','times','fontsize',32)
% set(h,'Rotation',-45);
xticklabels(cellstr(SectorNames))
xtickangle(60)
yticklabels(cellstr(SectorNames))
title('Imapct map')
set(gca,'fontsize',14,'fontname','times')
ylabel('stress','interpreter','latex','fontsize',20)
xlabel('response','interpreter','latex','fontsize',20)

%%
figure
plot([0 2],[0 2],'m')
hold on
plot(-vImpacted,-vImpacting,'oc','markersize',18) % impact of the other on a given sector
%plot(-Impacted,-Impacting,'.b') % impact of the other on a given sector
text(-vImpacted,-vImpacting,SectorNames,'fontname','times','fontsize',14)
xlim([0.02 .25])
ylim([0.04 .105])
set(gca,'fontsize',24)
xlabel('Average losses on the sector (VaR$_{0.95}$ stress)','interpreter','latex','fontsize',20)
ylabel('Average losses from the sector (Var$_{0.95}$ stress)','interpreter','latex','fontsize',20)

%%
figure
im=imagesc((I+I'-diag(diag(nan(10)))))
set(im,'AlphaData',~isnan((I+I'-diag(diag(nan(10))))))
colormap('hot')
colorbar
hold on
[i,j]=find(isnan((I+I'-diag(diag(nan(10))))));
plot(i,j,'xb','markersize',30)
xticklabels(cellstr(SectorNames))
xtickangle(60)
yticklabels(cellstr(SectorNames))
%
title('Mutual Information')
set(gca,'fontsize',14,'fontname','times')


%% angle rotation  
% large negative correlation with I ~ 70%
M=mod180(Angle-diag(diag(Angle)))'; corr(M(triu(ones(size(vLYX,1)),1)==1),I(triu(ones(size(vLYX,1)),1)==1),'type','Pearson') 
M=mod180(Angle-diag(diag(Angle))); corr(M(triu(ones(size(vLYX,1)),1)==1),I(triu(ones(size(vLYX,1)),1)==1),'type','Pearson')
M=mod180(Angle-diag(diag(Angle)))';M1=-(vLYX-diag(diag(nan(10))))';corr(M(triu(ones(size(vLYX,1)),1)==1),M1(triu(ones(size(vLYX,1)),1)==1)),corr(M(tril(ones(size(vLYX,1)),-1)==1),M1(tril(ones(size(vLYX,1)),-1)==1))
% reading: 
% Impact of c -> on r :: columns -> rows
% BUT THE FIGURE IS TRANSPOSED!!!!!
figure
M = mod180(Angle-diag(diag(Angle)))';
im=imagesc(M)
set(im,'AlphaData',~isnan(M))
colormap('hot')
colorbar
hold on
[i,j]=find(isnan(M));
plot(i,j,'xb','markersize',30)
xticklabels(cellstr(SectorNames))
xtickangle(60)
yticklabels(cellstr(SectorNames))
title('Principal axis rotation angle')
set(gca,'fontsize',14,'fontname','times')
ylabel('stress','interpreter','latex','fontsize',20)
xlabel('response','interpreter','latex','fontsize',20)

%% eigenvalue change , 
% large negative correlation with I ~ 60%
M=DXY;  corr(M(triu(ones(size(vLYX,1)),1)==1),I(triu(ones(size(vLYX,1)),1)==1),'type','Pearson') 
M=DXY'; corr(M(triu(ones(size(vLYX,1)),1)==1),I(triu(ones(size(vLYX,1)),1)==1),'type','Pearson')
% reading: 
% Impact of c -> on r :: columns -> rows
% BUT THE FIGURE IS TRANSPOSED!!!!!
figure
M = DXY';
im=imagesc(M)
set(im,'AlphaData',~isnan(M))
colormap('hot')
colorbar
hold on
[i,j]=find(isnan(M));
plot(i,j,'xb','markersize',30)
xticklabels(cellstr(SectorNames))
xtickangle(60)
yticklabels(cellstr(SectorNames))
title('Max eigenvalue change')
set(gca,'fontsize',14,'fontname','times')
ylabel('stress','interpreter','latex','fontsize',20)
xlabel('response','interpreter','latex','fontsize',20)


%% Total Variance
M=TvXY'; corr(M(triu(ones(size(vLYX,1)),1)==1),I(triu(ones(size(vLYX,1)),1)==1),'type','Pearson') 
M=TvXY; corr(M(triu(ones(size(vLYX,1)),1)==1),I(triu(ones(size(vLYX,1)),1)==1),'type','Pearson') 
% reading: 
% Impact of c -> on r :: columns -> rows
% BUT THE FIGURE IS TRANSPOSED!!!!!
figure
M = TvXY';
im=imagesc(M)
set(im,'AlphaData',~isnan(M))
colormap('hot')
colorbar
hold on
[i,j]=find(isnan(M));
plot(i,j,'xb','markersize',30)
xticklabels(cellstr(SectorNames))
xtickangle(60)
yticklabels(cellstr(SectorNames))
title('Total Variance log-change ')
set(gca,'fontsize',14,'fontname','times')

%% Mahlanobis impact factor
corr((dX./pX-1)',sum(I+I')') % string correlation with cumulate I ~ 50%
figure
plot(dX./pX,'o-')%(dX-pX)./pX,'o-')
xticklabels(cellstr(SectorNames))
xtickangle(60)
title('Mahlanobis impact factor')
set(gca,'fontsize',14,'fontname','times')
xlabel('conditioning variables','interpreter','latex','fontsize',20)
ylabel('$d_{\mathbf x_q}/p_{\mathbf X}$','interpreter','latex','fontsize',20)
