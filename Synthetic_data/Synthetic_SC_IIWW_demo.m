%% IIWLasso demo
%% Zhao 
clear
list={'datasets','functions'};
for ii=1:length(list)
    addpath(genpath(list{ii}),'-begin');
end
%% initialization
% Experimental setup 
N = 100;
fig=1;figno=fig:1:fig+15;
caselist={'Lasso';'IILasso';'IIWW0.7'};
caseolist={'Lasso';'IILasso';'IIWW0.7';'Orig'};
Diccase='CohGaussian';%'Gaussian'
minsim=0.2;% Minimum coherence for dictionary atoms
%0 for D0.15; 0.2 for D0.5
p = [1,1,0.7]; q=0.7;
casenumber=size(caselist,1);
eachtrailplot=0;eachlambdaplot=1;
% Data generation paramters
m = 30;           % Size of A
n = 50;
divp=[3]; %number of nonzeros
numtrail=size(divp,2);
SNRdB = 20; % SNR of generated data (Btrain), 0 for no noise
maxiter =100;
lambdaratio=1e-5;%2e-6;
lambda=[0.013894954943731, 0.013894954943731,0.007196856730012];
gamma=[0, 0.007196856730012, 1.000000000000000e-03];
Record.gamma=gamma;
convergetol = 1e-5; rzero=1e-3;
%fprintf(caselist);
%bestZerr=inf*ones(numtrail,casenumber);
%% Generate Dictionary
switch Diccase
    case 'Gaussian'
        % normal Gaussian Dic
        D = randn(m,n);
        D = D*diag(1./sqrt(sum(D.*D)));
        min sim Dic
    case 'CohGaussian'
        % min sim Dic
        D= randn(m,n); count=2; repeat=1;
        D = D*diag(1./sqrt(sum(D.*D)));SD=abs(D'*D);
        while min(min(SD))<minsim
            if min(min(SD(1:count,1:count)))<minsim
                D(:,count) = D(:,1)++(rand+1)/6*randn(m,1);%(1-minsim)*(rand+1)/2*randn(m,1);
                D(:,count) = D(:,count)/sqrt(D(:,count)'*D(:,count));
                SD=abs(D'*D);%-eye(n);
            end
            if count==n
                count=1;repeat=repeat+1;
            end
            count=count+1;
        end
end
%similarity of Dic
DTD=abs(D'*D);
Diag1=diag(ones(size(DTD,1),1));
DTD=DTD.*(1-Diag1)+Diag1;
Sim=DTD(logical(triu(ones(size(DTD)))));
Sim(Sim==1)=[];
cohermean=sum(Sim)/size(Sim,1);
barinternal=0.01;barSimD=zeros(1/barinternal,2);
for i=1:(1/barinternal)
    barSimD(i,1)=(i-1)*barinternal;
    sim1=Sim(Sim<(i*barinternal));
    sim1=sim1(sim1>=((i-1)*barinternal));
    barSimD(i,2)=size(sim1,1);
end
figure(figno(12));clf
%subplot(2,1,2)
bar(barSimD(:,1),barSimD(:,2),1);
title(sprintf('Dic column similarity distribution, mean: %.2f',cohermean))
minx=1;
%% trail begin
tic
% Generate data and coef
[X, Zorig, Xorig] = gererateDataColNnD(D, N, divp, minx, SNRdB);
%[X, Zorig, Xorig] = gererateDataBnD(D, N, divp(trail), minx, SNRdB);
posorig=abs(Zorig(:))>0;posorignon0=double(posorig);
posorignon0(posorignon0==0)=2;
% R setting
r=abs(D'*D)/norm(D,'fro').^2*n;
R2=r.*(ones(n,n)-eye(n));
R3=r./(1-r+eye(n)).*(ones(n,n)-eye(n));
count=1;
counta=1;
best.Zerr(count,:)=inf*ones(1,casenumber);
%% Run the experiments
% initialize alpha
eigv=eig(D'*D);
alpha=max(eigv(:))*1.01;%max(max(eigv(:)),0.5);
% initialize H W t
%W=D'/alpha;
H=eye(n)-D'*D/alpha;
B=D'*X/alpha;
tic
nancount=0;
Z=zeros(n,N,casenumber);Zk=zeros(n,N,casenumber);
t=lambda/alpha;%.*ones(n,N,casenumber)
cin=B.*ones(n,N,casenumber);
%% -- Algorithm begins here --
for iter = 1:maxiter
    % update Z proximal splitting
    for j=1:casenumber
        switch caselist{j}
            case 'ISTA'
                [Z(:,:,j),~,~]=softthl1(cin(:,:,j),t(:,:,j));
                cin(:,:,j)=B+H*(Z(:,:,j));
            case 'IILasso'
                for k=1:n
                    dk=D(:,k); zk=Z(k,:,j);
                    Xk=X-D*Z(:,:,j)+dk*zk;
                    Znj=Z(:,:,j);Znj(k,:)=0;
                    Rknk=R2(k,:);Rknk(:,k)=0;
                    tII=gamma(j)+lambda(j).*( Rknk*abs(Znj));
                    [Z(k,:,j)]=softthl1(dk'*Xk,tII);
                end
            case 'Lasso'
                for k=1:n
                    dk=D(:,k); zk=Z(k,:,j);
                    Xk=X-D*Z(:,:,j)+dk*zk;
                    [Z(k,:,j)]=softthl1(dk'*Xk,lambda(j));
                end
                
            case 'WLasso'
                for k=1:n
                    dk=D(:,k); zk=Z(k,:,j);
                    Xk=X-D*Z(:,:,j)+dk*zk;
                    if iter>=1
                        Z1=Z(k,:,j); Z1(Z1==0)=rzero;
                        tL=abs(Z1).^(p(j)-1);
                    else
                        tL=1;
                    end
                    [Z(k,:,j)]=softthl1(dk'*Xk,lambda(j).*tL);
                end
            case 'II-ISTA'
                t(:,:,j)=gamma(j)+lambda(j).*( R2*abs(Z(:,:,j)));
                [Z(:,:,j),~,~]=softthl1(cin(:,:,j),t(:,:,j));
                cin(:,:,j)=B+H*(Z(:,:,j));
            case 'IIWLasso'
                for k=1:n
                    dk=D(:,k); zk=Z(k,:,j);
                    Xk=X-D*Z(:,:,j)+dk*zk;
                    Rknk=R2(k,:);Rknk(:,k)=0;
                    
                    Znj=Z(:,:,j);Znj(k,:)=0;
                    if iter>=1
                        Z1=Z(k,:,j); Z1(Z1==0)=rzero;
                        tL=abs(Z1).^(p(j)-1);
                    else
                        tL=1;
                    end
                    tIW=gamma(j)*tL+lambda(j)*( Rknk*abs(Znj));
                    [Z(k,:,j)]=softthl1(dk'*Xk,tIW);
                end
            case {'IIWR0.7','IIWR0.9'}
                for k=1:n
                    dk=D(:,k); zk=Z(k,:,j);
                    Xk=X-D*Z(:,:,j)+dk*zk;
                    Znj=Z(:,:,j);Znj(k,:)=0;
                    Rknk=R2(k,:);Rknk(:,k)=0;
                    if iter>=1
                        Z1=Z(k,:,j); Z1(Z1==0)=rzero;
                        tL=abs(Z1).^(p(j)-1);
                    else
                        tL=1;
                    end
                    tIW=gamma(j)+lambda(j).*( Rknk*abs(Znj).^(p(j)).*tL);
                    [Z(k,:,j)]=softthl1(dk'*Xk,tIW);
                end
            case {'IIWW0.7'}
                for k=1:n
                    dk=D(:,k); zk=Z(k,:,j);
                    Xk=X-D*Z(:,:,j)+dk*zk;
                    Znj=Z(:,:,j);Znj(k,:)=0;
                    Rknk=R2(k,:);Rknk(:,k)=0;
                    if iter>=1
                        Z1=Z(k,:,j); Z1(Z1==0)=rzero;
                        tL=abs(Z1).^(p(j)-1);
                    else
                        tL=1;
                    end
                    tIW=tL.*(gamma(j)+ lambda(j).*Rknk*abs(Znj).^(q));
                    [Z(k,:,j)]=softthl1(dk'*Xk,tIW);
                end
                
            case 'II2/3PO'
                t(:,:,j)=gamma(j)+lambda(j).*( R2*abs(Z(:,:,j)).^(p(j)));
                [Z(:,:,j),~,~]=PO23(cin(:,:,j),t(:,:,j));
                cin(:,:,j)=B+H*(Z(:,:,j));
                
        end
        % Calculate err
        Zerr(j,iter)=mean(diag((Z(:,:,j)-Zorig)'*(Z(:,:,j)-Zorig))./diag(Zorig'*Zorig));
        resid = D*Z(:,:,j) - X;
        datafitting(j,iter)=0.5*norm(resid,'fro')^2;
        regularization(j,iter)=sum(sum(abs(Z(:,:,j)).^p(j)));
        % Calculate numerosity/diversity
        diversity_hoyer = numerosity_hoyer(Z(:,:,j));
        avgdiversity_hoyer(j,iter) = mean(diversity_hoyer);
        % Calculate Support error
        Z1=Z(:,:,j);
        pos=abs(Z1(:))>0;posnon0=double(pos);
        posnon0(posnon0==0)=-2;
        accn0=sum(posnon0==posorignon0);posnum=max(sum(pos),sum(posorig));
        supporterr(j,iter)=(posnum-accn0)/posnum;
        accuracy(j,iter)=(n*N-sum(pos==posorig))/n*N;
    end
    if iter>2
        for j=1:casenumber
            Zchange(j,iter)=(norm(Z(:,:,j)-Zk(:,:,j),'fro'))/norm(Zk(:,:,j),'fro');
            if isnan(Zchange(j,iter))==1
                nancount=nancount+1;
                nancase=j;
            end
        end
        if nancount>=10
            break
        end
        if max(Zchange(:,iter))<convergetol%*N
            break
        end
    end
    Zk=Z;
end
%calculate
totaliter=iter;
diversity_hoyer0 = numerosity_hoyer(Zorig);
avgdiversity_hoyer0= mean(diversity_hoyer0);
avgdiversity_hoyer0v=avgdiversity_hoyer0*ones(1, totaliter);

if eachlambdaplot==1
    
    figure(figno(1)); clf
    subplot(2, 1, 1);
    hold on
    title(sprintf([caseolist{1} ': %.2e; ' caseolist{2} ': %.2e; ' caseolist{3} ': %.2e'],Zerr(1,totaliter),Zerr(2,totaliter),Zerr(3,totaliter)));
    for j=1:casenumber
        plot(Zerr(j,1:totaliter));
    end
    axis([0 totaliter 0 1]);
    xlabel('Iteration');
    ylabel('Relative norm error');
    legend(caselist);
    subplot(2, 1, 2);
    hold on
    axis([0 totaliter 0 1]);
    title(sprintf('Original Hoyer = %.3g',  avgdiversity_hoyer0) );
    xlabel('Iteration');
    ylabel('Average Hoyer Sparsity');
    xlabel('Iteration');
    for j=1:casenumber
        plot(avgdiversity_hoyer(j,1:totaliter));
    end
    plot(avgdiversity_hoyer0v,'Color',[1 0 0],'LineWidth',2);
    legend(caseolist);
    
    resid=D*Zorig-X;
    regularization0=sum(sum(abs(Zorig).^p(1)));
    datafitting0=0.5*norm(resid,'fro')^2;
    figure(figno(2)); clf;
    subplot(2, 1, 1);
    hold on
    plot(datafitting0(ones(1,totaliter)),'Color',[1 0 0],'LineWidth',1);
    title(sprintf('Iteration %d |D|_F=%5.3f, A = (%d x %d), p = %5.3f', ...
        iter, norm(D,'fro'), m, n , p(1) ));
    ylabel('datafitting');
    for j=1:casenumber
        plot(datafitting(j,1:totaliter));
    end
    legend(caseolist);
    % xlabel('Iteration');
    subplot(2,1,2)
    hold on
    ylabel('regularization');
    for j=1:casenumber
        plot(regularization(j,1:totaliter));
    end
    plot(regularization0(ones(1,totaliter)),'Color',[1 0 0],'LineWidth',1);
    legend(caseolist);
    
    figure(figno(3));clf
    subplot(2, 1, 1);
    hold on
    for j=1:casenumber
        plot(supporterr(j,1:totaliter));
    end
    legend(caselist);
    ylabel('Support Error');
    subplot(2, 1, 2);
    hold on
    for j=1:casenumber
        plot(accuracy(j,1:totaliter));
    end
    legend(caselist);
    ylabel('Accuracy');
    
    figure(figno(4));clf
    pos=find(abs(Zorig)>0);
    hold on
    for j=1:casenumber
        Zp=Z(:,:,j);
        pos=find(abs(Zp)>0);
        stem(pos,Zp(pos));
    end
    stem(pos,Zorig(pos),'red','LineWidth',1);
    legend(caseolist);
    axis([0 n -1 1]);
    hold off
end
