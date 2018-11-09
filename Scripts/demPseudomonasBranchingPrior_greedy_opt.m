function [Fin] = demPseudomonasBranchingPrior_greedy_opt(batchi);

%Fit branching process to the Arabidopsis/pseudomonas dataset (Lewis et al., 2015).
%Split the analysis into batches of 1000 for running on the cluster. 

addpath(genpath('../'))
 
D1 = importdata('./Pa13-Combo2.txt');

try %Try to resume analysis 
    load(['./results/pseudomonas/PseudomonasResults_Feb_2018_5_' num2str(batchi) '_greedy.mat'])
    startind = length(Ps)+1;
    endind   = (batchi)*1000;
catch
    startind = (batchi-1)*1000 + 1;
    endind   = (batchi)*1000;
end

endind = min(endind,size(D1.data,1));

%Measurement times. Permute these for different branching structures. 
tt = [0,2,3,4,6,7,8,10,11,12,14,16,17.5];
X1 = [repmat(tt,1,12); ones(1,52),2*ones(1,52),3*ones(1,52)]';
X2 = [repmat(tt,1,12); ones(1,52),ones(1,52),2*ones(1,52)]';
X3 = [repmat(tt,1,12); ones(1,52),2*ones(1,52),2*ones(1,52)]';

%Prediction time points.
Xstar1 = [repmat(linspace(0,17.5,50),1,3);ones(1,50),2*ones(1,50),3*ones(1,50)]';
Xstar2 = [repmat(linspace(0,17.5,50),1,3);ones(1,50),1*ones(1,50),2*ones(1,50)]';
Xstar3 = [repmat(linspace(0,17.5,50),1,3);ones(1,50),2*ones(1,50),2*ones(1,50)]';
tic
for i = startind:endind

pcp1     = {@priorGamma,2,2};    %Mean 4, std 8 
pcp1p2   = {@priorGamma,4,2};    %Mean 8, std 16 
pcp2     = {@priorGauss,0,0.5};  %Quick transitions        
pctheta1 = {@priorGauss,3,1};    %Length-scale              
pctheta2 = {@priorGauss,1,1};    %Variance             
pcmean   = {@priorGauss,0,1};    %Zero, unit-variance
pclik    = {@priorGauss,0.1,0.5};
    
Y1 = D1.data(i,:)'; %Mock hrp DC
Y2 = D1.data(i,[1:52,2*52+1:3*52,52+1:2*52])'; %Mock DC hrpA

%Initialise some parameters
l1 = log(3); l2 = log(3); lg = log(3); v1 = log(3); v2 = log(3); vg = log(std(Y1));

%Joint GP (not DE). 
hyp.cov = [l1;v1]; hyp.mean = mean(Y1(1:52,1)); hyp.lik = 2;
prior.mean = {[]};  prior.cov  = {[];[]}; prior.lik = {pclik};
im = {@infPrior,@infExact,prior};                % inference method
par1a = {'meanConst','covSEiso','likGauss',X1(1:52,1),Y1(1:52)};
par1b = {'meanConst','covSEiso','likGauss',X1(1:52,1),Y1(1:52),Xstar1(:,1)};
hyp_pN7 = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         % optimise
[L7 dL7] = feval(@gp,hyp_pN7, im, par1a{:});         % optimise
[ymu7 ys27 fmu7 fs27   ]= feval(@gp, hyp_pN7, im, par1b{:});

%Now try merging ...    
hyp.cov = [4;1;hyp_pN7.cov;hyp_pN7.cov];hyp.mean = mean(Y1(1:2*52,1));hyp.lik = hyp_pN7.lik;
prior.mean = {[]};  prior.cov  = {pcp1;pcp2;[];[];[];[]}; prior.lik = {[]};
im = {@infPrior,@infExact,prior};                % inference method
par1a = {'meanConst','covBranchingProcess_2A','likGauss',X1(1:2*52,:),Y1(1:2*52)};
par1b = {'meanConst','covBranchingProcess_2A','likGauss',X1(1:2*52,:),Y1(1:2*52),Xstar2};
par1c = {'meanConst','covBranchingProcess_2A','likGauss',X1(1:2*52,:),Y1(1:2*52),hyp,sq_dist(X1(1:2*52,1)')};
hyp_pN3a = feval(@minimize, hyp, @gpfunc_2A, -20000, im, par1c{:});
[L3a dL3a] = feval(@gpfunc_2A,hyp_pN3a, im, par1c{:});         % optimise
[ymu3a ys23a fmu3a fs23a   ]= feval(@gp,hyp_pN3a, im, par1b{:});

hyp.cov = [4;1;hyp_pN7.cov;hyp_pN7.cov];hyp.mean = mean(Y2(1:2*52,1));hyp.lik = hyp_pN7.lik;
prior.mean = {[]};  prior.cov  = {pcp1;pcp2;[];[];[];[]}; prior.lik = {[]};
im = {@infPrior,@infExact,prior};                % inference method
par1a = {'meanConst','covBranchingProcess_2A','likGauss',X1(1:2*52,:),Y2(1:2*52)};
par1b = {'meanConst','covBranchingProcess_2A','likGauss',X1(1:2*52,:),Y2(1:2*52),Xstar2};
par1c = {'meanConst','covBranchingProcess_2A','likGauss',X1(1:2*52,:),Y2(1:2*52),hyp,sq_dist(X1(1:2*52,1)')};
hyp_pN3b = feval(@minimize, hyp, @gpfunc_2A, -20000, im, par1c{:});
[L3b dL3b] = feval(@gpfunc_2A,hyp_pN3b, im, par1c{:});         % optimise
[ymu3b ys23b fmu3b fs23b   ]= feval(@gp,hyp_pN3b, im, par1b{:});
    
%Joint GP (not DE). 
hyp = hyp_pN7;
prior.mean = {[]};  prior.cov  = {[];[]}; prior.lik = {pclik};
im = {@infPrior,@infExact,prior};                % inference method
par1a = {'meanConst','covSEiso','likGauss',X1(1:2*52,1),Y1(1:2*52)};
par1b = {'meanConst','covSEiso','likGauss',X1(1:2*52,1),Y1(1:2*52),Xstar1(:,1)};
hyp_pN3c = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         % optimise
[L3c dL3c] = feval(@gp,hyp_pN3c, im, par1a{:});         % optimise
[ymu3c ys23c fmu3c fs23c   ]= feval(@gp,hyp_pN3c, im, par1b{:});

hyp = hyp_pN7;
prior.mean = {[]};  prior.cov  = {[];[]}; prior.lik = {pclik};
im = {@infPrior,@infExact,prior};                % inference method
par1a = {'meanConst','covSEiso','likGauss',X1(1:2*52,1),Y2(1:2*52)};
par1b = {'meanConst','covSEiso','likGauss',X1(1:2*52,1),Y2(1:2*52),Xstar1(:,1)};
hyp_pN3d = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         % optimise
[L3d dL3d] = feval(@gp,hyp_pN3d, im, par1a{:});         % optimise
[ymu3d ys23d fmu3d fs23d   ]= feval(@gp,hyp_pN3d, im, par1b{:});
   
%Use ML to distuinguish 
Liks = [L3a,L3b,L3c,L3d];

Liks2 = - 2*Liks + [8,8,4,4]*log(size(X1,1));

ind = find(Liks2==max(Liks2));

if ind==4 %
 
hyp.cov = [4;1;hyp_pN3d.cov;hyp_pN3d.cov];hyp.mean = hyp_pN3d.mean;hyp.lik = hyp_pN3d.lik;
prior.mean = {[]};  prior.cov  = {pcp1;pcp2;[];[];[];[]}; prior.lik = {[]};
im = {@infPrior,@infExact,prior};                % inference method
par1a = {'meanConst','covBranchingProcess_2A','likGauss',X2,Y2};
par1b = {'meanConst','covBranchingProcess_2A','likGauss',X2,Y2,Xstar2};
par1c = {'meanConst','covBranchingProcess_2A','likGauss',X2,Y2,hyp,sq_dist(X2')};
hyp_pN4a = feval(@minimize, hyp, @gpfunc_2A, -20000, im, par1c{:});
[L4a dL4a] = feval(@gpfunc_2A,hyp_pN4a, im, par1c{:});  
[ymu4a ys24a fmu4a fs24a   ]= feval(@gp,hyp_pN4a, im, par1b{:});        

%Joint    
hyp = hyp_pN3d;    
prior.mean = {[]};  prior.cov  = {[];[]}; prior.lik = {pclik};
im = {@infPrior,@infExact,prior};              
par1a = {'meanConst','covSEiso','likGauss',X1(:,1),Y1};
par1b = {'meanConst','covSEiso','likGauss',X1(:,1),Y1,Xstar1(:,1)};
hyp_pN4b = feval(@minimize, hyp, @gp, -20000, im, par1a{:});    
[L4b dL4b] = feval(@gp,hyp_pN4b, im, par1a{:});
[ymu4b ys24b fmu4b fs24b   ]= feval(@gp,hyp_pN4b, im, par1b{:});
    
Lo = [L4a,L4b];

Lo2 = - 2*Liks + [8,4]*log(size(X1,1));


elseif ind==3
    
hyp.cov = [4;1;hyp_pN3c.cov;hyp_pN3c.cov];hyp.mean = hyp_pN3c.mean;hyp.lik = hyp_pN3c.lik;
prior.mean = {[]};  prior.cov  = {pcp1;pcp2;[];[];[];[]}; prior.lik = {[]};
im = {@infPrior,@infExact,prior};                % inference method
par1a = {'meanConst','covBranchingProcess_2A','likGauss',X2,Y1};
par1b = {'meanConst','covBranchingProcess_2A','likGauss',X2,Y1,Xstar2};
par1c = {'meanConst','covBranchingProcess_2A','likGauss',X2,Y1,hyp,sq_dist(X2')};
hyp_pN5a = feval(@minimize, hyp, @gpfunc_2A, -20000, im, par1c{:});
[L5a dL5a] = feval(@gpfunc_2A,hyp_pN5a, im, par1c{:});  
[ymu5a ys25a fmu5a fs25a   ]= feval(@gp,hyp_pN5a, im, par1b{:});        

%Joint    
hyp = hyp_pN3c;    
prior.mean = {[]};  prior.cov  = {[];[]}; prior.lik = {pclik};
im = {@infPrior,@infExact,prior};              
par1a = {'meanConst','covSEiso','likGauss',X1,Y1};
par1b = {'meanConst','covSEiso','likGauss',X1,Y1,Xstar1(:,1)};
hyp_pN5b = feval(@minimize, hyp, @gp, -20000, im, par1a{:});    
[L5b dL5b] = feval(@gp,hyp_pN5b, im, par1a{:});
[ymu5b ys25b fmu5b fs25b   ]= feval(@gp,hyp_pN5b, im, par1b{:});

Lo = [L5a,L5b];
Lo2 = - 2*Liks + [8,4]*log(size(X1,1));

elseif ind==1
        
hyp.cov = [hyp_pN3a.cov];hyp.mean = hyp_pN3a.mean;hyp.lik = hyp_pN3a.lik;
prior.mean = {[]};  prior.cov  = {pcp1;pcp2;[];[];[];[]}; prior.lik = {[]};
im = {@infPrior,@infExact,prior};                % inference method
par1a = {'meanConst','covBranchingProcess_2A','likGauss',X2,Y2};
par1b = {'meanConst','covBranchingProcess_2A','likGauss',X2,Y2,Xstar2};
par1c = {'meanConst','covBranchingProcess_2A','likGauss',X2,Y2,hyp,sq_dist(X2')};
hyp_pN6a = feval(@minimize, hyp, @gpfunc_2A, -20000, im, par1c{:});
[L6a dL6a] = feval(@gpfunc_2A,hyp_pN6a, im, par1c{:});  
[ymu6a ys26a fmu6a fs26a   ]= feval(@gp,hyp_pN6a, im, par1b{:});

hyp.cov = hyp_pN3a.cov;hyp.mean = hyp_pN3a.mean;hyp.lik = hyp_pN3a.lik;
prior.mean = {[]};  prior.cov  = {pcp1;pcp2;[];[];[];[]}; prior.lik = {[]};
im = {@infPrior,@infExact,prior};                % inference method
par1a = {'meanConst','covBranchingProcess_2A','likGauss',X3,Y2};
par1b = {'meanConst','covBranchingProcess_2A','likGauss',X3,Y2,Xstar2};
par1c = {'meanConst','covBranchingProcess_2A','likGauss',X3,Y2,hyp,sq_dist(X3')};
hyp_pN6b = feval(@minimize, hyp, @gpfunc_2A, -20000, im, par1c{:});
[L6b dL6b] = feval(@gpfunc_2A,hyp_pN6b, im, par1c{:});  
[ymu6b ys26b fmu6b fs26b   ]= feval(@gp,hyp_pN6b, im, par1b{:});
        
%3 branch structure. Mock->hrp->DC %Rearranging the hyperparameters here ....
hyp.cov  = [4;1;hyp_pN3a.cov(1:2);hyp_pN3a.cov(3:4);hyp_pN3a.cov(3:6)]; hyp.mean = hyp_pN3a.mean; hyp.lik  = hyp_pN3a.lik;
prior.mean = {[]};  prior.cov  = {pcp1;pcp2;pcp1p2;pcp2;[];[];[];[];[];[]}; prior.lik = {[]};
im = {@infPrior,@infExact,prior};                % inference method
par1a = {'meanConst','covBranchingProcess_3A','likGauss',X1,Y1};
par1b = {'meanConst','covBranchingProcess_3A','likGauss',X1,Y1,Xstar1};
par1c = {'meanConst','covBranchingProcess_3A','likGauss',X1,Y1,hyp,sq_dist(X1(:,1)')};
hyp_pN6c = feval(@minimize, hyp, @gpfunc_3A, -20000, im, par1c{:});
[L6c dL6c] = feval(@gpfunc_3A,hyp_pN6c, im, par1c{:});         % optimise
[ymu6c ys26c fmu6c fs26c   ]= feval(@gp,hyp_pN6c, im, par1b{:});
    
%3 branch structure. Mock->hrp->DC %Rearranging the hyperparameters here ....
hyp.cov  = [4;1;hyp_pN3a.cov(1:2);hyp_pN3a.cov(3:4);hyp_pN3a.cov(3:6)]; hyp.mean = hyp_pN3a.mean; hyp.lik  = hyp_pN3a.lik;
prior.mean = {[]};  prior.cov  = {pcp1;pcp2;pcp1p2;pcp2;[];[];[];[];[];[]}; prior.lik = {[]};
im = {@infPrior,@infExact,prior};                % inference method
par1a = {'meanConst','covBranchingProcess_3A','likGauss',X1,Y2};
par1b = {'meanConst','covBranchingProcess_3A','likGauss',X1,Y2,Xstar1};
par1c = {'meanConst','covBranchingProcess_3A','likGauss',X1,Y2,hyp,sq_dist(X1(:,1)')};
hyp_pN6d = feval(@minimize, hyp, @gpfunc_3A, -20000, im, par1c{:});
[L6d dL6d] = feval(@gpfunc_3A,hyp_pN6d, im, par1c{:});         % optimise
[ymu6d ys26d fmu6d fs26d   ]= feval(@gp,hyp_pN6d, im, par1b{:});

Lo = [L6a,L6b,L6c,L6d];
Lo2 = - 2*Liks + [8,8,12,12]*log(size(X1,1));

else
    
hyp.cov = [hyp_pN3b.cov];hyp.mean = hyp_pN3b.mean;hyp.lik = hyp_pN3b.lik;
prior.mean = {[]};  prior.cov  = {pcp1;pcp2;[];[];[];[]}; prior.lik = {[]};
im = {@infPrior,@infExact,prior};                % inference method
par1a = {'meanConst','covBranchingProcess_2A','likGauss',X2,Y2};
par1b = {'meanConst','covBranchingProcess_2A','likGauss',X2,Y2,Xstar2};
par1c = {'meanConst','covBranchingProcess_2A','likGauss',X2,Y2,hyp,sq_dist(X2')};
hyp_pN8a = feval(@minimize, hyp, @gpfunc_2A, -20000, im, par1c{:});
[L8a dL8a] = feval(@gpfunc_2A,hyp_pN8a, im, par1c{:});  
[ymu8a ys28a fmu8a fs28a   ]= feval(@gp,hyp_pN8a, im, par1b{:});

hyp.cov = hyp_pN3b.cov;hyp.mean = hyp_pN3b.mean;hyp.lik = hyp_pN3b.lik;
prior.mean = {[]};  prior.cov  = {pcp1;pcp2;[];[];[];[]}; prior.lik = {[]};
im = {@infPrior,@infExact,prior};                % inference method
par1a = {'meanConst','covBranchingProcess_2A','likGauss',X3,Y2};
par1b = {'meanConst','covBranchingProcess_2A','likGauss',X3,Y2,Xstar2};
par1c = {'meanConst','covBranchingProcess_2A','likGauss',X3,Y2,hyp,sq_dist(X3')};
hyp_pN8b = feval(@minimize, hyp, @gpfunc_2A, -20000, im, par1c{:});
[L8b dL8b] = feval(@gpfunc_2A,hyp_pN8b, im, par1c{:});  
[ymu8b ys28b fmu8b fs28b   ]= feval(@gp,hyp_pN8b, im, par1b{:});
        
%3 branch structure. Mock->hrp->DC %Rearranging the hyperparameters here ....
hyp.cov  = [4;1;hyp_pN3b.cov(1:2);hyp_pN3b.cov(3:4);hyp_pN3b.cov(3:6)]; hyp.mean = hyp_pN3b.mean; hyp.lik  = hyp_pN3b.lik;
prior.mean = {[]};  prior.cov  = {pcp1;pcp2;pcp1p2;pcp2;[];[];[];[];[];[]}; prior.lik = {[]};
im = {@infPrior,@infExact,prior};                % inference method
par1a = {'meanConst','covBranchingProcess_3A','likGauss',X1,Y1};
par1b = {'meanConst','covBranchingProcess_3A','likGauss',X1,Y1,Xstar1};
par1c = {'meanConst','covBranchingProcess_3A','likGauss',X1,Y1,hyp,sq_dist(X1(:,1)')};
hyp_pN8c = feval(@minimize, hyp, @gpfunc_3A, -20000, im, par1c{:});
[L8c dL8c] = feval(@gpfunc_3A,hyp_pN8c, im, par1c{:});         % optimise
[ymu8c ys28c fmu8c fs28c   ]= feval(@gp,hyp_pN8c, im, par1b{:});
    
%3 branch structure. Mock->hrp->DC %Rearranging the hyperparameters here ....
hyp.cov  = [4;1;hyp_pN3b.cov(1:2);hyp_pN3b.cov(3:4);hyp_pN3b.cov(3:6)]; hyp.mean = hyp_pN3b.mean; hyp.lik  = hyp_pN3b.lik;
prior.mean = {[]};  prior.cov  = {pcp1;pcp2;pcp1p2;pcp2;[];[];[];[];[];[]}; prior.lik = {[]};
im = {@infPrior,@infExact,prior};                % inference method
par1a = {'meanConst','covBranchingProcess_3A','likGauss',X1,Y2};
par1b = {'meanConst','covBranchingProcess_3A','likGauss',X1,Y2,Xstar1};
par1c = {'meanConst','covBranchingProcess_3A','likGauss',X1,Y2,hyp,sq_dist(X1(:,1)')};
hyp_pN8d = feval(@minimize, hyp, @gpfunc_3A, -20000, im, par1c{:});
[L8d dL8d] = feval(@gpfunc_3A,hyp_pN8d, im, par1c{:});         % optimise
[ymu8d ys28d fmu8d fs28d   ]= feval(@gp,hyp_pN8d, im, par1b{:});

Lo = [L8a,L8b,L8c,L8d];
Lo2 = - 2*Liks + [8,8,12,12]*log(size(X1,1));
end

%Store likelihoods
L = -[L7,L3a,L3b,L3c,L3d];
Lfin = Lo;
BICfin = Lo2;    

%L4a,L4b,L5a,L5b,L6a,L6b,L6c,L6d,L7a,L7b,L7c,L7d];
%AIC = 2*[12,12,8,8,4,8] - 2*L;
%BIC = - 2*L + [12,12,8,8,4,8]*log(size(X1,1));

Output.L = L;
Output.Lfin = Lfin;
Output.BICfin = BICfin;
%Output.AIC = AIC;
%Output.BIC = BIC;
Output.H7 = hyp_pN7;
Output.H3a = hyp_pN3a;
Output.H3b = hyp_pN3b;
Output.H3c = hyp_pN3c;
Output.H3d = hyp_pN3d;

if ind==4
Output.H4a = hyp_pN4a;
Output.H4b = hyp_pN4b;
Output.fmu4a = fmu4a;
Output.fmu4b = fmu4b;
Output.fs24a = fs24a;
Output.fs24b = fs24b;
elseif ind==3
Output.H5a = hyp_pN5a;
Output.H5b = hyp_pN5b;
Output.fmu5a = fmu5a;
Output.fmu5b = fmu5b;
Output.fs25a = fs25a;
Output.fs25b = fs25b;
elseif ind==1
Output.H6a = hyp_pN6a;
Output.H6b = hyp_pN6b;
Output.H6c = hyp_pN6c;
Output.H6d = hyp_pN6d;

Output.fmu6a = fmu6a;
Output.fmu6b = fmu6b;
Output.fmu6c = fmu6c;
Output.fmu6d = fmu6d;
Output.fmu6a = fmu6a;
Output.fmu6b = fmu6b;
Output.fmu6c = fmu6c;
Output.fmu6d = fmu6d;
Output.fs26a = fs26a;
Output.fs26b = fs26b;
Output.fs26c = fs26c;
Output.fs26d = fs26d;

else
Output.H8a = hyp_pN8a;
Output.H8b = hyp_pN8b;
Output.H8c = hyp_pN8c;
Output.H8d = hyp_pN8d;

Output.fmu8a = fmu8a;
Output.fmu8b = fmu8b;
Output.fmu8c = fmu8c;
Output.fmu8d = fmu8d;
Output.fmu8a = fmu8a;
Output.fmu8b = fmu8b;
Output.fmu8c = fmu8c;
Output.fmu8d = fmu8d;
Output.fs28a = fs28a;
Output.fs28b = fs28b;
Output.fs28c = fs28c;
Output.fs28d = fs28d;
end


%Output.H4 = hyp_pN4;
%Output.H5 = hyp_pN5;
%Output.H6 = hyp_pN6;
%Output.H7 = hyp_pN7;
%Output.H8 = hyp_pN8;
%Output.fmu1 = fmu1;
%Output.fmu2 = fmu2;
%Output.fmu3 = fmu3;
%Output.fmu4 = fmu4;
%Output.fmu5 = fmu5;%
%Output.fmu6 = fmu6;
%Output.fmu7 = fmu7;
%Output.fmu8 = fmu8;
%Output.fs21 = fs21;
%Output.fs22 = fs22;
%Output.fs23 = fs23;
%Output.fs24 = fs24;
%Output.fs25 = fs25;
%Output.fs26 = fs26;
%Output.fs27 = fs27;
%Output.fs28 = fs28;

Ps{i,1} = Output;
%save(['./results/pseudomonas/PseudomonasResults_Feb_5_' num2str(batchi) '.mat'],'Ps')

     %end

end
toc
save(['./results/pseudomonas/PseudomonasResults_Feb_5_2018_' num2str(batchi) '_greedy.mat'],'Ps')

Fin = 1;
