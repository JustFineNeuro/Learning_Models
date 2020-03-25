clear all
close all
clc;
for rr=1:20
nt=240; %Length of trial

f_fname=@f_twomemorymodel; %evoluation funciton
g_fname=@g_multi_memory_obs; %Observation function (This assumes identity mapping);


dim.n_theta = 4; 
dim.n=3; 

% % Build priors for model inversion
priors.muX0 = [400 ;-400;0]
priors.SigmaX0 = 1e-0*eye(3);
priors.muTheta = 0*ones(4,1);
priors.SigmaTheta = 1e-1*eye(4);
priors.a_alpha = 1e0;
priors.b_alpha = 1e0;
priors.a_sigma = 1e0;
priors.b_sigma = 1e0;
options.inG.statemap = sign(ut);
ut(find(mx<0))=ut(find(mx<0)).*-1;
ut=ut.*-1;
options.in.active=sign(ut);
options.priors=priors;
% Build options and dim structures for model inversion
% options.inF = in;

dim.n_theta = 4;%Number of evolution equation paramters.
dim.n_phi = 0;
dim.n = 3;%number of latentstates in evolution equation
dim.p = 1; disp('not sure what this is, but its for inversion?')
x0 = zeros(dim.n,1);
[posterior,out] = VBA_NLStateSpaceModel(y',ut',f_fname,g_fname,dim,options);
%Evolution equation parameters for simulation
% theta=[0.9;0.3;1/2;1/2];  %Can Extend to MFX by using these as populatio means
theta=[0.9195+(randn(1))/10.5
    0.7249+(randn(1))/10.5
    -0.2418+(randn(1))/10.5
    -0.488+(randn(1))/10.5]

ut=[zeros(40,1);10*ones(50,1);zeros(30,1);-10*ones(120,1)];

alpha   = [5]; %Precision on evolution parameters (alpha = hyperparameter)
sigma   = [3]; %Precision on observation parameters (alpha = hyperparameter)

[y,x,x0,eta,e] = simulateNLSS(nt,f_fname,g_fname,theta,[],ut',10,5,options,x0);

[posterior,out] = VBA_NLStateSpaceModel(y,ut',f_fname,g_fname,dim,options);
post{rr}=posterior;
end
[ehat,v_e,etahat,v_eta] = VBA_getNoise(posterior,out);
[haf,hf,hp] = plotUncertainTimeSeries(ehat,VBA_getVar(v_e));
%%  Hierarchical model fitting -- subjects fits provide a prior.

%Number of subjects 
ns=10;

dim.n_t = 240; % trials (within-subject)
nt=240; %Length of trial

f_fname=@f_twostatemodel; %evoluation funciton
g_fname=@Two_State_Obs; %Observation function (This assumes identity mapping);


dim.n_theta = 4; 
dim.n=2; 

% Build priors for model inversion
priors.muX0 = zeros(2,1);
priors.SigmaX0 = 1e-0*eye(2);
priors.muTheta = 0*ones(4,1);
priors.SigmaTheta = 1e-1*eye(4);
priors.a_alpha = 1e0;
priors.b_alpha = 1e0;
priors.a_sigma = 1e0;
priors.b_sigma = 1e0;

    options{1}.DisplayWin = 0;
    options{1}.verbose = 1;
    options{1}.dim = dim;
    options{1}.binomial = 0;
    options{1}.inG.statemap = [1;1];
    options{1}.priors=priors;




for i=1:ns

theta=[0.9195+(randn(1))/8.5
    0.8249+(randn(1))/8.5
    -0.2418+(randn(1))/8.5
    -0.3888+(randn(1))/8.5]

ut{i}=[zeros(40,1);10*ones(50,1);zeros(30,1);-10*ones(120,1)];
ut=[zeros(40,1);10*ones(50,1);zeros(30,1);-10*ones(120,1)];

alpha   = [9+randn(1)]; %Precision on evolution parameters (alpha = hyperparameter)
sigma   = [5+randn(1)]; %Precision on observation parameters (alpha = hyperparameter)

[y{i},x,x0,eta,e] = simulateNLSS(nt,f_fname,g_fname,theta,[],ut',alpha,sigma,options,x0);
end
clear options
for i=1:ns
    options{i}.DisplayWin = 0;
    options{i}.verbose = 1;
    options{i}.dim = dim;
    options{i}.binomial = 0;
    options{i}.inG.statemap = [1;1];
    options{i}.priors=priors;

end


[p_sub,o_sub,p_group,o_group] = VBA_MFX(y,ut,f_fname,g_fname,dim,options);%,priors_group);


for i=1:ns
    % with MFX-type priors 
    Theta(:,i,1) = p_sub{i}.muTheta;
%     Phi(:,i,1) = p_sub{i}.muPhi;
    X0(:,i,1) = p_sub{i}.muX0;
    % without MFX-type priors
    Theta(:,i,2) = o_group.initVBA.p_sub{i}.muTheta;
%     Phi(:,i,2) = o_group.initVBA.p_sub{i}.muPhi;
    X0(:,i,2) = o_group.initVBA.p_sub{i}.muX0;
end



hf = figure('name','btw-subject variability','color',[1 1 1]);
col = getColors(max([dim.n;dim.n_theta;dim.n_phi]));
tistr = {' (with MFX priors)',' (without MFX priors)'};
for j=1:2
    ha = subplot(2,3,(j-1)*3+1,'parent',hf,'nextplot','add');
    plot(ha,theta',Theta(:,:,j)','.')
    leg = cell(0);
    for i=1:dim.n_theta
        hp = plot(ha,theta(i,:),Theta(i,:,j),'.','color',col(i,:));
        [oa] = addBestLinearPredictor(hp,0);
        if oa.pv<0.05
            leg{i} = ['dim #',num2str(i),': *'];
        else
            leg{i} = ['dim #',num2str(i)];
        end
    end
    axis(ha,'tight')
    title(ha,['evolution params',tistr{j}])
    xlabel(ha,'simulated')
    ylabel(ha,'estimated')
    legend(ha,leg)
    ha = subplot(2,3,(j-1)*3+2,'parent',hf,'nextplot','add');
    plot(ha,phi',Phi(:,:,j)','.')
    leg = cell(0);
    for i=1:dim.n_phi
        hp = plot(ha,phi(i,:),Phi(i,:,j),'.','color',col(i,:));
        [oa] = addBestLinearPredictor(hp,0);
            leg{i} = ['dim #',num2str(i)];
    end
    title(ha,['observation params',tistr{j}])
    xlabel(ha,'simulated')
    ylabel(ha,'estimated')
    legend(ha,leg)
    ha = subplot(2,3,(j-1)*3+3,'parent',hf,'nextplot','add');
    plot(ha,x0',X0(:,:,j)','.')
    leg = cell(0);
    for i=1:dim.n
        hp = plot(ha,x0(i,:),X0(i,:,j),'.','color',col(i,:));
        [oa] = addBestLinearPredictor(hp,0);
        
            leg{i} = ['dim #',num2str(i)];
    end
    title(ha,['initial conditions',tistr{j}])
    xlabel(ha,'simulated')
    ylabel(ha,'estimated')
    legend(ha,leg)
end

%%%%%%%%%%%%%%%%%%%%%%%

%% Lets do it on real data now (Matt's study)
ns=24;

load('ALLDATA.mat')
len=size(Force_Data.Reduced,1);
psycho_idx=1:45;
psycho=[];
ss=[1:24];
clear options

for i=1:ns
    datforce{i}=inpaint_nans(transpose(Force_Data.Reduced.Full_Data{ss(i)}(46:445)));
    temp=inpaint_nans(transpose(Force_Data.Reduced.Full_Data{ss(i)}(46:445)));
    target{i}=(Force_Data.Reduced.Target_Profile{ss(i)}(46:445));
    temptar=(Force_Data.Reduced.Target_Profile{ss(i)}(46:445));
    
    c1{i}=smooth(abs(temp(1:40)-temptar(1:40)'),3);
    c2{i}=smooth(abs(temp(101:40)-temptar(101:40)'),3);
    c3{i}=smooth(abs(temp(201:40)-temptar(201:40)'),3);
    c4{i}=smooth(abs(temp(301:40)-temptar(301:40)'),3);
    
end

 gpA=[1 2 4 5 7 10 12 14 16 18 22 24];
 gpB=find(~ismember(1:24,gpA));

 c1A=smooth(mean(c1(gpA,1:40),1),3); c2A=smooth(mean(c2(gpA,1:40),1),3);c3A=smooth(mean(c3(gpA,1:40),1),3);  c4A=smooth(mean(c4(gpA,1:40),1),3);
 c1B=smooth(mean(c1(gpB,1:40),1),3); c2B=smooth(mean(c2(gpB,1:40),1),3);c3B=smooth(mean(c3(gpB,1:40),1),3);c4B=smooth(mean(c4(gpB,1:40),1),3);

%% Lets use the toolbox to fit exponentials use vba_mfx
c2a=c2(gpA);
c2b=c2(gpB);
f_fname=@f_DoubleExponent; %evoluation funciton
g_fname=@g_ID_Obs; %Observation function (This assumes identity mapping);
% Build options and dim structures for model inversion
% options.inF = in;
dim.n_t=100;
dim.n_theta = 5; 
dim.n=1; 
dim.n_phi = 0;
dim.p = 1; 
x0 = zeros(dim.n,1);
ns=12;
for i=1:ns
    options{i}.isYout=isnan(c1{i});
    options{i}.DisplayWin = 0;
    options{i}.verbose = 1;
    options{i}.dim = dim;
    options{i}.binomial = 0;
    options{i}.backwardLag = 1 ;
end

[p_sub_c1,o_sub_c1,p_group_c1,o_group_c1] = VBA_MFX(c1a,[],f_fname,g_fname,dim,options);%,priors_group);
[p_sub_c2,o_sub_c2,p_group_c2,o_group_c2] = VBA_MFX(c2b,[],f_fname,g_fname,dim,options);%,priors_group);



%% Now will fit the two state learnign error-based model
priors.a_alpha = 1;
priors.b_alpha = 1;
priors.muTheta=([0.9 0.7 -0.1 -0.4])'
priors.SigmaTheta=eye(dim.n_theta);
% priors_group.QPhi = 0.*eye(dim.n_phi);
% priors_group.QTheta = 0.*eye(dim.n_theta);
% priors_group.QX0 = 0.*eye(dim.n);
% priors_group.QPhi(2,2) = 0; % ffx
% priors_group.SigmaPhi = eye(dim.n_phi);
% priors_group.SigmaPhi(2,2) = 0; % fix population mean to 0


dim.n_t =400;
for i=1:ns
    options{i}.isYout=isnan(datforce{i});
    options{i}.DisplayWin = 0;
    options{i}.verbose = 1;
    options{i}.dim = dim;
    options{i}.binomial = 0;
    options{i}.inG.statemap = [1;1];
    options{i}.backwardLag = 3 ;
    options{i}.priors=priors;

end

%% Inversion time across everyone
[p_sub,o_sub,p_group,o_group] = VBA_MFX(datforce,target,f_fname,g_fname,dim,options);%,priors_group);

%% Invert now per learning group and use the previously learned priors as group priors --keep individual level as constraints.
priors_group.MuTheta = p_group.muTheta;
prior_group.SigmaTheta=p_group.SigmaTheta;

gpA=[1 2 4 5 7 10 12 14 16 18 22 24];
gpB=find(~ismember(1:24,gpA));

clearvars options datforce ut
ns=gpA;
for i=1:length(gpA)
    datforce{i}=inpaint_nans(transpose(Force_Data.Reduced.Full_Data{gpA(i)}(46:445)));
    target{i}=(Force_Data.Reduced.Target_Profile{gpA(i)}(46:445));
end

dim.n_t =400;

for i=1:length(gpA)
    options{i}.isYout=isnan(datforce{i});
    options{i}.DisplayWin = 0;
    options{i}.verbose = 1;
    options{i}.dim = dim;
    options{i}.binomial = 0;
    options{i}.inG.statemap = [1;1];
    options{i}.backwardLag = 3 ;
    options{i}.priors=priors;

end

[p_sub_A,o_sub_A,p_group_A,o_group_A] = VBA_MFX(datforce,target,f_fname,g_fname,dim,options,priors_group);



gpB=find(~ismember(1:24,gpA));

clearvars options datforce ut
ns=gpB;
for i=1:length(gpB)
    datforce{i}=inpaint_nans(transpose(Force_Data.Reduced.Full_Data{gpB(i)}(46:445)));
    target{i}=(Force_Data.Reduced.Target_Profile{gpB(i)}(46:445));
end

dim.n_t =400;

for i=1:length(gpB)
    options{i}.isYout=isnan(datforce{i});
    options{i}.DisplayWin = 0;
    options{i}.verbose = 1;
    options{i}.dim = dim;
    options{i}.binomial = 0;
    options{i}.inG.statemap = [1;1];
    options{i}.backwardLag = 3 ;
    options{i}.priors=priors;

end

[p_sub_B,o_sub_B,p_group_B,o_group_B] = VBA_MFX(datforce,target,f_fname,g_fname,dim,options,priors_group);




%%%%%%%%% Do other analyses of Matt's data
%1. psychophysics analysis
%2. use different noise for different states (compare models through wihtin-subject BMA)
%3. Calculate change in motor reponse (n-1 and n for error in n)