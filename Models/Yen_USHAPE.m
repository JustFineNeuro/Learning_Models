%Models:
%A. Single -state (garbage, could never work without extraneous assumption)

%B. 1 Fast + 2 context-specific slow states - Essentially same as ETM Howard;
%Schweighoffer, J NEurosci; Fine; NeuroImage

%C. Transfer model. (1 fast, + 1 slow state). Assumption is observer
%component essentially sign flips either (the slow alone or the slow +
%fast) --Note a more complicated model would use the history of error
%directions to potentially gate the effective contribution of a fast
%process.

%D. Generalization uses errors in one slow or slow + fast context to update
%other model state


%%%% 1. FIT MODELS




f_fname=@f_twomemorymodel; %evoluation funciton
g_fname=@g_single_obs; %Observation function (This assumes identity mapping);

ut=500*ones(95,1);
mx=peak_mx.peakMx;
y=mx;
ut(find(mx<0))=ut(find(mx<0)).*-1;
ut=ut.*-1;
dim.n_theta = 4; 
dim.n=3; 
 A=ones(95,1)
A(find(sign(ut)==-1))=2;
ut(:,1)=ut(:,1).*-1;
A(:,2)=3
ut(:,2:3)=A;
% % Build priors for model inversion
priors.muX0 = [-300 ;300;0]
priors.SigmaX0 = 1e-0*eye(3);
priors.muTheta = 0*ones(4,1);
priors.SigmaTheta = 1e-1*eye(4);
priors.a_alpha = 1e0;
priors.b_alpha = 1e0;
priors.a_sigma = 1e0;
priors.b_sigma = 1e0;

% 
% options.inG.statemap = sign(ut);
% 
% options.inF.active=sign(ut);
% options.in.active=sign(ut);

options.priors=priors;
% Build options and dim structures for model inversion
% options.inF = in;

dim.n_theta = 4;%Number of evolution equation paramters.
dim.n_phi = 0;
dim.p = 1; disp('not sure what this is, but its for inversion?')
x0 = zeros(dim.n,1);
[posterior,out] = VBA_NLStateSpaceModel(y',ut',f_fname,g_fname,dim,options);

VBA_sigmoid(posterior.muTheta./sqrt((1+0.3*diag(posterior.SigmaTheta))))





%%%%%%%
%BEst model here -- two slow state (same learning rate)
f_fname=@f_twomemorymodel; %evoluation funciton
g_fname=@g_multi_memory_obs; %Observation function (This assumes identity mapping);

ut=500*ones(95,1);
mx=peak_mx.peakMx;
y=mx;
ut(find(mx<0))=ut(find(mx<0)).*-1;
ut=ut.*-1;
dim.n_theta = 4; 
dim.n=3; 
 A=ones(95,1)
A(find(sign(ut)==-1))=2;
ut(:,1)=ut(:,1).*-1;
A(:,2)=3
ut(:,2:3)=A;
% % Build priors for model inversion
priors.muX0 = [-300 ;300;0]
priors.SigmaX0 = 1e-0*eye(3);
priors.muTheta = 0*ones(4,1);
priors.SigmaTheta = 1e-1*eye(4);
priors.a_alpha = 1e0;
priors.b_alpha = 1e0;
priors.a_sigma = 1e0;
priors.b_sigma = 1e0;

% 
% options.inG.statemap = sign(ut);
% 
% options.inF.active=sign(ut);
% options.in.active=sign(ut);

options.priors=priors;
% Build options and dim structures for model inversion
% options.inF = in;

dim.n_theta = 4;%Number of evolution equation paramters.
dim.n_phi = 0;
dim.p = 1; disp('not sure what this is, but its for inversion?')
x0 = zeros(dim.n,1);
[posterior,out] = VBA_NLStateSpaceModel(y',ut',f_fname,g_fname,dim,options);

VBA_sigmoid(posterior.muTheta./sqrt((1+0.3*diag(posterior.SigmaTheta))))


%% Model using same slow + fast state; production occurs through sign-switching of slow state only.

f_fname=@f_twomemorymodel_B; %evoluation funciton
g_fname=@g_multi_memory_obs_B; %Observation function (This assumes identity mapping);

ut=500*ones(95,1);
mx=peak_mx.peakMx;
y=mx;
ut(find(mx<0))=ut(find(mx<0)).*-1;
ut=ut.*-1;
dim.n_theta = 4; 
dim.n=2; 
%  A=ones(95,1)
% A(find(sign(ut)==-1))=2;
% ut(:,1)=ut(:,1).*-1;
% A(:,2)=3
% ut(:,2:3)=A;
% % Build priors for model inversion
priors.muX0 = [-300 ;0]
priors.SigmaX0 = 1e-0*eye(2);
priors.muTheta = 0*ones(4,1);
priors.SigmaTheta = 1e-1*eye(4);
priors.a_alpha = 1e0;
priors.b_alpha = 1e0;
priors.a_sigma = 1e0;
priors.b_sigma = 1e0;

% 

options.priors=priors;
% Build options and dim structures for model inversion
% options.inF = in;

dim.n_theta = 4;%Number of evolution equation paramters.
dim.n_phi = 0;
dim.p = 1; disp('not sure what this is, but its for inversion?')
x0 = zeros(dim.n,1);
[posterior,out] = VBA_NLStateSpaceModel(y',ut',f_fname,g_fname,dim,options);

VBA_sigmoid(posterior.muTheta./sqrt((1+0.3*diag(posterior.SigmaTheta))))

%% MODEL C
%Yen modeling

f_fname=@f_twomemorymodel_C; %evoluation funciton
g_fname=@g_multi_memory_obs; %Observation function (This assumes identity mapping);

ut=500*ones(95,1);
mx=peak_mx.peakMx;
y=mx;
ut(find(mx<0))=ut(find(mx<0)).*-1;
ut=ut.*-1;
dim.n_theta = 5; 
dim.n=3; 
 A=ones(95,1)
A(find(sign(ut)==-1))=2;
ut(:,1)=ut(:,1).*-1;
A(:,2)=3
ut(:,2:3)=A;
% % Build priors for model inversion
priors.muX0 = [-300 ;300;0]
priors.SigmaX0 = 1e-0*eye(3);
priors.muTheta = 0*ones(5,1);
priors.SigmaTheta = 1e-1*eye(5);
priors.a_alpha = 1e0;
priors.b_alpha = 1e0;
priors.a_sigma = 1e0;
priors.b_sigma = 1e0;

% 
% options.inG.statemap = sign(ut);
% 
% options.inF.active=sign(ut);
% options.in.active=sign(ut);

options.priors=priors;
% Build options and dim structures for model inversion
% options.inF = in;

dim.n_phi = 0;
dim.p = 1; disp('not sure what this is, but its for inversion?')
x0 = zeros(dim.n,1);
[posterior,out] = VBA_NLStateSpaceModel(y',ut',f_fname,g_fname,dim,options);

VBA_sigmoid(posterior.muTheta./sqrt((1+0.3*diag(posterior.SigmaTheta))))

%% SAME MODELS WITH DECAY ONLY DURING ACTIVE CONTEXT
%A2. 

%B2. 

%C2. 

%D2.

%%%% 2. TEST MODEL IDENTIFIABILITY TO ASSESS IDEAL PARADIGMS: 
%TODO. find average parameters across all subjects and use as simulation...
%%%%%%%
