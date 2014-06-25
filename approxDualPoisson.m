function [alpha, sW, L, nlZ, dnlZ] = approxDualPoisson(hyper, covfunc, x, y)

% Approximation to the posterior Gaussian Process by minimization of the 
% KL-divergence. The function takes a specified covariance function (see 
% covFunction.m) and likelihood function (see likelihoods.m), and is designed to
% be used with binaryGP.m. See also approximations.m.
%
% Written by Hannes Nickisch, 2007-03-29

n = size(x,1);
K = feval(covfunc{:}, hyper.cov, x);                % evaluate the covariance matrix

% a) simply start at zero
%alla_init{1} = zeros(2*n,1);                       % stack alpha/lambda together
alla_init{1} = [ones(n,1).*0.5];                       % stack alpha/lambda together
alla_init=alla_init([1]);

for alla_id = 1:length(alla_init)              % iterate over initial conditions

    alla = alla_init{alla_id};

    use_pinv=false; check_cond=true;
    nlZ_old = Inf; nlZ_new = 1e100; it=0;      % make sure the while loop starts

	[alla nlZ_new] = lbfgs(alla, K, y,hyper);  %using L-BFGS to find the opt alla

    % save results
    alla_result{alla_id} = alla;
    nlZ_result( alla_id) = nlZ_new;
end

alla_id = find(nlZ_result==min(nlZ_result)); alla_id = alla_id(1);
alla    = alla_result{alla_id};                            % extract best result

%y_con=y>0
y_con=y;
lambda = alla(1:end  ,1);
alpha  = y_con-lambda;

W  = lambda

% recalculate L
sW = sqrt(W);                     
L  = chol(eye(n)+sW*sW'.*K)                             % L'*L=B=eye(n)+sW*K*sW 

%targe
%tmp=L'\(diag(sW));
%rees=-tmp'*tmp

% bound on neg log marginal likelihood

nlZ = margLik_log(alpha,lambda,K,y,hyper)


%estimate the hpyer parameter
% do we want derivatives?
if nargout >=4                                     

    dnlZ = zeros(size(hyper.cov));                  % allocate space for derivatives
    % parameters after optimization
    
    A = inv(eye(n)+K*W )


	Sigma = A*K

	tt=Sigma-inv(K)

	alpha
	mu = K*alpha

    v=abs(diag(A*K))
	[a,dm,dC] = a_related2(K*alpha,v,y,hyper);

    for j=1:length(hyper.cov)
        dK = feval(covfunc{:},hyper.cov,x,j);
		%           from the paper 
        %           -alpha'*dK*dm +(alpha'*dK*alpha)/2 -diag(A*dK*A')'*dC 
        %           -trace(A'*diag(lambda)*dK) +trace(A*dK*diag(lambda)*A)
		%           Note that lambda == dC
        AdK  = A*dK;
        dnlZ(j) = -(alpha'*dK*(dm-alpha/2) +sum(A.*AdK,2)'*dC                ...
                    +(diag(AdK)'-sum(A'.*AdK,1))*lambda);
    end


	dnlZ = hyper.cov;                                   % allocate space for derivatives

	for j=1:length(hyper.cov)                                    % covariance hypers
		dK = feval(covfunc{:},hyper.cov,x,j)
		%dK = feval(cov{:},hyp.cov,x,[],j);
		AdK = A*dK;
		tmp1=sum(A.*AdK,2)
		tmp2=sum(A'.*AdK,1)

		z = diag(AdK) + sum(A.*AdK,2) - sum(A'.*AdK,1)';
		%dnlZ(j) = alpha'*dK*(alpha/2-df) - z'*dv;
		dnlZ(j) = alpha'*dK*(alpha/2-dm) - z'*dC;
	end


   % dnlZ_lik=zeros(size(hyper.lik));
	%for j=1:length(hyper.lik)                                    % likelihood hypers
		%lp_dhyp = likKL(v,lik,hyper.lik,y,K*alpha,[],[],j);
		%dnlZ_lik(j) = -sum(lp_dhyp);
	%end
	%disp('dnlZ_lik=')
	%sprintf('%.15f\n',dnlZ_lik)

  %for j=1:length(hyp.mean)                                         % mean hypers
	%dm_t = feval(mean{:}, hyp.mean, x, j);
	%dnlZ.mean(j) = -alpha'*dm_t;
  %end

end




%% evaluation of current negative log marginal likelihood depending on the
%  parameters alpha (al) and lambda (la)
function [nlZ] = margLik_log(alpha,lambda,K,y,hpyer)
    % extract single parameters
    % dimensions
    n  = length(y);


    % original variables instead of alpha and la
    VinvK = inv(eye(n)+K*diag(lambda));                          % A:=V*inv(K)
    V     = VinvK*K; V=(V+V')/2;                              % enforce symmetry
    v     = abs(diag(V));             % abs prevents numerically negative values
    m     = K*alpha;

    % calculate alpha related terms we need
    [a] = a_related2(m,v,y,hpyer);

	%res1=trace(VinvK)
	%W = abs(-2*lambda);
    %sW = sqrt(W); L = chol(eye(n)+sW*sW'.*K); 
	%L_inv=L\eye(n);
	%res2=trace(L_inv'*L_inv)
	%Note res1==res2

    %negative Likelihood
    nlZ = -logdet(VinvK)/2 -n/2 +(alpha'*K*alpha)/2 +trace(VinvK)/2;
    nlZ = nlZ-a;




function [alla2 nlZ] = lbfgs(alla, K, y,  hyper)
	optMinFunc = struct('Display', 'FULL',...
    'Method', 'lbfgs',...
    'DerivativeCheck', 'off',...
    'LS_type', 1,...
    'MaxIter', 1000,...
	'LS_interp', 1,...
    'MaxFunEvals', 1000000,...
    'Corr' , 100,...
    'optTol', 1e-15,...
    'progTol', 1e-15);
	[alla2, nlZ] = minFunc(@dual_lik, alla, optMinFunc, K, y, hyper);


function [nlZ,dnlZ] = dual_lik(alla,K,y,hpyer)
	lambda=alla(1:end,1);
	assert(min(lambda)>0);
	%assert(max(lambda)<1);

	raw_y=y;
	%y=(y>0);

	assert(min(y)>=0)

	%a=sum(lambda.*log(lambda)+(1-lambda).*log(1-lambda));
	a=sum(lambda.*(log(lambda)-1));

	A_lambda=inv(K)+diag(lambda);

	n=size(K,1);

	nlZ=0.5*((lambda-y)'*K*(lambda-y)-logdet(A_lambda))+a



	alpha1=lambda-y;
	kk = margLik_log(-alpha1,lambda,K,raw_y,hpyer)


	h=-K*(alpha1);
	V=inv(A_lambda);
	p=diag(V);

	trace1=trace(K\V);
	exact=sum(y.*h-exp(h+0.5.*p)-gammaln(y+1));

	kk=0.5*(n-logdet(K)-trace1+logdet(V)-alpha1'*K*alpha1)+exact;
	kk=-kk



	if nargout>1
		%d_lambda=log(lambda)-log(1-lambda);
		d_lambda=log(lambda);

		A_lambda_inv=inv(A_lambda);
		dnlZ=K*(lambda-y)-0.5*diag(A_lambda_inv)+d_lambda;
	end

%% log(det(A)) for det(A)>0
function y = logdet(A)
    % 1) y=det(A); if det(A)<=0, error('det(A)<=0'), end, y=log(y); 
    %   => naive implementation, not numerically stable
    % 2) U=chol(A); y=2*sum(log(diag(U))); 
    %   => fast, but works for symmetric p.d. matrices only
    % 3) det(A)=det(L)*det(U)=det(L)*prod(diag(U)) 
    %   => logdet(A)=log(sum(log(diag(U)))) if det(A)>0
    [L,U]=lu(A); 
    u=diag(U); 
    if prod(sign(u))~=det(L)
        error('det(A)<=0')
    end
    y=sum(log(abs(u))); % slower, but no symmetry needed 
    % 4) d=eig(A); if prod(sign(d))<1, error('det(A)<=0'), end
    %    y=sum(log(d)); y=real(y); 
    %   => slowest

%% compute all terms related to a
% derivatives w.r.t diag(V) and m, 2nd derivatives w.r.t diag(V) and m
function [a,dm,dV,d2m,d2V,dmdV]=a_related2(m,v,y,hyper)
	[a,dm,dV] = likKL(v, hyper.lik,y,m);
	a = sum(a);

%add

function [ll,df,dv] = likKL(v, varargin)
  f = varargin{3};                               % obtain location of evaluation
  y = varargin{2}
  
  elg = exp(f+0.5*v);
  ll = f.*y - elg - gammaln(y+1);

  df= y-elg;
  dv= 0.5*(-elg);

