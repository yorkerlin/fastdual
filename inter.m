%% compute all terms related to a
% derivatives w.r.t diag(V) and m, 2nd derivatives w.r.t diag(V) and m
function [a,dm,dV,d2m,d2V,dmdV]=inter(m,v,y,lik)
	N = 20;                                % number of hermite quadrature points
	[f,w]  = gauher(N)            % location and weights for hermite quadrature
    %[f,w]  = hermquad(N)            % location and weights for hermite quadrature
    f_dV   = f.^2-1;
    f_dm   = f;
    f_d2V  = f.^4-6*f.^2+3;
    f_dmdV = f.^3-3*f;
    SumLog_lam = zeros(size(f)); % init accumulator
    if nargout>2, dV  = zeros(size(m)); dm   = dV;  end            % init result
    if nargout>5, d2V = zeros(size(m)); dmdV = d2V; end            % init result
	
	%out2(ok) = -log(1+exp(-yf(ok))); 
	save =[]
    for i=1:length(y)
        [dummy,log_lam] = feval( lik, y(i), (sqrt(v(i))*f+m(i)) );   % squashing 
		log_lam
        SumLog_lam = SumLog_lam+log_lam;                          % accumulation
		save = [save log_lam];
        if nargout>2                                            % do integration
            dV(i)   = (w'*(log_lam.*f_dV  )) / (2*v(i)      );
            dm(i)   = (w'*(log_lam.*f_dm  )) /    v(i)^(1/2) ;
            if nargout>5
                d2V(i)  = (w'*(log_lam.*f_d2V )) / (4*v(i)^2    );
                dmdV(i) = (w'*(log_lam.*f_dmdV)) / (2*v(i)^(3/2));
            end
        end
    end
    a = w'*SumLog_lam;                                          % do integration
	aa = save'*w;
	sprintf('%.13f\n',aa)
    if nargout>5, d2m= 2*dV; end
