function [t_new]=check_valid(x,t,d)
	t_new=t;
	x_check=x+t*d;

	%for unconstrained problems
	%lb=-inf;
	%ub=inf;

	%for logit
	lb=0;
	ub=1;

	left_idx=(x_check<=lb);
	right_idx=(x_check>=ub);

	left=x(left_idx);
	right=x(right_idx);

	%if it is the strictly bound, set lambda=1e-5, else lambda=0
	lambda=1e-5;
	if length(left)>0
		d_left=abs(d(left_idx));
		assert(min(d_left)>0)
		step=min((left-lb)./d_left);
		%strictly bound
		t_new=(1.0-lambda)*min(t_new,step);
	end

	if length(right)>0
		d_right=abs(d(right_idx));
		assert(min(d_right)>0)
		step=min((ub-right)./d_right);
		t_new=(1.0-lambda)*min(t_new,step);
	end

end

