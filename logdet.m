function y = logdet(A)
	[L,U] = lu(A); u = diag(U); 
	if prod(sign(u))~=det(L), error('det(A)<=0'), end
	y = sum(log(abs(u)));
