function u = Constant_Extrapolation(u,phi,m,n,dx,dy,iterations)

% Similar to the third order extrapolation method's implementation, we call
% upon the iterative solver. However in this case, we are only solving the
% first PDE so that we get a first order accurate result or constant
% extrapolation. See Third Order Extrap for further details on H_0

rhs  = zeros(m,n);
H_0  = ones(m,n);

for j = 1:n
    for i = 1:m
        if phi(i,j) < 0
            H_0(i,j) = 0;
        end
    end
end

u = iterative_solver(u,rhs,H_0,phi,dx,dy,m,n,iterations);

end