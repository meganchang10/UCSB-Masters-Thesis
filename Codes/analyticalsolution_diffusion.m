function exact = analyticalsolution(m,x0,xf,y0,yf,tfinal,exact_sol,phi_initial)
n = m;
exact = zeros(m,n);
x = linspace(x0,xf,m);
y = linspace(y0,yf,n);
t = tfinal;

phi = zeros(m,n);
for j = 1:n
    for i = 1:m
        phi(i,j) = phi_initial(x(i),y(j),t);
    end
end

for j = 1:n
    for i = 1:m
         if phi(i,j) <= 0
            exact(i,j)= exact_sol(x(i),y(j),t);
         end
    end
end

end