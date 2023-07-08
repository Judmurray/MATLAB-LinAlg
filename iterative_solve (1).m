function iterative_solve(tol)
 T = [-4, 1, 0, 0, 0;
1, -5, 1, 0, 0;
0, 1, -4, 1, 0;
0, 0, 1, -5, 1;
0, 0, 0, 1, -4]; 
b = [0 6 3 1 1 6 -1 4 0 6]'; 
x0 = b; 
a = [T, eye(5); eye(5), T]; 
n = length(b);

m = diag(a);
M = diag(m);
N = M-a; 
spec_rad_jacobi = max(abs(eig(inv(M)*N)))

for j = 1 : n
     x(j) = ((b(j) - a(j,[1:j-1,j+1:n]) * x0([1:j-1,j+1:n])) / a(j,j)); 
end
x1 = x';
k = 1;
while norm(x1-x0,2) > tol
    for j = 1 : n
     x_j(j) = ((b(j) - a(j,[1:j-1,j+1:n]) * x1([1:j-1,j+1:n])) / a(j,j));
    end
    x0 = x1;
    x1 = x_j';
    k = k + 1;
end
disp('iterations and x from jacobi  =')
k
x1


 %gauss seidel starts 
d=b;
es=5;
for i=1:n
    d(i)=b(i)/a(i,i);
end
C=a;

for i=1:n
    C(i,:)=1/a(i,i)*a(i,:);
    C(i,i)=0;
end

ea=100;
x_i=zeros(n,1);
while ea>es
    xold=x_i;
    for i=1:n
        x_i(i)=d(i)-C(i,:)*x_i;
        
    end

ea=norm(x_i-xold,2)/norm(x_i,2)*100;
x2 = x_i;
end
 disp('x from Gauss-Seidel = ')
disp(x2)
 disp('actual error as percent = ')
 disp(ea)
 
 P = tril(a); 
 Q = P-a;
 spec_radius_GS = max(abs(eig(inv(P)*Q)))
 
 d=b;
es=5;
for i=1:n
    d(i)=b(i)/a(i,i);
end
C=a;

for i=1:n
    C(i,:)=1/a(i,i)*a(i,:);
    C(i,i)=0;
end

ea=100;
x_i=zeros(n,1);
w = 1.1174 % from the SOR graph 
while ea>es
    xold=x_i;
    for i=1:n
        x_i(i)=d(i)-C(i,:)*x_i;
        x_i(i)=w*x_i(i)+(1-w)*xold(i);
        
    end

ea=norm(x_i-xold,2)/norm(x_i,2)*100;
x2 = x_i;
end
 disp('x from Gauss-Seidel = ')
disp(x2)
 disp('actual error as percent = ')
 disp(ea)
 
 P = tril(a); 
 Q = P-a;
 spec_radius_GS = max(abs(eig(inv(P)*Q)))
end

