% loading initial image in 
load gatlin.mat;
colormap(map);
image(X);

[U,S,V] = svd(X); 
sing_vals = diag(S);
seq = 1:480; 
plot(seq,log(sing_vals)) 
title("x Values vs Log of Singular Values of X matrix")
xlabel("x values") 
ylabel('log of singular values')



[m,n] = size(X); 
 % setting the values of k for approximation
% error_table = zeros([size(values,2),6]);initializing table to have 6 spots

abs_error = zeros(5,1);
relative_error = zeros(5,1);
sigma1_rat = zeros(5,1);
sigmak_rat = zeros(5,1);
comp_rat = zeros(5,1); 

for i = 1:8
    values = [1,5,10,20,40,60, 100, 200];
    [A_k, approx_error, sigma] = bestApprox(X,values(i));
    abs_error(i) = approx_error;
    relative_error(i) = approx_error/norm(X,2);
   sigma1_rat(i) = sigma/S(1,1);
   sigmak_rat(i) = sigma/S(values(i),values(i)); 
    comp_rat(i) = (m + n + 1)*values(i)/(m*n);
    colormap(map);
    figure; 
    %imshow(A_k);
    image(A_k); 
end 




error_table = [abs_error, relative_error, sigma1_rat, sigmak_rat, comp_rat];
disp(array2table(error_table))


 
 