function [Qhat, Rhat] = mgs(A);

    [m, n] = size(A);
    Qhat = zeros(m,n);
    Rhat = zeros(n,n);
    for i = 1:n
        v = A(:,i);
        for j = 1:i-1
            Rhat(j,i) = Qhat(:,j)'*v;
            v = v - Qhat(:,j)*Rhat(j,i);
        end
        Rhat(i,i) = norm(v);
        Qhat(:,i) = v/Rhat(i,i);
    end
    
    accuracy = norm(A - Qhat*Rhat);
    disp(accuracy);
    
    deviation = norm(Qhat*Qhat' - eye(m));
    disp(deviation);
    
end