function [relative_action] = get_relative_action(X)
%

N = size(X,2);
relative_action = zeros(2,2,N,N);
for j=1:N
    xj = X(:,j);
    for l=(j+1):N
        xl = X(:,l);
        c = reverese(circcorr(xj, xl));
     %  cc = abs(c)/max(abs(c));
        [~, indm] = max(abs(c));
    end
end
end

