function [val, grad] = LS_cost_grad(x, rho, mu, M, lambda, sigma)
%
% Cost function that includes gradient for Matlab's solver
% rho and x are the variables
% mu, lambda, sigma, and  P are the parameters
%
%============== Dihedral Version ==============

% main objects
L = length(x);
x = x(:);
p = rho(1:L);
q = rho((1+L):end);

% auxiliary matrices
circ = @(v) toeplitz([v(1); v(end:-1:2)], v)';
C_x = circ(x);
M_p = circ(p);                      % C_x*p  ==  M_p*x
M_q = reverse(eye(L))*(circ(q))';   % C_x'*q ==  M_q*x

% first moment
inner1 = C_x*p  + C_x'*q - mu(:);  
term1  = norm(inner1)^2;
% second moment
inner2 = triu(C_x*diag(p)*C_x' + C_x'*diag(q)*C_x - M + sigma^2*eye(L));
term2  = norm(inner2,'fro')^2;
val    = term1 + lambda*term2;

% gradient calculation (if needed)
if nargout>1
    % first moment part
    rho_p_part1(1:L) = 2*C_x'*(inner1); 
    rho_q_part1(1:L) =  2*C_x*(inner1);
    x_part1   =  2*(M_p'+M_q')*(inner1); 
    
    % second moment part
    rho_p_part2 = zeros(1,L);
    rho_q_part2 = zeros(1,L);
    x_part2   = 0;
    pp = p(end:-1:1);  
    for j=1:L
        for m=j:L
            rho_q_part2 = rho_q_part2 + 2*inner2(j,m)*( (C_x(:,j)).*(C_x(:,m)) )';
            rho_p_part2 = rho_p_part2 + 2*inner2(j,m)*( (C_x(j,:).').*(C_x(m,:).') )';           
            
            % xpart -- auxiliaries from the quadratic form
            AA = C_x(:,mod(j-m,L)+1);
            BB = C_x(:,mod(m-j,L)+1);                       
            pAA = circshift(pp,mod(j-1,L)+1);
            pBB = circshift(pp,mod(m-1,L)+1);                  
            qAA = circshift(q,mod(-m,L)+1);
            qBB = circshift(q,mod(-j,L)+1);            
            added_part_q = qAA.*AA+qBB.*BB;
            added_part_p = pAA.*AA+pBB.*BB;
            % the final derivative
            x_part2 = x_part2 + 2*inner2(j,m)*(added_part_p)+ 2*inner2(j,m)*(added_part_q);            
        end
    end    
    % summary
    rho_part1 = [rho_p_part1(:); rho_q_part1(:)];
    rho_part2 = [rho_p_part2(:); rho_q_part2(:)];
    % unifing the two parts of cost function
    rho_part  = rho_part1 + lambda*rho_part2(:);
    x_part    = x_part1 + lambda*x_part2(:);
    % unifing the gradient (signal and dist.)
    grad = [x_part(:) ; rho_part(:)];
end
end
