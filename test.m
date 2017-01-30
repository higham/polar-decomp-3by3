function test
%TEST   Test the function polar_quaternion.

fprintf('Testing function polar_quaternion\n')
fprintf('Test statistics should all be of order %8.2e\n\n', eps)

fprintf('Test (5.1) from paper:\n')
A = [0.1 0.2 0.3
     0.1 0.1 0.0
     0.3 0.2 0.1];
[Q,H] = polar_quaternion(A);
check(A,Q,H)

fprintf('Test (5.2) from paper:\n')
for y = sqrt([1 1e-4 1e-8 1e-12 1e-16])
    A = (1/1275)*(...
                  [720 650 710
                  396 -145 178
                  972 610 -529]*y + ...
                  [-25 300 300
                  70 -840 -840
                  -10 120 120]);
    [Q,H] = polar_quaternion(A);
    check(A,Q,H)
    % Check relative error in Q.
    [U,S,V] = svd(A); Q1 = U*V';
    kappa = sqrt( (1+2*y^2)/(3*y^2) ); % Condition numner of U.
    fprintf('Scaled rel error in Q  = %8.2e\n',...
             norm(Q - Q1)/(norm(Q1)*kappa));
end
end

function check(A,Q,H)
res = norm(A - Q*H)/norm(A);
orth = norm(Q'*Q-eye(size(Q)));
fprintf('rel resid = %8.2e,  orthogonality = %8.2e\n', res, orth)
end
