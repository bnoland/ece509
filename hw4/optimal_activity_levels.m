%%
A = [1 2 0 1;
     0 0 3 1;
     0 3 1 1;
     2 1 2 5;
     1 0 3 2];

c_max = repmat(100, 5, 1);
p = [3 2 7 6]';
p_disc = [2 1 4 2]';
q = [4 10 5 10]';

n = 4;

%%
cvx_begin
  variables x(n) u(n);
  maximize sum(u);
  subject to
    x >= 0;
    A * x <= c_max;
    for j = 1:n
      p(j) * x(j) >= u(j);
      p(j) * q(j) + p_disc(j) * (x(j) - q(j)) >= u(j);
    end
cvx_end

%%
disp('Optimal activity levels:')
disp(x)
disp('Associated revenues:')
disp(u)
disp(['Total revenue: ', num2str(sum(u)), newline])
disp('Average price per unit:')
disp(u ./ x)
