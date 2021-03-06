
%%
% Read in the image data and select the input points on the image.

M = 10;  % Number of segments
K = 5;  % Number of input points per segment (constant for simplicity)
N = K * M;  % Number of input points

image_data = imread('curvedriver.jpg');
image(image_data)
hold on

[x, y] = ginput(N);
plot(x, y, '-wx')
title('Cubic spline fit')
xlabel('x')
ylabel('y')

%%
% Split the input points into collections corresponding to segments.

C = mat2cell([x, y], repmat(K, 1, M));
celldisp(C);

%%
% Solve the optimization problem.

t = linspace(0, 1, K)';
A = [ones(K, 1), t, t .^ 2, t .^ 3];

% Model the x coordinates.
cvx_begin
  variable a(4,M);
  
  loss = 0;
  for m = 1:M
    x = C{m}(:,1);
    loss = loss + norm(x - A * a(:,m));
  end
  
  minimize loss
  subject to
    % Smoothness constraints.
    for m = 1:(M-1)
      a(1,m) + a(2,m) + a(3,m) + a(4,m) == a(1,m+1);
      a(2,m) + 2 * a(3,m) + 3 * a(4, m) == a(2,m+1);
      2 * a(3, m) + 6 * a(4,m) == 2 * a(4,m+1);
    end
cvx_end

% Model the y coordinates.
cvx_begin
  variable b(4,M);
  
  loss = 0;
  for m = 1:M
    y = C{m}(:,2);
    loss = loss + norm(y - A * b(:,m));
  end
  
  minimize loss
  subject to
    % Smoothness constraints.
    for m = 1:(M-1)
      b(1,m) + b(2,m) + b(3,m) + b(4,m) == b(1,m+1);
      b(2,m) + 2 * b(3,m) + 3 * b(4, m) == b(2,m+1);
      2 * b(3, m) + 6 * b(4,m) == 2 * b(4,m+1);
    end
cvx_end

%%
% Plot the resulting cubic spline.

t = linspace(0, 1, 100)';
A = [ones(size(t, 1), 1), t, t .^ 2, t .^ 3];

for m = 1:M
  x = A * a(:,m);
  y = A * b(:,m);
  plot(x, y)
end
