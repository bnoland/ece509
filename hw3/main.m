%%

%%%%%%%%%%%%%%%%%%
% Image correction
%%%%%%%%%%%%%%%%%%

clear

%%
% Read in the image data and display it.

image_data = double(imread('cguitar.tif'));
image(image_data)
colormap('gray')
title('Uncorrected image')
xlabel('x')
ylabel('y')

%%
% The optimization problem will be solved using the upper left corner of
% the image, so extract it.

x_bound = 49;
y_bound = 249;

upper_left_corner = image_data(1:y_bound, 1:x_bound);

%%
% Solve the optimization problem.

i_mat = repmat((1:y_bound)', 1, x_bound);
j_mat = repmat(1:x_bound, y_bound, 1);

cvx_begin
  variables a b c;
  R = a * i_mat + b * j_mat + c;
  
  minimize norm(upper_left_corner - R * 255, 'fro');
  subject to
    abs(R) <= 1;
cvx_end

disp('Optimal values:');
disp(['a = ', num2str(a)]);
disp(['b = ', num2str(b)]);
disp(['c = ', num2str(c)]);

%%
% Compute the corrected image and display it.

[y_bound, x_bound] = size(image_data);

i_mat = repmat((1:y_bound)', 1, x_bound);
j_mat = repmat(1:x_bound, y_bound, 1);

R = a * i_mat + b * j_mat + c;
corrected_image = image_data ./ R;

figure;
image(corrected_image)
colormap('gray')
title('Corrected image')
xlabel('x')
ylabel('y')

%%

%%%%%%%%%%%%%%%%
% Spline fitting
%%%%%%%%%%%%%%%%

clear

%%
% Read in the image data and select the input points on the image.

M = 10;  % Number of segments
K = 5;  % Number of input points per segment (constant for simplicity)
N = K * M;  % Number of input points

image_data = imread('curvedriver.jpg');
fig = figure;
image(image_data)
hold on

figure(fig);  % Ensure this figure is focused for ginput()
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

