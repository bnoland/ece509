
%%
% Select the input points on the image.

M = 1;  % Number of segments
K = 10;  % Number of input points per segment (constant for simplicity)
N = K * M;  % Number of input points

image_data = imread('curvedriver.jpg');
image(image_data)
hold on

[x, y] = ginput(N);
plot(x, y, '-wx')
%hold off

%%
% Split the input points into collections corresponding to segments.

C = mat2cell([x, y], repmat(K, 1, M));
celldisp(C);

%%
% Solve the optimization problem.

% XXX: Solve for a single segment for now.

x = C{1}(:,1);
x_start = x(1);
x_end = x(K);
t = linspace(0, 1, K)';
X = [ones(K, 1), t, t .^ 2, t .^ 3];

cvx_begin
  variables a(4);
  minimize norm(x - X * a);
cvx_end

y = C{1}(:,2);
y_start = y(1);
y_end = y(K);
t = linspace(0, 1, K)';
X = [ones(K, 1), t, t .^ 2, t .^ 3];

cvx_begin
  variables b(4);
  minimize norm(y - X * b);
cvx_end

%%
t = linspace(0, 1, 100)';
X = [ones(100, 1), t, t .^ 2, t .^ 3];
x = X * a;
y = X * b;

plot(x, y)
