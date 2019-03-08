
%%
% Select the input points on the image.

M = 1;  % Number of segments
N = 4 * M;  % Number of input points (multiple of M for simplicity)

image_data = imread('curvedriver.jpg');
image(image_data)
hold on

[x, y] = ginput(N);
plot(x, y, '-wx')
hold off

%%
C = mat2cell([x, y], repmat(N / M, 1, M));
celldisp(C);

t = 0:0.01:1;

cvx_begin
  variables a(4) b(4);
  x_fit = a(4) * t.^3 + a(3) * t.^2 + a(2) * t + a(1);
  y_fit = b(4) * t.^3 + b(3) * t.^2 + b(2) * t + b(1);
cvx_end
