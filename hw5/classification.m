%%
object_type = 'hand';  % either 'hand' or 'apple'
file = strcat(object_type, '.jpg');
image_data = imread(file);

%%
image(image_data);
title('Original image');

%%
N = 3;  % Number of training patches to sample from background.

bg_patches = cell(N, 1);
for i = 1:N
  rect = ginput(2);
  bg_patches{i} = get_patch(image_data, rect);
end

X = cell2mat(bg_patches);

% Collect a single training patch from the object we're trying to detect.
rect = ginput(2);
Y = get_patch(image_data, rect);

%%
N = size(X, 1);
M = size(Y, 1);

gamma = 0.1;

% Train the support vector machine.
cvx_begin
  variables a(5) b u(N) v(M);
  
  minimize norm(a) + gamma * (sum(u) + sum(v));
  subject to
    for i = 1:N
      X(i,:) * a - b >= 1 - u(i);
    end
    for i = 1:M
      Y(i,:) * a - b <= -(1 - v(i));
    end
    u >= 0;
    v >= 0;
cvx_end

%%
[width, height, depth] = size(image_data);
rect = [1 1; height width];
Z = get_patch(image_data, rect);
N = size(Z, 1);

%%
% Classify each pixel in the image using the SVM.
class = zeros(height, width);
for i = 1:N
  obs = Z(i,:);
  x = obs(1);
  y = obs(2);
  class(y, x) = get_svm_class(obs, a, b);
end

figure;
image(class);
title('SVM pixel classification');
colormap(gray(256));

%%
% Classify each pixel in the image using color thresholding.
class = zeros(height, width);
for i = 1:N
  obs = Z(i,:);
  x = obs(1);
  y = obs(2);
  class(y, x) = get_color_threshold_class(obs, object_type);
end

figure;
image(class);
title('Color thresholding pixel classification');
colormap(gray(256));

%%
% Project the data onto an appropriate 2D subspace and plot it.

perp = null(a');
v = perp(:,1);  % get a vector perpendicular to a
basis = [a v];  % we want to project onto space with this basis
P = basis * ((basis' * basis) \ basis'); % projection matrix

N = size(X, 1);
M = size(Y, 1);
W = P * [X; Y]';

coeffs = basis \ W; % coefficients in the basis
coeffs = coeffs';

figure;
scatter(coeffs(1:N,1), coeffs(1:N,2), 5, 'red', 'filled');
hold on;
scatter(coeffs(N+1:N+M,1), coeffs(N+1:N+M,2), 5, 'blue', 'filled');

x = b / norm(a)^2 * [1 1];
y = ylim;
line(x, y);

title('2D projection of training data');

%%
% Extracts a patch of data from the image `image_data' defined by the
% rectangle `rect', along with the associated pixel coordinates, and puts
% the data in a format suitable for working with the classifier.
function output = get_patch(image_data, rect)
  lim(1) = min(floor(rect(1)), floor(rect(2)));  % x min
  lim(2) = min(floor(rect(3)), floor(rect(4)));  % y min
  lim(3) = max(ceil(rect(1)), ceil(rect(2)));  % x max
  lim(4) = max(ceil(rect(3)), ceil(rect(4)));  % y max

  y_coords = lim(2):lim(4);
  x_coords = lim(1):lim(3);
  patch_size = length(y_coords) * length(x_coords);

  patch = image_data(y_coords, x_coords, :);
  patch = permute(patch, [3 1 2]);
  patch = double(patch);
  patch = reshape(patch, 3, patch_size);

  output = [
    repelem(x_coords, length(y_coords))
    repmat(y_coords, 1, length(x_coords))
    patch
  ]';
end

% Returns the class of an observation `obs` as computed by the support
% vector machine defined by `a` and `b`.
function output = get_svm_class(obs, a, b)
  if obs * a - b > 0
    output = 0;  % background
  else
    output = 255;  % object
  end
end

% Returns the class of an observation `obs` as computed using color
% thresholding based on the object type `object_type` (either 'hand' or
% 'apple').
function output = get_color_threshold_class(obs, object_type)
  r = obs(3);
  g = obs(4);
  b = obs(5);
  
  if strcmp(object_type, 'apple')
    if r < 150
      output = 0; % background
    else
      output = 255; % object
    end
  elseif strcmp(object_type, 'hand')
    if abs(r - 245) <= 20 && abs(g - 245) <= 20 && abs(b - 220) <= 20
      output = 255; % object
    else
      output = 0; % background
    end
  end
end
