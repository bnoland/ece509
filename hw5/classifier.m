%%
file = "apple.jpg";
image_data = imread(file);

%%
image(image_data);

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

% Compute the support vector classifier.
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

class = zeros(height, width);

%%
N = size(Z, 1);

% Classify each pixel in the image.
for i = 1:N
  obs = Z(i,:);
  x = obs(1);
  y = obs(2);
  class(y, x) = get_class(obs, a, b);
end

image(class);

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
% vector classifier defined by `a` and `b`.
function output = get_class(obs, a, b)
  if obs * a - b > 0
    output = 0;  % background
  else
    output = 255;  % object
  end
end