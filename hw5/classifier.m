%%
file = "hand.jpg";
image_data = imread(file);

%%
image(image_data);

%%
N = 2;  % Number of training patches to sample from background.

training_patches = cell(N+1, 1);
for i = 1:N
  training_patches{i} = input_training_patch(image_data);
end

% Collect a single training patch from the object we're trying to detect.
training_patches{N+1} = input_training_patch(image_data);

training_data = cell2mat(training_patches);

function training_data = input_training_patch(image_data)
  rect = ginput(2);

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

  training_data = [
    repmat(y_coords, 1, length(x_coords))
    repmat(x_coords, 1, length(y_coords))
    patch
  ]';
end
