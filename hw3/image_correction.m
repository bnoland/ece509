
%%
% Read in the image data and display it.

image_data = double(imread('cguitar.tif'));
image(image_data)
colormap('gray')

%%
% The optimization problem will be solved using the upper left corner of
% the image, so extract it.

x_bound = 49;
y_bound = 249;

upper_left_corner = image_data(1:y_bound, 1:x_bound);
%image(upper_left_corner)

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

image(corrected_image)
colormap('gray')
