%RUN_DEMO run through this project. 
%   input image is set in `generate_init_data` file

generate_init_data
write_kernel(k, 'true');

verbose = true;

% matlab builtin regularized filter
[Ireg, K] = deblur( Nd, B, true, 'reg', verbose );
write_kernel(K, 'estimated');
imwrite(Ireg, 'images/deblurred_reg.jpg');

% call deconv directly since we already have K

% matlab builtin Richardson-Lucy
Ilucy = deconv(double(Nd), double(B), double(K), 'lucy', verbose);
imwrite(Ilucy, 'images/deblurred_lucy.jpg');

% residual RL algorithm
Ir = deconv(double(Nd), double(B), double(K), 'resRL', verbose);
imwrite(Ir, 'images/deblurred_residual_RL.jpg');

% gain-controlled RL
Ig = deconv(double(Nd), double(B), double(K), 'gcRL', verbose);
imwrite(Ig, 'images/deblurred_gain_controlled_RL.jpg');

% gain-controlled Rl with detail
% no need to call deconv again, because we already have Ir and Ig
% Idetailed = deconv(double(Nd), double(B), double(K), 'detailedRL', verbose);
disp('Calculating Ibar with joint/cross bilateral filtering...');
Ibar = zeros(size(Ir), 'double');
[~, ~, d] = size(Ir);
for i = 1:d
    % joint/cross bilateral filtering
    Ibar(:,:,i) = jbfilter2(Ir(:,:,i), Ig(:,:,i), 5, [1.6, 0.08]);
end
Id = Ir - Ibar; % detail layer Id
imwrite(Id+0.8, 'images/detail_layer.jpg');
Idetailed = Ig + Id;    % final result
imwrite(Idetailed, 'images/deblurred_detail_RL.jpg');


disp('ALL CLEAR.')

%%
figure
m=3;n=4;a=201;b=307;c=179;d=284;
subplot(m,n,1)
imshow(I(a:b,c:d),[]);title('Original Image');
subplot(m,n,1+n*1)
imshow(B(a:b,c:d),[]);title('Blurred Image');
subplot(m,n,1+n*2)
imshow(k,[]);title('true kernel');

subplot(m,n,2)
imshow(N(a:b,c:d),[]);title('Noised Image');
subplot(m,n,2+n*1)
imshow(Nd(a:b,c:d),[]);title('Denoised Image');
subplot(m,n,2+n*2)
imshow(K,[]);title('estimated kernel');

subplot(m,n,3)
imshow(Ireg(a:b,c:d),[]);title('deblurred reg');
subplot(m,n,3+n*1) 
imshow(Ilucy(a:b,c:d),[]);title('deblurred lucy');
subplot(m,n,3+n*2)
imshow(Ir(a:b,c:d),[]);title('deblurred residual RL');

subplot(m,n,4) 
imshow(Ig(a:b,c:d),[]);title('deblurred gain controlled RL');
subplot(m,n,4+n*1)
imshow(Idetailed(a:b,c:d),[]);title('deblurred detail RL');
subplot(m,n,4+n*2)
imshow(Id(a:b,c:d),[]);title('detail image');
