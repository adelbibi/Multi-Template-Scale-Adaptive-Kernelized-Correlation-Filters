function [positions, rect_results, time] = tracker(video_path, img_files, pos, target_sz, ...
    padding, kernel, lambda, output_sigma_factor, interp_factor, cell_size, ...
    features, show_visualization,mu,maxitr,mu_inc,step_sc)
%TRACKER Kernelized/Dual Correlation Filter (KCF/DCF) tracking.
%   This function implements the pipeline for tracking with the KCF (by
%   choosing a non-linear kernel) and DCF (by choosing a linear kernel).
%
%   It is meant to be called by the interface function RUN_TRACKER, which
%   sets up the parameters and loads the video information.
%
%   Parameters:
%     VIDEO_PATH is the location of the image files (must end with a slash
%      '/' or '\').
%     IMG_FILES is a cell array of image file names.
%     POS and TARGET_SZ are the initial position and size of the target
%      (both in format [rows, columns]).
%     PADDING is the additional tracked region, for context, relative to
%      the target size.
%     KERNEL is a struct describing the kernel. The field TYPE must be one
%      of 'gaussian', 'polynomial' or 'linear'. The optional fields SIGMA,
%      POLY_A and POLY_B are the parameters for the Gaussian and Polynomial
%      kernels.
%     OUTPUT_SIGMA_FACTOR is the spatial bandwidth of the regression
%      target, relative to the target size.
%     INTERP_FACTOR is the adaptation rate of the tracker.
%     CELL_SIZE is the number of pixels per cell (must be 1 if using raw
%      pixels).
%     FEATURES is a struct describing the used features (see GET_FEATURES).
%     SHOW_VISUALIZATION will show an interactive video if set to true.
%
%   Outputs:
%    POSITIONS is an Nx2 matrix of target positions over time (in the
%     format [rows, columns]).
%    TIME is the tracker execution time, without video loading/rendering.
%
%   Joao F. Henriques, 2014


%if the target is large, lower the resolution, we don't need that much
%detail
scales = [0.97 0.98 0.99 1 1.01 1.02 1.03];
num_scale = length(scales);
scale_weights = normpdf(linspace(-step_sc,step_sc,num_scale),0,1);
best_scale = 1;

resize_image = (sqrt(prod(target_sz)) >= 100);  %diagonal size >= threshold
if resize_image,
    pos = floor(pos / 2);
    target_sz = floor(target_sz / 2);
end

%window size, taking padding into account
window_sz = floor(target_sz * (1 + padding));
%create regression labels, gaussian shaped, with a bandwidth
%proportional to target size
output_sigma = sqrt(prod(target_sz)) * output_sigma_factor / cell_size;
yf = fft2(gaussian_shaped_labels(output_sigma, floor(window_sz / cell_size)));
y = ifft2(yf);
%store pre-computed cosine window
cos_window = hann(size(yf,1)) * hann(size(yf,2))';

if show_visualization,  %create video interface
    update_visualization = show_video(img_files, video_path, resize_image);
end

%note: variables ending with 'f' are in the Fourier domain.
time = 0;  %to calculate FPS
positions = zeros(numel(img_files), 2);  %to calculate precision
rect_results = zeros(numel(img_files), 4);  %to calculate
size_ver = size(cos_window,1);
size_hor = size(cos_window,2);
    
for frame = 1:numel(img_files),
    %load image
    im = imread([video_path img_files{frame}]);
    if size(im,3) > 1,
        im = rgb2gray(im);
    end
    if resize_image,
        im = imresize(im, 0.5);
    end
    
    tic()
    
    if(frame ==1)
        patch = get_subwindow(im, pos, window_sz);
        xf   = fft2(get_features(patch, features, cell_size, cos_window));
        
        switch kernel.type
            case 'gaussian',
                kf = gaussian_correlation(xf, xf, kernel.sigma);
            case 'polynomial',
                kf = polynomial_correlation(xf, xf, kernel.poly_a, kernel.poly_b);
            case 'linear',
                kf = linear_correlation(xf, xf);
        end
        alphaf = yf ./ (kf + lambda);   %equation for fast training
        model_alphaf = alphaf;
        model_xf = xf;        
    else
        [best_scale,pos,response] = tracker_multi_KCF(im,pos,kernel, ...
            cell_size, features,window_sz, ...
            cos_window,model_xf,model_alphaf,num_scale,scales,scale_weights);

        if (best_scale <0.3)
            best_scale = 0.3;
        elseif (best_scale>3)
            best_scale =3;
        end

        %Update scales
        scales = scales*best_scale/scales(round(end/2));
        [vert_delta, horiz_delta] = find(response == max(response(:)), 1);
        if vert_delta > size_ver / 2,  %wrap around to negative half-space of vertical axis
            vert_delta = vert_delta - size_ver;
        end
        if horiz_delta > size_hor / 2,  %same for horizontal axis
            horiz_delta = horiz_delta - size_hor;
        end
        
        pos = pos + (cell_size * [vert_delta - 1, horiz_delta - 1])*best_scale;
        patch = get_subwindow(im, pos, round(window_sz*best_scale));     
        patch = imresize(patch,window_sz,'nearest');
        
        %training at newly estimated target position
        xf = fft2(get_features(patch, features, cell_size, cos_window));

        switch kernel.type
            case 'gaussian',
                kf = gaussian_correlation(xf, xf, kernel.sigma);
            case 'polynomial',
                kf = polynomial_correlation(xf, xf, kernel.poly_a, kernel.poly_b);
            case 'linear',
                kf = linear_correlation(xf, xf);
        end
        alphaf = yf ./ (kf + lambda);   %equation for fast training
        [~, alphaf_2_new] = retrain_fixed_pt_method(model_xf,model_alphaf,xf,alphaf,y,kernel,lambda,mu,maxitr,mu_inc);
        model_alphaf = alphaf_2_new;
        model_xf = (1 - interp_factor) * model_xf + interp_factor * xf;
        
    end
    
    %save position and timing
    positions(frame,:) = pos;
    time = time + toc();
    box = [pos([2,1]) - target_sz([2,1])/2, target_sz([2,1])];
    rect_results(frame,:)=box;
    
    %visualization
    if show_visualization,
        box = [pos([2,1]) - target_sz([2,1])/2, target_sz([2,1])];
        stop = update_visualization(frame, box);
        if stop, break, end  %user pressed Esc, stop early        
        drawnow
        % 			pause(0.05)  %uncomment to run slower
    end
    
end

if resize_image,
    positions = positions * 2;
    rect_results = rect_results *2;
end
end

