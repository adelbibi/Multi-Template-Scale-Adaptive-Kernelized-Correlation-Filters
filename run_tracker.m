%
%  High-Speed Tracking with Kernelized Correlation Filters
%
%  Joao F. Henriques, 2014
%  http://www.isr.uc.pt/~henriques/
%
%  Main interface for Kernelized/Dual Correlation Filters (KCF/DCF).
%  This function takes care of setting up parameters, loading video
%  information and computing precisions. For the actual tracking code,
%  check out the TRACKER function.
%
%  RUN_TRACKER
%    Without any parameters, will ask you to choose a video, track using
%    the Gaussian KCF on HOG, and show the results in an interactive
%    figure. Press 'Esc' to stop the tracker early. You can navigate the
%    video using the scrollbar at the bottom.
%
%  RUN_TRACKER VIDEO
%    Allows you to select a VIDEO by its name. 'all' will run all videos
%    and show average statistics. 'choose' will select one interactively.
%
%  RUN_TRACKER VIDEO KERNEL
%    Choose a KERNEL. 'gaussian'/'polynomial' to run KCF, 'linear' for DCF.
%
%  RUN_TRACKER VIDEO KERNEL FEATURE
%    Choose a FEATURE type, either 'hog' or 'gray' (raw pixels).
%
%  RUN_TRACKER(VIDEO, KERNEL, FEATURE, SHOW_VISUALIZATION, SHOW_PLOTS)
%    Decide whether to show the scrollable figure, and the precision plot.
%
%  Useful combinations:
%  >> run_tracker choose gaussian hog  %Kernelized Correlation Filter (KCF)
%  >> run_tracker choose linear hog    %Dual Correlation Filter (DCF)
%  >> run_tracker choose gaussian gray %Single-channel KCF (ECCV'12 paper)
%  >> run_tracker choose linear gray   %MOSSE filter (single channel)
%
function [precision, overlap, fps, fn] = run_tracker(dataset, video, kernel_type, feature_type, show_visualization, show_plots,...
    padding,lambda, output_sigma_factor,interp_factor,kernel_sigma,cell_size,hog_orientations,mu,maxitr,mu_inc,step_sc)
	
	%path to the videos (you'll be able to choose one with the GUI).
	base_path = strcat(pwd,'\');
	result_path = 'results/';
	tracker_name = 'KCF_MTSA';
	thresholdPrecision = 20;
	thresholdOverlap = 0.5;

%default settings
if nargin < 1, dataset = 'videos'; end
if nargin < 2, video = 'choose'; end
if nargin < 3, kernel_type = 'gaussian'; end
if nargin < 4, feature_type = 'hog'; end
if nargin < 5, show_visualization = ~strcmp(video, 'all'); end
if nargin < 6, show_plots = ~strcmp(video, 'all'); end
if nargin < 7, padding = 1.7; end %extra area surrounding the target
if nargin < 8, lambda = 1e-3; end %regularization
if nargin < 9, output_sigma_factor = 0.1; end %spatial bandwidth (proportional to target)
if nargin < 16, mu_inc = 2; end
if nargin < 17, step_sc = 1e-2; end

%parameters based on the chosen kernel or feature type
kernel.type = kernel_type;
features.gray = false;
features.hog = false;
features.hogcolor = false;

	
	switch feature_type
    case 'gray',
        if nargin < 10, interp_factor = 0.075; end
        %interp_factor = 0.075;  %linear interpolation factor for adaptation
        
        if nargin < 11, kernel_sigma = 0.2; end %gaussian kernel bandwidth
        kernel.sigma = kernel_sigma;
        
        kernel.poly_a = 1;  %polynomial kernel additive term
        kernel.poly_b = 7;  %polynomial kernel exponent
        
        if nargin < 12, cell_size = 1; end
        
        hog_orientations = 0;
        if nargin < 14, mu = 1e-6; end
        if nargin < 15, maxitr = 10; end

        
        features.gray = true;
        
    case 'hog',
        if nargin < 10, interp_factor = 0.01; end
        
        if nargin < 11, kernel_sigma = 0.5; end %gaussian kernel bandwidth
        kernel.sigma = kernel_sigma;
        
        kernel.poly_a = 1;
        kernel.poly_b = 9;
        
        if nargin < 12, cell_size = 4; end
        
        if nargin < 13, hog_orientations = 9; end
        features.hog_orientations = hog_orientations;

        if nargin < 14, mu = 1e-6; end
        if nargin < 15, maxitr = 10; end
        
        features.hog = true;
        
    case 'hogcolor',
        if nargin < 10, interp_factor = 0.01; end
        
        if nargin < 11, kernel_sigma = 0.5; end %gaussian kernel bandwidth
        kernel.sigma = kernel_sigma;
        
        kernel.poly_a = 1;
        kernel.poly_b = 9;
        
        if nargin < 12, cell_size = 4; end
        
        if nargin < 13, hog_orientations = 9; end
        features.hog_orientations = hog_orientations;

        if nargin < 14, mu = 1e-6; end
        if nargin < 15, maxitr = 10; end
        
        features.hogcolor = true;
        
    otherwise
        error('Unknown feature.')
        
end


	assert(any(strcmp(kernel_type, {'linear', 'polynomial', 'gaussian'})), 'Unknown kernel.')


	switch video
    
    case 'choose',
        %ask the user for the video, then call self with that video name.
        video = choose_video([base_path dataset '/']);
        if ~isempty(video),
            [precision, overlap, fps, fn] = run_tracker(dataset, video, kernel_type, feature_type, show_visualization, show_plots,...
                padding,lambda, output_sigma_factor,interp_factor,...
                kernel_sigma,cell_size,hog_orientations,mu,maxitr,mu_inc,step_sc)
        end
        
    case 'all',
        %all videos, call self with each video name.
        
        %Print parameters to screen
        fprintf('----------------------------------------------------------------------------------\n');
        fprintf([' Dataset ......................... ', dataset '\n']);
        fprintf([' Kernel .......................... ', kernel_type '\n']);
        fprintf([' Feature ......................... ', feature_type  '\n']);
        fprintf(' Padding ......................... %.3g\n', padding);
        fprintf(' Lambda .......................... %.1e\n', lambda);
        fprintf(' Interpolation Factor............. %.2g\n', interp_factor);
        fprintf(' Output Sigma Factor ............. %.2g\n', output_sigma_factor);
        fprintf(' Kernel Sigma .................... %.2g\n', kernel_sigma);
        fprintf(' Cell Size ....................... %u\n', cell_size);
        fprintf(' Hog Orientations ................ %u\n', hog_orientations);
        fprintf(' Mu ................ %u\n', mu); 
        fprintf(' Max Iterations ................ %u\n', maxitr);
        fprintf(' Increase in Mu ................ %u\n', mu_inc);
        fprintf(' Scale Step ................ %u\n', step_sc); 		
        fprintf('----------------------------------------------------------------------------------\n');
        
        %only keep valid directory names
        dirs = dir([base_path dataset '/']);
        videos = {dirs.name};
        videos(strcmp('.', videos) | strcmp('..', videos) | ...
            strcmp('anno', videos) | ~[dirs.isdir]) = [];
        
        %the 'Jogging' sequence has 2 targets, create one entry for each.
        %we could make this more general if multiple targets per video
        %becomes a common occurence.
        if (strcmp(dataset, 'OTB100'))
            videos(strcmpi('Jogging', videos)) = [];
            videos(end+1:end+2) = {'Jogging-1', 'Jogging-2'};
            videos(strcmpi('Human4', videos)) = [];
            videos(end+1) = {'Human4.2'};
            videos(strcmpi('Skating2', videos)) = [];
            videos(end+1:end+2) = {'Skating2.1', 'Skating2.2'};
        end
        
        if (strcmp(dataset, 'OTB50'))
            videos(strcmpi('Jogging', videos)) = [];
            videos(end+1:end+2) = {'Jogging-1', 'Jogging-2'};
        end
        
        all_precisions = zeros(numel(videos),1);
        all_overlaps = zeros(numel(videos),1);
        all_fps = zeros(numel(videos),1);
        all_fn = zeros(numel(videos),1);
        
        % 		if ~exist('matlabpool', 'file'),
        %no parallel toolbox, use a simple 'for' to iterate
        for k = 1:numel(videos),
            [all_precisions(k), all_overlaps(k), all_fps(k), all_fn(k)] = run_tracker(dataset, videos{k}, kernel_type, feature_type, show_visualization, show_plots,...
                padding,lambda, output_sigma_factor,interp_factor,...
                kernel_sigma,cell_size,hog_orientations,mu,maxitr,mu_inc,step_sc);
        end
        %compute average precision, overlap, and FPS
        mean_precision_otb = mean(all_precisions);
        mean_overlap_otb = mean(all_overlaps);
        mean_fps_otb = mean(all_fps);
        mean_precision = sum(all_precisions.*all_fn)/sum(all_fn);
        mean_overlap = sum(all_overlaps.*all_fn)/sum(all_fn);
        mean_fps = sum(all_fps.*all_fn)/sum(all_fn);
        
        %Write results to file
        result_filename = fullfile(sprintf(['%.3f_' '%.3f_' '%.3g' 'fps_' dataset '_' tracker_name '_' kernel_type '_' feature_type '_padding%.3g' '_lambda%.1e' '_interp%.2g' 'cell%u' 'hog%u' 'Mu%u' 'MaxItr%u' 'IncMu%u' 'numSC%u'  'ScaleStep%u' '.txt'],...
            mean_precision_otb,mean_overlap_otb,mean_fps_otb,padding,lambda,interp_factor,cell_size,hog_orientations,mu,maxitr,mu_inc,step_sc));
        fout = fopen([result_path result_filename],'w');
        fprintf(fout,'----------------------------------------------------------------------------------\n');
        fprintf(fout,[' Dataset ......................... ', dataset '\n']);
        fprintf(fout,[' Kernel .......................... ', kernel_type '\n']);
        fprintf(fout,[' Feature ......................... ', feature_type  '\n']);
        fprintf(fout,' Padding ......................... %.3g\n', padding);
        fprintf(fout,' Lambda .......................... %.1e\n', lambda);
        fprintf(fout,' Interpolation Factor............. %.2g\n', interp_factor);
        fprintf(fout,' Output Sigma Factor ............. %.2g\n', output_sigma_factor);
        fprintf(fout,' Kernel Sigma .................... %.2g\n', kernel_sigma);
        fprintf(fout,' Cell Size ....................... %u\n', cell_size);
        fprintf(fout,' Hog Orientations ................ %u\n', hog_orientations);
        fprintf(fout,' Mu ................ %u\n', mu); 
        fprintf(fout,' Max Iterations ................ %u\n', maxitr);
        fprintf(fout,' Increase in Mu ................ %u\n', mu_inc);  
        fprintf(fout,' Step Scale ................ %u\n', step_sc);  
        fprintf(fout,'----------------------------------------------------------------------------------\n');
        
        for k = 1:numel(videos)
            fprintf(fout,' %12s - Precision (20px): %.3f, Overlap (0.5): %.3f, FPS: %.4g, FN: %u\n', videos{k}, all_precisions(k), all_overlaps(k), all_fps(k), all_fn(k));
        end
        fprintf(fout,'----------------------------------------------------------------------------------\n');
        fprintf(fout,' Average (Weighted) - Precision(%upx): %.3f, Overlap(%u%%): %.3f, FPS: %.4g \n', thresholdPrecision, mean_precision, thresholdOverlap*100, mean_overlap, mean_fps);
        fprintf(fout,' Average (OTB)      - Precision(%upx): %.3f, Overlap(%u%%): %.3f, FPS: %.4g \n', thresholdPrecision, mean_precision_otb, thresholdOverlap*100, mean_overlap_otb, mean_fps_otb);
        fclose(fout);
        
        %Print to screen
        fprintf('----------------------------------------------------------------------------------\n');
        fprintf(' Average (Weighted) - Precision(%upx): %.3f, Overlap(%u%%): %.3f, FPS: %.4g \n', thresholdPrecision, mean_precision, thresholdOverlap*100, mean_overlap, mean_fps);
        fprintf(' Average (OTB)      - Precision(%upx): %.3f, Overlap(%u%%): %.3f, FPS: %.4g \n', thresholdPrecision, mean_precision_otb, thresholdOverlap*100, mean_overlap_otb, mean_fps_otb);
        
        if (nargout > 0),
            precision = mean_precision_otb;
        end
        
        
    otherwise
        %we were given the name of a single video to process.
        
        %get image file names, initial state, and ground truth for evaluation
        [img_files, pos, target_sz, ground_truth, video_path] = load_video_info(base_path, dataset,  video);
        [positions, rects, time] = tracker(video_path, img_files, pos, target_sz, ...
            padding, kernel, lambda, output_sigma_factor, interp_factor, ...
            cell_size, features, show_visualization,mu,maxitr,mu_inc,step_sc);
        
        %calculate precision, overlap and fps
        if (~exist ('rects','var'))
            target_sz_vec = repmat(fliplr(target_sz),[length(positions) 1]);
            rects = [fliplr(positions) - target_sz_vec/2, target_sz_vec];
        end
        
        [aveErrCoverage, aveErrCenter,errCoverage, errCenter] = calcSeqErrRobust(rects, ground_truth);
        overlap = sum(errCoverage > thresholdOverlap)/length(errCoverage);
        precision = sum(errCenter <= thresholdPrecision)/length(errCenter);
        fn = numel(img_files);
        fps = fn / time;
        fprintf('%12s - Precision(%upx): %.3f, Overlap(%u%%): %.3f, FPS: %.4g,\tFN: %u \n', video, thresholdPrecision, precision, thresholdOverlap*100, overlap, fps, fn)
end
end