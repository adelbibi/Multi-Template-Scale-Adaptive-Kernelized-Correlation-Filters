function grid_search (task_id)

dataset = 'OTB100'; %OTB100, UAV
kernel_type = 'gaussian';
feature_type = 'hog';

padding = [1.5 1.7 1.8 2];
lambda = [1e-2 1e-3 1e-4 1e-5];
output_sigma_factor = [0.1];
interp_factor = [0.01 0.02 0.1];
kernel_sigma = [0.5];
cell_size = [4];
hog_orientations = [9];
mu = [1e-4 1e-5 1e-6];
maxitr = [30];
mu_inc = [1 2];
step_sc = [0.01 0.02 0.03 0.2 0.1];

parameters = combvec(padding,lambda,output_sigma_factor,interp_factor,kernel_sigma, cell_size, hog_orientations,mu,maxitr,mu_inc,step_sc);




tracker_attempts = 2;
for tracker_counter = 1:tracker_attempts
    try
        run_tracker(dataset,'all', kernel_type,feature_type,0, 0,...
            parameters(1,task_id+1),parameters(2,task_id+1),parameters(3,task_id+1),parameters(4,task_id+1),...
            parameters(5,task_id+1),parameters(6,task_id+1),parameters(7,task_id+1),parameters(8,task_id+1),parameters(9,task_id+1),parameters(10,task_id+1),parameters(11,task_id+1));
        break;
        
     catch
        if (tracker_counter < tracker_attempts)
            fprintf('Tracker failed! (Attempt [%u]).\n', tracker_counter);
		else
			fprintf('Tracker failed! (Last Attempt [%u]).\n', tracker_counter);
        			run_tracker(dataset,'all', kernel_type,feature_type,0, 0,...
            				parameters(1,task_id+1),parameters(2,task_id+1),parameters(3,task_id+1),parameters(4,task_id+1),...
            				parameters(5,task_id+1),parameters(6,task_id+1),parameters(7,task_id+1),parameters(8,task_id+1),parameters(9,task_id+1),parameters(10,task_id+1),parameters(11,task_id+1));
        end
    end
end



end
