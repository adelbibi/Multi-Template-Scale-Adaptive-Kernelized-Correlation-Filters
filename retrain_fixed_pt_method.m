function [alphaf_1_new,alphaf_2_new] = retrain_fixed_pt_method(xf_1,alphaf_1,xf_2,alphaf_2,y,kernel,lambda,mu,maxitr,mu_inc)

k=mu/(lambda + mu); 
c = lambda*k + mu*(k-1);
yf = fft2(y);

switch kernel.type 
    case 'gaussian',
        kf_1_1 = gaussian_correlation(xf_1, xf_1, kernel.sigma);
        kf_2_2 = gaussian_correlation(xf_2, xf_2, kernel.sigma);
        kf_2_1 = gaussian_correlation(xf_1, xf_2, kernel.sigma);
        kf_1_2 = gaussian_correlation(xf_2, xf_1, kernel.sigma);
    case 'polynomial',
        kf_1_1 = polynomial_correlation(xf_1, xf_1, kernel.poly_a, kernel.poly_b);
        kf_2_2 = polynomial_correlation(xf_2, xf_2, kernel.poly_a, kernel.poly_b);
        kf_2_1 = polynomial_correlation(xf_1, xf_2, kernel.poly_a, kernel.poly_b);
        kf_1_2 = polynomial_correlation(xf_2, xf_1, kernel.poly_a, kernel.poly_b);

    case 'linear',
        kf_1_1 = linear_correlation(xf_1, xf_1); 
        kf_2_2 = linear_correlation(xf_2, xf_2);    
        kf_2_1 = linear_correlation(xf_1, xf_2);
        kf_1_2 = linear_correlation(xf_2, xf_1);
end

i=1;
while(true)  
%Force filter of second training example to look similar to first one
%% new (Update Filter2)
num = conj(yf) - ( ((c./(kf_2_2+eps)) + k) .* conj(kf_1_2 .* alphaf_1) );
den = kf_2_2 + (lambda+mu);
alphaf_2_new = conj(num ./ den);
%% new (Update Filter1)
num = conj(yf) - ( ((c./(kf_1_1+eps)) + k) .* conj(kf_2_1 .* alphaf_2) );
den = kf_1_1 + (lambda+mu);
alphaf_1_new = conj(num ./ den);
%% Both Filter Update
alphaf_1 = alphaf_1_new;
alphaf_2 = alphaf_2_new;
%% Stopping Crieterion
data_cost(i,:) = real(sum(sum(fft2((kf_2_2 .* conj(alphaf_2_new)) - conj(yf)).^2)));
if(i>5)
    if(std(data_cost(i:-1:i-4,:)) <=10^-5 || i>maxitr)
        break;
    end
end
i=i+1;
mu = mu * mu_inc;
end
end