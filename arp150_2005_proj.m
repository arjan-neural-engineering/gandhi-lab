%% Arjan Singh Puniani 
%% BIOE 2005 - Computer Proj 1
% we will be calculating the dielectric constant from experimentally
% measured reflection coefficients

%% method: newton
gammas = [-0.1 -0.5 -0.9];
alphas = [0.1 0.01 0.001 0.0001];
iterations = 100000000;
for i=1:numel(gammas)
    for j=1:numel(alphas)
        for theta_i=pi/9:pi/9:(4*pi/9)
              epsilon = newt(gammas(i),iterations,theta_i,alphas(j));
              disp("for gamma " + gammas(i) + " with alpha " + alphas(j))
              disp("our \epsilon at " + (180/pi)*theta_i + " is " + epsilon)
        end 
    end
end


%% method: naive guessing
% gamma_meas = -0.1;
% alpha = 0.01; % tolerance, from 0.1 to 0.0001
% for theta_i=pi/9:pi/9:(4*pi/9) % try each incident theta (in radians)
%    epsilon_r = 1.000001; % each incident theta gets its own candidate epsilon
%     for trial=1:1e11 % 100 billion trials
%        theta_t = asin(sin(theta_i)/(sqrt(epsilon_r))); % answer in radians
%        gamma_exact = (cos(theta_t)-cos(theta_i)*sqrt(epsilon_r))/(cos(theta_t)+cos(theta_i)*sqrt(epsilon_r));
%        %% check if your analytical (exact) solution (gamma_exact) is within the tolerance 
%        if (abs(gamma_meas-gamma_exact)<alpha)
%            disp("found it on trial #" + trial + " with \Theta_i = " + (180/pi)*theta_i)
%            disp( "\Theta_t = " + (180/pi)*theta_t + " has epsilon of " + epsilon_r)
%           break  
%        end
%        %% update epsilon
%        epsilon_r = epsilon_r + 0.00000001;
%     end
% end



%     theta_t1 = asin(sin(theta_i)/(sqrt(x1)));
%     theta_t2 = asin(sin(theta_i)/(sqrt(x2)));
  

   
    % no need to code y since we know its 0
%     y1 = gamma_meas - (cos(theta_t1)-cos(theta_i)*sqrt(x1))/(cos(theta_t1)+cos(theta_i)*sqrt(x1));
%     y2 = gamma_meas - (cos(theta_t2)-cos(theta_i)*sqrt(x2))/(cos(theta_t2)+cos(theta_i)*sqrt(x2));
    
    % tolerance, from 0.1 to 0.0001
%     if (abs(y1)<alpha)
%         disp('we done')
%         disp("our \epsilon is " + x1)
%     elseif (abs(y2)<alpha)
%         disp('we done')
%         disp("our \epsilon is " + x2)
%     else
%         x3 = abs(((x2-x1)*(-y1)/(y2-y1))+x1);
%        % theta_t = asin(sin(theta_i)/(sqrt(x3)));      % answer in radians
%     end



% while (abs(gamma_meas-gamma_exact)>alpha)
%     epsilon_r=epsilon_r+0.00000001    
% end
% theta_i = 0;
% epsilon_r = 1.01;
% theta_t = asin(sin(theta_i)/(sqrt(epsilon_r)))
       %theta_t = asin(sin(theta_i)/(sqrt(epsilon_r)));
%          theta_t = (360*asin(sin(theta_i)/(sqrt(epsilon_r)))/(2*pi))+90;
%gamma_exact = (cos(theta_t*180/pi)-cos(theta_i*180/pi)*sqrt(epsilon_r))/(cos(theta_t*180/pi)+cos(theta_i*180/pi)*sqrt(epsilon_r));