function [lm_mv,lm_mp,lm_fwhm,el_t_LM]=...
    LM_fit2(samps,samps_n,yk,yk_n,xk,sg,step,inter,num_p,c,n_valm,sp)


%s_path  = '\\tsclient\C\Users\jil\Desktop\Kaivos\Thesis\Essay_and_Img\';
s_path  = 'C:\Users\jil\Desktop\Kaivos\Thesis\Essay_and_Img\';
sf = num2str(num_p/(step*inter));
time = (inter/num_p);

% Matrices for storing parameters
lm_mv = zeros(sg,1); lm_mp = lm_mv; lm_fwhm = lm_mv; el_time_LM = lm_mv;

% Find points in the sampled values that are higher than a threshold
% for orig gen pulse.
[slocsr,slocsc] = find(samps >= 10^-3);

pos = 1;
for i=1:sg
    nofs =  sum(slocsc(:)==i);      % number of sampled datapoints chosen from grid (this will decrease as sampling frequency decreases)
    % devides the plotting range into intervals with the step size
    %% Parameters for plotting original generated pulse
    %%{
    %This part becomes usefull when plotting
    gpp1 = find(yk == samps(slocsr(pos,1),i));           % gpp1 = generated pulse starting position
    if length(gpp1) >1
        gpp1 = gpp1(1);
    end
 
    gpp2 =find(yk == samps(slocsr(pos+nofs-1,1),i));    % generated pulse ending position
    if length(gpp2) >1
        gpp2 = gpp2(end);
    end
    slocs =time*((slocsr(pos:pos+nofs-1)-1)*step+sp(i));
    %%}
    %% Levenberg Marquardt Fitting Procedure
    %'xx' is the x-data the model uses to make the fit
    %'yy' is the y-data
    %{
    [~,mp] =max(samps(:,i));
    start =mp-50;
    stop = mp+50;
    xx=linspace(time*(start-1)*step+sp(i),time*(stop-1)*step+sp(i),(stop-start));
    yy=samps_n(start:stop-1,i);
    size(xx)
    size(yy)
    %}
    
    xx= linspace(time*((slocsr(pos,1)-1)*step+sp(i)),...
    time*((slocsr(pos+nofs-1,1)-1)*step+sp(i)),nofs);
    xx=xx.';
    thisisxx = size(xx)
    yy = samps_n(slocsr(pos,1):slocsr(pos+nofs-1,1),i)
    
    %yy = yy-n_valm;  %% THIS LINE SHOULD BE UNCOMMENTED
    
    % Guess function is the initial guessing function used by the LM-alg
    % the parameters of which the function guesses. A gaussian function 
    % is used for guessing.
    % Here x represents the gaussian parameters that are changed
    % The initial guesses for 'x' are given in 'x0'
    %if step < 31;
        % Iitial Guess for amplitude
        A_guess = max(yy);
        pos_guess = find(yy == A_guess);
        %Initial Guess for sigma 
        hm = A_guess/2;
        idx1 = find(yy>hm,1) +[-1 0];
        idx2 = find(yy>hm,1,'last') + [0 1];
        x1 = interp1(yy(idx1),xx(idx1),hm);
        x2 = interp1(yy(idx2),xx(idx2),hm);
        fwhm= x2-x1;
        sig_guess = fwhm/(2*sqrt(2*log(2)));
        x0 =[A_guess, xx(pos_guess,1), sig_guess];
    %elseif step >31 && step <81   
    %    x0 =[600, 150, 0.45];
    %end
    
    guess_fun = @(x) x(1).*exp((-(xx-x(2)).^2)./(2.*x(3).^2));
    obj=@(x) yy-guess_fun(x);
    
    
    %Lower and upper bounds for the initial guesses varied during the 
    %fitting process
    %lb = []; ub=[];
    lb = [550, 140, 0.3]; ub=[630,160,0.6];
    
    %options for lsqcurvefit telling which algorithm it should use
    %opt.Display ='iter';
    %opt.title = 'fitting Gaussian curve';
    opt.Jacobian = 'romberg';
    
    tic
    x = LevenbergMarquardt(obj,x0,lb,ub,opt);
    el_time_LM(i,1) = toc;
    
    lm_mv(i,1)= x(1);
    lm_mp(i,1)= x(2);
    lm_fwhm(i,1) = 2*sqrt(2*log(2))*x(3);
    
    % Gerate plot data from LM-retreived parameters
    xlm = xk(gpp1:gpp2);
    ylm = x(1).*exp((-(xlm-x(2)).^2)/(2*x(3).^2));
    
    
    %Plotting
    %%{
    % Colors for different plotting
    clrs = {[0.1, 0.8, 0.7],[1, 69/255, 0],[0.2, 0.2, 0.8],[0.6, 0.8, 0.1 ],[0.8, 0.1, 0.4]};
    if c ==5 && i==1
        figure(940+i);
        %%%%%%%%%%%%%%%%%%%%%%%%%% plotting origianl pulse %%%%%%%%%%%%%%%%
        %gx = linspace(xx(1),xx(length(xx)),length(yk(gpp1:gpp2)));
        %yx =  yk(gpp1:gpp2);
        %plot(gx,yx,'m','LineWidth',3)
        plot(xk(gpp1:gpp2),yk(gpp1:gpp2),'Color',clrs{1},'LineWidth',2)                            % plotting the original pulse
        hold on
        plot(xk(gpp1:gpp2),yk_n(gpp1:gpp2),'Color',clrs{2},'LineWidth',2)      % plotting noisy orig pulse
        %%%%%%%%%%%%%%%%%%%%%%%%% Plotting spline fitted data %%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%% and the sampled data-points %%%%%%%%%%%%%%%%%
        % plot for orig.gen pulse
        %plot(slocs(:,1),samps(slocsr(pos,1):slocsr(pos+nofs-1,1),i),...
        %    'o',xx,yy,'g','LineWidth',3);
        % plot for noisy gen pulse
        plot(slocs(:,1),samps_n(slocsr(pos,1):slocsr(pos+nofs-1,1),i),'o',...
            'MarkerFaceColor',clrs{5},'MarkerEdgeColor',clrs{5})
        plot(xlm,ylm,'Color',clrs{3},'LineWidth',2);
        %grid on
        axis([xk(gpp1) xk(gpp2) 0 inf])
        title(['LM fit, ' sf 'GHz'],'FontSize',20)
        xlabel('Time (ns)','FontSize',18)
        ylabel('Relative Intensity','FontSize',18)
        set(gca,'FontSize',16)
       
            % Saving Figure
            fig = gcf;
            fig.PaperPosition = [1 1 20 20];
            print([s_path 'LM_pulse_' sf 'GHz'],'-dpng', '-r0')
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    %%}
    
    
    pos=pos+nofs;
end 
el_t_LM =mean(el_time_LM);
end