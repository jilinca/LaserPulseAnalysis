%% LEVENBERG-MARQUARDT FIT

% This function (as the gaussian parametrization function), uses the 
% Gaussain equation to parametrize the and plot the waveform. In case the 
% measured data contains values that are equal to zero some values produced
% by the parametrization alagorithm will be Inf.s or NaN.s. This is due to 
% division by zero occuring. Thus in the function below there are some
% implementations used to remove these zeros marked by an * symbol

% **FUNCTIONALITY FOR GUESSING SIGMA
% This functionality does not work properly due to the many different
% variations in the data that need to be taken into accoutnt so that 
% each FWHM for each waveform can be calculated correctly 
% (It is not strictly needed but would aid the levenberg marquardt
% to do processing faster due to a better initial guess given by it)

function [amp_err, pos_err, WF_amp_M_err, WF_amp_std] = ...
           FRD_LM(dat,n_dat,x,NoCh,MinLoc,flag,fig_s_path)
       
   % Vectors for storing error variables for each pulse (Can be used as error matrices)
   amp_err = zeros(1, NoCh); pos_err = amp_err;
   % Vectors storing the mean error for all differnce between actual
   % datapoint values and fitted data point values in the same spatio-
   % temporal location.
   WF_amp_M_err = amp_err; WF_amp_std= amp_err; 
   % Order of polynomial
       
   for i=1:NoCh
        % Initial Guess for amplitude
        %A_guess = dat(i+1,ceil(n_dat/2));
        A_guess = max(dat(i+1,:));
        [~,position] = find(dat(i+1,:)== A_guess);
        %Inital Guess for peak position
        pos_guess = dat(1,position);
        
        % FUNCTION FOR GUESSING SIGMA** (See note on top of function**)
        %{
        %Initial Guess for sigma 
        hm = A_guess/2; %Half Maximum
        idx1 = find(dat(i+1,:)>hm,1) +[-1 0];
        idx2 = find(dat(i+1,:)>hm,1,'last') + [0 1];
        flag
        %n_dat
        if idx1(1,1) == 0 && idx2(1,2) > n_dat 
           if dat(i+1,1) == dat(i+1,n_dat)
               fwhm = dat(1,1) - dat(1,n_dat);
           elseif dat(i+1,1) < dat(i+1,n_dat)
                disp ('kek')
                idx = find(dat(i+1,:)>dat(i+1,n_dat),1)
                x2 = dat(i+1,n_dat);
                x1 = interp1(dat(i+1,[idx-1 idx]),...
                    dat(1,[idx-1 idx]),dat(i+1,n_dat))
                fwhm = x2-x1;
           elseif dat(i+1,1) > dat(i+1,n_dat)
                disp('kok')
                idx = find(dat(i+1,:)>dat(i+1,1),1,'last');
                x1 = dat(i+1,1);
                x2 = interp1(dat(i+1,[idx idx+1]),...
                    dat(1,[idx idx+1]),dat(i+1,1));
                %x2 = interp1(dat(i+1,[n_dat-1 n_dat]),...
                %    dat(1,[n_dat-1 n_dat]),dat(i+1,1))
                fwhm = x2-x1;  
           end
        elseif idx1(1,1) == 0 && idx2(1,2) <= n_dat
           idx2 = find(dat(i+1,:)>dat(i+1,1),1,'last') + [0 1];
           x2 = interp1(dat(i+1,idx2),dat(1,idx2),dat(i+1,1));
           fwhm = x2 - dat(1,1);
        elseif idx2(1,2) > n_dat && idx1(1,1) ~= 0
           val = find(dat(i+1,:)>dat(i+1,n_dat),1)
           idx1 = find(dat(i+1,:)>dat(i+1,n_dat),1) + [-1 0];
           x1 = interp1(dat(i+1,idx1),dat(1,idx1),dat(i+1,n_dat));
           fwhm = dat(1,n_dat)-x1;
        else
            %idx1 = find(dat(i+1,:)>hm,1) +[-1 0];
            %idx2 = find(dat(i+1,:)>hm,1,'last') + [0 1];
            x1 = interp1(dat(i+1,idx1),dat(1,idx1),hm);
            x2 = interp1(dat(i+1,idx2),dat(1,idx2),hm);
            fwhm= x2-x1;
        end
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %sig_guess = fwhm/(2*sqrt(2*log(2)))
        sig_guess = 0.5;
        x0 =[A_guess, pos_guess, sig_guess];
        
        % Remove zeros from data*
        if any(dat(i+1,:)== 0)
           %Time 
           dat1 = dat(1,:);
           dat1(1) =[];
           dat1(n_dat-1) =[];
           %amplitude
           dat2 = dat(i+1,:);
           dat2(1) = [];
           dat2(n_dat-1) = [];
           guess_fun = @(z) z(1).*exp((-(dat1(1,:)'-z(2)).^2)./(2.*z(3).^2));
           obj=@(z) dat2(1,:)'-guess_fun(z);
        else    
           guess_fun = @(z) z(1).*exp((-(dat(1,:)'-z(2)).^2)./(2.*z(3).^2));
           obj=@(z) dat(i+1,:)'-guess_fun(z);
        end   
    
    
        %Lower and upper bounds for the initial guesses
        lb = [A_guess-5, pos_guess-0.5, sig_guess-0.4];
        ub = [A_guess+5, pos_guess+0.5, sig_guess+0.4];
        %options for lsqcurvefit telling which algorithm it should use
        %opt.Display ='iter';
        %opt.title = 'fitting Gaussian curve';
        opt.Jacobian = 'romberg';
    
        z = LevenbergMarquardt(obj,x0,lb,ub,opt);
    
        %amp = z(1);
        %pos = z(2);
        %lm_fwhm = 2*sqrt(2*log(2))*z(3);
        
        % Generate full-waveform from LM parameters
        WFy = z(1).*exp((-(x-z(2)).^2)/(2*z(3).^2));
       
       %% Error Calculations 
       %amp_err(1,i) = abs(max(WFy)-max(dat(i+1,:)));
       % Percentual Error
       
       
       %amp_err(1,i) = abs(1-abs(max(WFy)./max(dat(i+1,:))));
       if max(dat(i+1,:)) > max(WFy)
             amp_err(1,i) = (1-max(WFy)/max(dat(i+1,:))) *100;
       elseif max(dat(i+1,:)) < max(WFy)
            amp_err(1,i) = abs(1-max(WFy)/max(dat(i+1,:)))*100;
       else
            amp_err(1,i) = 0;
       end    
       
       
       % Calculate postion Error
       [~,p] = max(WFy);
       pos_err(1,i) =abs(x(1,p)-dat(1,ceil(n_dat/2)));
       
       
       
       % Calculate WF error
       %err = abs(WFy(MinLoc)-dat(i+1,:));
       % Take into account the case of removed zeros*
       if any(dat(i+1,:)==0)
         MinLoc2 = MinLoc;
         MinLoc2(1) = [];
         MinLoc2(n_dat-1) = [];
         
         err = zeros(1,length(dat1));
         for g=1:length(dat1)
             if dat2(g) > WFy(1,MinLoc2(g))
                  err(1,g) = 100*(1-WFy(1,MinLoc2(g))/dat2(g));
             elseif dat2(g) < WFy(1,MinLoc2(g))    
                  err(1,g) = abs(1-WFy(1,MinLoc2(g))/dat2(g))*100;
             else
                  err(1,g) = 0;  
             end    
         end
         %err = abs(1-abs(WFy(MinLoc2)./dat2(1,:)));
         WF_amp_M_err(1,i) = mean(err);
         % n_dat-3 because 2 data-points were removed
         WF_amp_std(1,i) = sqrt(sum((err-WF_amp_M_err(1,i)).^2) /(n_dat-3));    
       
       else
         %err = abs(1-abs(WFy(MinLoc)./dat(i+1,:)));
         err = zeros(1,length(dat(1,:)));
         for g=1:length(dat(1,:))
            if dat(i+1,g) > WFy(1,MinLoc(1,g))
                err(1,g) = 100*(1-WFy(1,MinLoc(1,g))/dat(i+1,g));
            elseif dat(i+1,g) < WFy(1,MinLoc(1,g))
                err(1,g) = abs(1-WFy(1,MinLoc(1,g))/dat(i+1,g))*100;
            else
                err(1,g) = 0;
            end    
         end
         WF_amp_M_err(1,i) = mean(err);
         WF_amp_std(1,i) = sqrt(sum((err-WF_amp_M_err(1,i)).^2) /(n_dat-1));
       end
       
       
       
       
      %% Plotting functions to check that everything is ok. Remove later
      %%{
      clrs = {[0.1, 0.8, 0.7],[1, 69/255, 0],[0.2, 0.2, 0.8],[0.6, 0.8, 0.1 ],[0.8, 0.1, 0.4],[0.4, 0.1, 0.2]};
      % One figure for each channel
      if (flag == 50 && i==1) %|| (flag == 3 && i == 2) ...
              %|| (flag == 7  && i == 3) || (flag == 9 && i == 4) 
          figure(flag)
          %plot(dat(1,:), dat(i+1,:), 'o', 'MarkerEdgeColor','m', 'MarkerFaceColor', 'm')
          hold on
          plot(x, WFy,'Color', clrs{6},'LineWidth',1.5)
          legend('GP (PG1)','GP (PG2)','GP (PG3)','CS','PLSQ','LM'); 
          %title('LM Fit');
          xlabel('time (ns)','FontSize',28)
          ylabel('Amplitude','FontSize',28)
          grid on
          set(gca,'FontSize',28)
          
    	  % Saving Figure
          %fig = gcf;
          %fig.PaperPosition = [1 1 20 20];
          %print([fig_s_path 'real_fit_5GHz'],'-dpdf')
      end
      %%}
   end     
   
end