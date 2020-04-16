%% Function to fit all channel recording versions of one Waveform using
%% the gaussian parametrization method
%FRD_Gauss: fit real data with Gaussian Parametrization Method
% See notes page 7.

% The flag parameter can and should be removed after testing the function
% is done

% This function (as the Levenberg-Marquardt fitting function), uses the 
% Gaussain equation to parametrize and plot the waveform. In case the 
% measured data contains values that are equal to zero some values produced
% by the parametrization alagorithm will be Inf.s or NaN.s. This is due to 
% division by zero occuring. Thus, in the function below there are some
% implementations used to remove these zeros marked by a '*' symbol

% Furthermore some the measured data values may be negative which leads to
% imaginary numbers being produced by the algorithm. Those numbers are
% removed and processed correctly in sections indicated by the symbol '**'

function [amp_err, pos_err, WF_amp_M_err, WF_amp_std, peak_pos, peak_val]= ...
    FRD_Gauss(dat,n_dat,x,MinLoc, flag,fig_s_path)
%% Parameters
% Arrays to store point group (PG) time and position values
pgp = zeros(3,3);
% Array for holding peak position and peak values for PG2 data
peak_pos = zeros(1,4); peak_val = peak_pos;
% Number of Channels | Number of PG.s | N of points in fit | length of data
NoCh=4;
NoPG=3;
% Check if number of data points is odd
if mod(n_dat,2) ~= 0
   cen = ceil(n_dat/2);
   
   quart = floor(cen/2);
   pgp(1,:) = cen-1:cen+1; %[1,cen-1:cen+1]
   pgp(2,:) = [quart,cen,cen+quart];
   pgp(3,:) = [1,cen,n_dat];
   
   
   % IF THE PULSES END UP LOOKING WIERD IN THE TIME DOMAIN THE ERROR WILL
   % BE FOUND HERE (switch to the other)
   % Spatio-temporal pos values of waveform data points for different PG.s 
   PG_t(1,:) = dat(1,pgp(1,:));
   PG_t(2,:) = dat(1,pgp(2,:));
   PG_t(3,:) = dat(1,pgp(3,:));
   
   %PG_t(1,:) = dat(pgp(1,:),1);
   %PG_t(2,:) = dat(pgp(2,:),1);
   %PG_t(3,:) = dat(pgp(3,:),1);
   
   % Vectors for storing error variables for each pulse (Can be used as error matrices)
   amp_err = zeros(1, NoCh*NoPG); pos_err = amp_err;
   
   % Vectors storing the mean error for all differnce between actual
   % datapoint values and fitted data point values in the same spatio-
   % temporal location.
   WF_amp_M_err = amp_err; WF_amp_std= amp_err;
      
   % track position (see notes on page 8.)
   trk = 1;
    
   pgp_temp = zeros(1,3);
   %Perform Gaussian fit for all channels and point groups.
   % Go through each channel one pointgroup at a time 
   for j=1:NoPG 
       for i=1:NoCh
           datt = dat(1,:);
           % Relevant raw amplitude data
           data = dat(i+1,:);
           % Extracted point group amplitude values
           a_vals = data(pgp(j,:));
           
           %% Check that the center is found correctly
           % if not find correct center value and re-establish pgp
           %(points in point-group)
           check1 = 1;
           maximum = max(data);
           
           if data(cen) ~= maximum                
                %check that last value is not the largest
                % This could be probably fixed by modifying the primary for
                % loop in process for fitting where the peak values are
                % extracted. (the condition below would not be needed, but this would reduce unnecceisarily the number of data points availble for some pulses)
                if data(n_dat) == maximum || data(n_dat)> data(n_dat-1)
                    check1 = 0;
                    pgpm = pgp;
                    [~,posi] = find(data== max(data(1:n_dat-1)));
                elseif data(1) == maximum || data(1)> data(2)    
                    check1 = 0;
                    pgpm = pgp;
                    [~,posi]= find(data == max(data(2:n_dat)));
                else
                    check1 = 0;
                    pgpm = pgp;
                    [~,posi]= find(data == max(data));  
                end
                quart = floor(cen/2);
                % condition for checking that all point groups can be
                % constructed
                %if posi+quart > || posi-quart < 1
                %
                %end   
                
                
                pgp(1,:) = posi-1:posi+1; 
                pgp(2,:) = [quart,cen,cen+quart];
                pgp(3,:) = [1,cen,n_dat];
                
                PG_t(1,:) = datt(pgp(1,:));
                PG_t(2,:) = datt(pgp(2,:));
                PG_t(3,:) = datt(pgp(3,:));
                a_vals = data(pgp(j,:));
           end
           
           %% check if any a_val is zero if so update fit-parameters*
           check = 1;
           if any(a_vals == 0) && j~=1  
              check =0;
              pgp_temp = pgp(j,:);
              pgp(j,1) = pgp(j,1)+1;
              pgp(j,3) = pgp(j,3)-1;
              a_vals = data(pgp(j,:));
              PG_t(j,:) = datt(pgp(j,:));
              n_dat2 = n_dat-2;
           end
           
           %% Check if any a_val is negative if so update fit-parameters**
           check2 = 1;           
           if any(a_vals < 0)
                check2 = 0;
                pgp_temp2 = pgp(j,:);
                %positions of the negative values in (dat(i+1,:))
                neg_pos = find(data < 0);
                % Some other value than the first or the last is negative
                if length(neg_pos) >= 2
                   %postition of center value
                   cntr = ceil(n_dat/2);
                   % Minimum 
                   mv = min(abs(cntr-neg_pos));
                   %Update amplitude and time values
                   pgp(j,:) = [cen-(mv-1) cen cen+(mv-1)];
                   a_vals = data(pgp(j,:));   
                   PG_t(j,:) = datt(pgp(j,:));
                   % update length of data array
                   n_dat2 = 2*(mv-1)+1;

                else
                    mv =cen-1;
                    pgp(j,1) = pgp(j,1)+1;
                    pgp(j,3) = pgp(j,3)-1;
                    a_vals = data(pgp(j,:));
                    PG_t(j,:) = datt(pgp(j,:));
                    
                    % update length of data array for this specific calc
                    n_dat2 = n_dat-2;                   
                end
           end
           
           %% Produce the fit
           al_vals = log(a_vals);
           
           % spatio temporal positions for different point groups
           x1 = PG_t(j,1);
           x2 = PG_t(j,2);
           x3 = PG_t(j,3);
           
           % Amplitude values
           lny1 = al_vals(1,1);
           y2 = a_vals(1,2); 
           lny2 = al_vals(1,2);
           lny3 = al_vals(1,3);
           
           % Parametrization            
           x0 = 0.5*((x1^2*(lny3-lny2)+x2^2*(lny1-lny3)+x3^2*(lny2-lny1))...
            /(x1*(lny3-lny2)+x2*(lny1-lny3)+x3*(lny2-lny1))); 
           sigma = sqrt(0.5 * ((x2-x0)^2 - (x1-x0)^2) / (lny1-lny2) );
           y0 = y2 * exp((x2-x0)^2 / (2*sigma^2));  
           %fwhm = 2*sqrt(2*log(2))*sigma;
           %vol = y0*sqrt(2*pi)*sigma;
           
           %% Error calculation
           % Error between measured peak-amp and fitted peak-amp
           %amp_err(1,trk) = abs(y0 - y2);
           % percentual error for amplitude
           if y2 > y0
               amp_err(1,trk) = y0/y2 *100;
           elseif y2 < y0
               amp_err(1,trk) = abs(1-y0/y2);
           else
               amp_err(1,trk) = 0;
           end
           
           % Error between measured peak-pos and fitted peak-pos
           pos_err(1,trk) = abs(x0 - x2);
           
           % Store fitted positional value for each channel from pg2
           if j ==2
                peak_pos(1,i) = x0;
                peak_val(1,i) = y0;
           end
           
           % Fit parametrized data using Gaussian equation and store mean
           % errors
           pg_fit = y0.*exp((-(x-x0).^2)./(2*sigma^2));
           
           % This if condition is here to avoid divisions by zero*
           if check == 0%%any(dat(i+1,:) == 0) %&& no dat values are negative 
               % update Amplitude
               dat2 = data;
               dat2(1) = [];
               dat2(n_dat-1) = [];
               % upate Location (x-axis values)
               MinLoc2 = MinLoc;
               MinLoc2(1) = [];
               MinLoc2(n_dat-1) = [];
               
               err = zeros(1,length(dat2));
          
               for g=1:length(dat2)
                  if dat2(g) > pg_fit(1,MinLoc2(g))
                       err(1,g) = 100*(1-pg_fit(1,MinLoc2(g))/dat2(g));
                  elseif dat2(g) < pg_fit(1,MinLoc2(g))    
                       err(1,g) = abs(1-pg_fit(1,MinLoc2(g))/dat2(g))*100;
                  else
                       err(1,g) = 0;  
                  end    
               end
               WF_amp_M_err(1,trk) = mean(err);
               WF_amp_std(1,trk) =sqrt(sum((err-WF_amp_M_err(1,trk)).^2)/(n_dat2-1));
           
           % Error processing for data containg negative values
           elseif check2 == 0 %any(dat(i+1,:)<0)
               % update location 
               MinLoc3 = MinLoc;
               MinLoc3(1:cen-mv)= [];
               MinLoc3(cen+mv:n_dat-(cen-mv))=[];
               % updat amplitude 
               dat_neg_A = data;
               dat_neg_A(1:cen-mv) = [];
               dat_neg_A(cen+mv:n_dat-(cen-mv)) = [];
               
               err = zeros(1,length(dat_neg_A));
               for g=1:length(dat_neg_A(1,:))
                 if dat_neg_A(g) > pg_fit(1,MinLoc3(1,g))
                     err(1,g) = 100*(1-pg_fit(1,MinLoc3(1,g))/dat_neg_A(g));
             
                 elseif dat_neg_A(g) < pg_fit(1,MinLoc3(1,g))
                     err(1,g) = abs(1-pg_fit(1,MinLoc3(1,g))/dat_neg_A(g))*100;
                 
                 else
                     err(1,g) = 0;
                 end    
               end
               WF_amp_M_err(1,trk) = mean(err);
               WF_amp_std(1,trk) =sqrt(sum((err-WF_amp_M_err(1,trk)).^2)/(n_dat2-1));
           
           % Normal (No zeros or negative numbers in the data)    
           else
              err = zeros(1,length(datt));
              for g=1:length(datt)
                 if data(g) > pg_fit(1,MinLoc(1,g))
                     err(1,g) = 100*(1-pg_fit(1,MinLoc(1,g))/data(g));
                 elseif data(g) < pg_fit(1,MinLoc(1,g))
                     err(1,g) = abs(1-pg_fit(1,MinLoc(1,g))/data(g))*100;
                 else
                     err(1,g) = 0;
                 end    
              end
    
              WF_amp_M_err(1,trk) = mean(err);
              WF_amp_std(1,trk) =sqrt(sum((err-WF_amp_M_err(1,trk)).^2)/(n_dat-1));   
           end    

           trk = trk +1;
           
           %% Some plotting to make sure things are going as planned           
           %%{
           clrs = {[0.1, 0.8, 0.7],[1, 69/255, 0],[0.2, 0.2, 0.8],[0.6, 0.8, 0.1 ],[0.8, 0.1, 0.4]};
           % One figure for each channel
            if flag == 15 && i ==1 && j==1 %&& any(WF_amp_M_err > 50)
                %dat(i+1,:)
                %pg_fit(1,MinLoc(1,:));
                %WF_amp_M_err
                %figure(i+i-1);
                figure(flag)
                plot(datt, data, 'o', 'MarkerEdgeColor',clrs{1}, 'MarkerFaceColor', clrs{1},'HandleVisibility','off','MarkerSize',10)
                hold on
                plot(x, pg_fit(1,:),'Color', clrs{j},'LineWidth',1.5);
                legend('GP (PG1)','GP (PG2)','GP (PG3)','CS','PLSQ','LM');  
                %plot(x(1,MinLoc(1,:)),pg_fit(1,MinLoc(1,:)),'o', 'MarkerEdgeColor','g', 'MarkerFaceColor', 'g')
               
                %{
                title('Gaussian');
                xlabel('time (ns)')
                ylabel('Amplitude')
                grid on
                %}
            end
           %%}
           %{
           if  (WF_amp_M_err > 1000)
                %dat(i+1,:)
                %pg_fit(1,MinLoc(1,:));
                WF_amp_M_err
                %figure(i+i-1);
                figure(10000+i)
                plot(x, pg_fit(1,:));
                hold on
                plot(x(1,MinLoc(1,:)),pg_fit(1,MinLoc(1,:)),'o', 'MarkerEdgeColor','g', 'MarkerFaceColor', 'g')
                plot(datt, data, 'o', 'MarkerEdgeColor','m', 'MarkerFaceColor', 'm')
                title('Gaussian');
                xlabel('time (ns)')
                ylabel('Amplitude')
                grid on
            end
            %}
            if check1 == 0
                pgp = pgpm;
            end
                
            if check == 0
                pgp(j,:) = pgp_temp;
            end
            
            if check2 == 0
                pgp(j,:) = pgp_temp2;
            end    
       end
   end    

else   
   disp('Number of data points in pulse is not odd')
   disp('See parameter n_sub in function process_for_fitting.m and notes')
   disp('on page 9')
end    

%end