%% Function for performing polynomial fitting

function [amp_err, pos_err, WF_amp_M_err, WF_amp_std] = ...
           FRD_Poly(dat,n_dat,x,NoCh,MinLoc,flag,fig_s_path)
       
   % Vectors for storing error variables for each pulse (Can be used as error matrices)
   amp_err = zeros(1, NoCh); pos_err = amp_err;
   % Vectors storing the mean error for all differnce between actual
   % datapoint values and fitted data point values in the same spatio-
   % temporal location.
   WF_amp_M_err = amp_err; WF_amp_std= amp_err; 
   % Order of polynomial
   n = length(dat(1,:))-2; %n=7; % Try changing this value    
   for i=1:NoCh
       [ply,~,mu] = polyfit(dat(1,:), dat(i+1,:),n);
       plyy = polyval(ply,x,[],mu);
    
       % Find maximum value from fit
       % Calculate Amplitude Error
       amp_err(1,i) = abs(max(plyy)-max(dat(i+1,:)));
       %Calulate percentage
       amp_err(1,i) = (amp_err(1,i)/max(dat(i+1,:)))*100;
       
       
       % Calculate postion Error
       [~,p] = max(plyy);
       pos_err(1,i) = abs(x(1,p)-dat(1,floor(n_dat/2)));
       
       % Calculate WF error
       %err = abs(plyy(MinLoc)-dat(i+1,:));
       %err = abs(1-abs(dat(i+1,:)./plyy(1,MinLoc(1,:))));
       
       % THIS PART IS NEW
       err = zeros(1,length(dat(i+1,:)));

       for g=1:length(dat(i+1,:))
           if dat(i+1,g) > plyy(1,MinLoc(1,g))
               err(1,g) = 100*(1-plyy(1,MinLoc(1,g))/dat(i+1,g));
             
           elseif dat(i+1,g) < plyy(1,MinLoc(1,g))
               err(1,g) = abs(1-plyy(1,MinLoc(1,g))/dat(i+1,g))*100;
           else
               err(1,g) = 0;
           end
       end
       %%%%%%%%%%%%%%%%%%%
       
       %err = abs(1-abs(plyy(MinLoc)./dat(i+1,:)));
       WF_amp_M_err(1,i) = mean(err);
       WF_amp_std(1,i) = sqrt(sum((err-WF_amp_M_err(1,i)).^2) /(n_dat-1));
       
      % Plotting functions to check that everything is ok. Remove later
      %%{
      % One figure for each channel
      clrs = {[0.1, 0.8, 0.7],[1, 69/255, 0],[0.2, 0.2, 0.8],[0.6, 0.8, 0.1 ],[0.8, 0.1, 0.4]};
      if (flag == 50 && i==1) %|| (flag == 3 && i == 2) ...
              %|| (flag == 7  && i == 3) || (flag == 9 && i == 4) 
          figure(flag)
          %plot(dat(1,:), dat(i+1,:), 'o', 'MarkerEdgeColor','m', 'MarkerFaceColor', 'm')
          hold on
          plot(x, plyy,'Color', clrs{5},'LineWidth',1.5)
          legend('GP (PG1)','GP (PG2)','GP (PG3)','CS','PLSQ','LM'); 
          %{
          title('Poly Fit');
          xlabel('time (ns)')
          ylabel('Amplitude')
          grid on
          %}
      end
      %%}
   end     
       
       
       
       
end       