%% PROCESS FOR FITTING PRIMARY
% The structure darras (change name to ch data) containes all the
% processable waveforms. 
warning('off','all')
fig_s_path = 'C:\Users\jil\Desktop\Kaivos\PARAINEN\Codes\Plots\';
% ret_pos
% The first column of each array stored
% in the structure contains the waveform number. The data appearing after
% the waveform number contains the positional values of the positive slopes
% (indicating a found pulse) of the waveform in question reffered to by the 
% waveform number. The postional values can be directly used to extract the
% pulse values for the rising side of the pulse from the channel data
% (chX_d)

%function [] = process_for_fitting(darrs,ret_pos)%% Parameters
[lr, ~] = size(darrs.ad1);
l_ret_pos = length(ret_pos(:,1));
trigs2 = trigs;
NoPiF = 200;
%Indicates on which row of ret_pos the primary return pulse of a waveform
%is found.
PrimaryReturnRows = zeros(l_ret_pos,1);


%% Get row indices for the primary returns from ret_pos
% sometimes the primary return is not the first return see notes on page 9.
i =1;
x =1;

while i <= l_ret_pos
   % waveform number
   wfn = ret_pos(i,1);
   % peak values per waveform
   pvpwf = [];
   %i2 = i; 
   % Count number of returns in wavefom (to make code faster an 
   % array of the number of pulses per waveform could be returned from 
   % find_return_peak_locations.m, this would need a bit more processing 
   % to be done in Align_and_multiply.m. So do it if you have time later)
   while (wfn == ret_pos(i,1))
       %last data point
       ldp = find(ret_pos(i,:)==0,1,'first')-1;
       % Get the maximum value from each found pulse in a waveform
       pvpwf = [pvpwf darrs.ad2(x,ret_pos(i,ldp))];  %#ok<AGROW>
       
       % implement check for last index
       if i == l_ret_pos
            [~,mp] = max(pvpwf);
            cc = length(pvpwf);
            PrimaryReturnRows(i-cc+(mp),1) = 1;
            break
       end  
       i = i+1; 
   end
   
   if i == l_ret_pos
       break
   else 
        % Here i is at the index of the next waveform
        [~,mp] = max(pvpwf);
        cc = length(pvpwf);
        PrimaryReturnRows(i-cc+mp-1,1) = 1;
   end
   x = x+1; 
end

disp("YEY, GOT HERE")
%{
for i=1:l_ret_pos
    if i == 1
        PrimaryReturnRows(i,1) = 1;
    else
    % Check form ret_pos if the waveform number is the same as the previous
    % if yes, change value to zero 
    
    % A minor reflection might also occur before the so called primary
    % return, to separate this we check which return has the highes peak
    % value and process this as the primary return (see page 9 of notes)
        if ret_pos(i,1) == ret_pos(i-1,1)
            % Process reflection befor primary return
            frstZ1 = find(ret_pos(i,:)==0,1,'first');
            frstZ2 = find(ret_pos(i-1,:)==0,1,'first');
            if i == 494 || i == 495
                frstZ1
                frstZ2
            end    
            if darrs.ad1(ret_pos(i-1,frstZ2-1)) < darrs.ad1(ret_pos(i,frstZ1-1))  
                PrimaryReturnRows(i-1,1) = 0;

            else    
                PrimaryReturnRows(i,1) = 0;
            end
        end
    end     
end    
%}

%Remove all but primary return pulses
%PRIMARY RETURN (PULSE) POSITIONS = prp 
if l_ret_pos == length(PrimaryReturnRows)
    prp = ret_pos.*PrimaryReturnRows;
    prp = prp(any(prp,2),:);
    disp('I went here')
else
    error('Error')
end

l_rpp = length(prp(:,1));

% Peak amplitude and postion errors and difference between measured and
% fitted value squared.
GA_err = zeros(l_rpp,12); GP_err = GA_err;
SA_err = zeros(l_rpp,4); SP_err = SA_err;
PA_err = SA_err; PP_err = SA_err;
LA_err = SA_err; LP_err = SA_err;

% Arrays to store peak value and position found using the PG2 of GP
peak_pos = zeros(l_rpp,4); peak_val = peak_pos;
CS_pval = peak_pos; CS_tval = peak_pos;
% The arrays below store the mean and standard deviation of the differences
% calculated betwen measured data points and fitted data points for all
% data points comprising the pulse. IE it considers a pulse as a whole, and
% not just the peak amplitude.

% GAME = Gaussian Amplitude Mean Error
% GASD = Gaussian Amplitude Standard Deviation
GAME = GA_err; GASD = GA_err;
SAME = SA_err; SASD = SA_err;
PAME = SA_err; PASD = SA_err;
LAME = SA_err; LASD = SA_err;

%% Extract first returns data for first returns and Call Fitting Functions
err_msg2 = 'The number of primary returns does not match the number of waveforms';
n_sub = 3;
n_dat =zeros(l_rpp,1);

% Array to store skipped pulse row numbers in the next for loop
skpd_r =[];

if l_rpp ~= length(darrs.td1(:,1))
    length(darrs.td1(:,1))
    error(err_msg2)    
else 
   for i=1:l_rpp  
        %% Extract Primary returns
        % find at which position the position values end
        frst = find(prp(i,:)==0,1,'first');
        % NOTE! The value of n_sub might have to be tweaked depending on
        % the data set due to the oscillations behind the peaks. See notes
        % in notebook page 5. on the logic behind the next 6 lines.
        dat = zeros(5,(frst-2)+frst-n_sub);
        % N. of data points used to represent the pulse being processed
        n_dat(i,1) = 2*frst-(2+n_sub);
        dat(1,:) = darrs.td1(i,prp(i,2):prp(i,(frst-1))+frst-n_sub);
        dat(2,:) = darrs.ad1(i,prp(i,2):prp(i,(frst-1))+frst-n_sub);
        dat(3,:) = darrs.ad2(i,prp(i,2):prp(i,(frst-1))+frst-n_sub);
        dat(4,:) = darrs.ad3(i,prp(i,2):prp(i,(frst-1))+frst-n_sub);
        dat(5,:) = darrs.ad4(i,prp(i,2):prp(i,(frst-1))+frst-n_sub);
        
        % Skip fittig if data contains seven or less data points and negative numbers
        % This causes problems for the gaussian function when trying to use different pg.s
        % and complicates the error analysis, especially if there are multiple negative values
        % multiple negative values. This is not a problem in a real
        % situation where different PG.s will not be used, or the code is
        % written to choose the most suitable PG automatically.
        
        %% 
        skip = 0;
        for q=2:length(dat(:,1))
            gg = dat(q,:);
            if length(dat(q,:)) <= 7 && length(gg(gg<0)) >= 1 
               skip =1;
               skpd_r = [skpd_r;i]; %#ok<*AGROW>
            end    
        end
        
        % Check that the largest value on each channel is not the first
        % or the last value.
        
        for q=2:length(dat(:,1))
           if max(dat(q,:)) == dat(q,1) || max(dat(q,:)) == dat(q,n_dat(i,1)) % length(dat(1,:))
              skip =1;
              % Add skipped pulse row number only if it not added already
              % in the previous if clause
              if isempty(skpd_r) == 0 % if not empty
                if skpd_r(length(skpd_r)) ~= i
                    skpd_r = [skpd_r;i];
                end
              else
                 skpd_r = [skpd_r;i];  
              end
           end
        end

        if skip == 0
        % X-axis time values for performing fit
        x = linspace(dat(1,1),dat(1,n_dat(i,1)),NoPiF);
        
        % Find indices of values closest to those of the dat
        % time values in the x-axis (time) values
        MinLoc = zeros(1,n_dat(i,1));            
        for k=1:n_dat(i,1)
            Subs = abs(dat(1,k)-x);
            % MinLoc holds the locations of the dat time values in x 
            minloc = find(Subs == min(Subs));
            if length(minloc) > 1
                MinLoc(1,k) = minloc(1,1);
            else
                MinLoc(1,k) = minloc;
            end
        end
        
        
        %% Call fitting functions
        % All fitted pulses will contain the same number of points.
        % The number is the same as for the simulated pulses.
        %NoPiF= 200; %Number of points in fit
        NoCh = 4;
        % REMEBER TO REMOVE LAST PARAEMTER i FROM FROM FUNCTION CALL BELOW
        
        % Gaussian Fitting
        [GA_err(i,:), GP_err(i,:), GAME(i,:), GASD(i,:),peak_pos(i,:), peak_val(i,:)] = ... 
           FRD_Gauss(dat, n_dat(i,1), x, MinLoc,i);
        % Cubic Spline Fitting
        [SA_err(i,:), SP_err(i,:), SAME(i,:),SASD(i,:),CS_pval(i,:),CS_tval(i,:)] = ...
           FRD_Spline(dat, n_dat(i,1),x, NoCh, MinLoc,i);
        %Polynomial 7th order
        [PA_err(i,:), PP_err(i,:), PAME(i,:), PASD(i,:)] = ...
           FRD_Poly(dat, n_dat(i,1), x,NoCh, MinLoc, i);
        % Levenberg-Marquardt
        [LA_err(i,:), LP_err(i,:), LAME(i,:), LASD(i,:)] = ...
           FRD_LM(dat, n_dat(i,1), x,NoCh, MinLoc, i);

        end       
        
        
        
   end
   
   %% %%%%%%%%%%%%%%%%%%% FINAL ERROR CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%
   
   % Calc mean error: The error data is stored so that data corresponding
   % To a particular channel, point group and fitting method is stored in
   % one column. 
   % NOTE: You can concatenate all the fitting methods later to a larger
   % array so that all mean amplitude errors can be calculated at once
   % Now the processing is done for the Gaussian Parametrization method
   
   %% Amp Mean Error
   dpx = data_pos_x; dpy = data_pos_y;
   A_ERRS = [GA_err SA_err PA_err LA_err];
   %{
   %%%%%%%%%%%REMOVING UNWANTED DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
   %% Find NANs and remove them from all arrays
   re2 = 1-any(isnan(A_ERRS),2);
   A_ERRS = A_ERRS.*re2;
   A_ERRS = A_ERRS(any(A_ERRS,2),:);
   
   
    %% Remove inf.s 
   remz2 = 1-any(isinf(A_ERRS),2);
   A_ERRS = A_ERRS.*remz2;
   A_ERRS = A_ERRS(any(A_ERRS,2),:);
   
   %% Reomve all errors larger than 50% These pulses were not fitted properly
   [re_r,~] = find(A_ERRS(:,1)>=50);
   for i =1:length(re_r)
       A_ERRS(i,:) = [];
   end    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   A_FIN_ERR = mean(A_ERRS);
   % Amp standard deviation
   A_FIN_STD = sqrt(sum((A_FIN_ERR-A_ERRS).^2)/(length(A_ERRS(:,1)-1)));
   
   % The arrays below are created to speed up the computuation of the 
   % combined standard deviation
   % mean of errors for each point group and fittign method (All 4 channels
   % combined)
   
   TEMP_A = zeros(1,6);
   TEMP_A(1,1) = mean(A_FIN_ERR(1:4));
   TEMP_A(1,2) = mean(A_FIN_ERR(5:8));
   TEMP_A(1,3) = mean(A_FIN_ERR(9:12));
   TEMP_A(1,4) = mean(A_FIN_ERR(13:16));
   TEMP_A(1,5) = mean(A_FIN_ERR(17:20));
   TEMP_A(1,6) = mean(A_FIN_ERR(21:24));
   
   % See formula for combined standard deviation
   % sum of std for each point group and fitting method. (All 4 channels
   % are summed for each fitting method)
   TEMP_AS=zeros(1,6);
   TEMP_AS(1,1) = sum((A_FIN_STD(1:4)).^2);
   TEMP_AS(1,2) = sum((A_FIN_STD(5:8)).^2);
   TEMP_AS(1,3) = sum((A_FIN_STD(9:12)).^2);
   TEMP_AS(1,4) = sum((A_FIN_STD(13:16)).^2);
   TEMP_AS(1,5) = sum((A_FIN_STD(17:20)).^2);
   TEMP_AS(1,6) = sum((A_FIN_STD(21:24)).^2);
   
   
   % Difference between mean value of means (taken from indivudal fits)
   % and the mean of individual fits squared and summed
   TEMP_A_DIF= zeros(1,6);
   TEMP_A_DIF(1,1) = sum((TEMP_A(1,1)-A_FIN_ERR(1:4)).^2);
   TEMP_A_DIF(1,2) = sum((TEMP_A(1,2)-A_FIN_ERR(5:8)).^2);
   TEMP_A_DIF(1,3) = sum((TEMP_A(1,3)-A_FIN_ERR(9:12)).^2);
   TEMP_A_DIF(1,4) = sum((TEMP_A(1,4)-A_FIN_ERR(13:16)).^2);
   TEMP_A_DIF(1,5) = sum((TEMP_A(1,5)-A_FIN_ERR(17:20)).^2);
   TEMP_A_DIF(1,6) = sum((TEMP_A(1,6)-A_FIN_ERR(21:24)).^2);
   % Number of fits per channel
   len = length(A_ERRS(:,1));
   
   %Calculation of Combined standard deviation for each fitting method
   A_FIN_STD = sqrt((4*(len-1).*TEMP_AS+4*len*TEMP_A_DIF)/(4*len-1));
   
   
   % Relative standard deviation of amplitude
   A_FIN_RSTD = (A_FIN_STD./TEMP_A)*100;
   
   %% Position Mean Error
   P_ERRS = [GP_err SP_err PP_err LP_err];
   % Remove all rows with zeros

   
   P_FIN_ERR = mean(P_ERRS);
   % Position standard deviation
   P_FIN_STD = sqrt(sum((P_FIN_ERR-P_ERRS).^2)/(length(P_ERRS(:,1)-1)));
   TEMP_P = zeros(1,6);
   TEMP_P(1,1) = mean(P_FIN_ERR(1:4));
   TEMP_P(1,2) = mean(P_FIN_ERR(5:8));
   TEMP_P(1,3) = mean(P_FIN_ERR(9:12));
   TEMP_P(1,4) = mean(P_FIN_ERR(13:16));
   TEMP_P(1,5) = mean(P_FIN_ERR(17:20));
   TEMP_P(1,6) = mean(P_FIN_ERR(21:24));
   
   % See formula for combined standard deviation
   TEMP_PS=zeros(1,6);
   TEMP_PS(1,1) = sum((P_FIN_STD(1:4)).^2);
   TEMP_PS(1,2) = sum((P_FIN_STD(5:8)).^2);
   TEMP_PS(1,3) = sum((P_FIN_STD(9:12)).^2);
   TEMP_PS(1,4) = sum((P_FIN_STD(13:16)).^2);
   TEMP_PS(1,5) = sum((P_FIN_STD(17:20)).^2);
   TEMP_PS(1,6) = sum((P_FIN_STD(21:24)).^2);
   
   TEMP_P_DIF= zeros(1,6);
   TEMP_P_DIF(1,1) = sum((TEMP_P(1,1)-P_FIN_ERR(1:4)).^2);
   TEMP_P_DIF(1,2) = sum((TEMP_P(1,2)-P_FIN_ERR(5:8)).^2);
   TEMP_P_DIF(1,3) = sum((TEMP_P(1,3)-P_FIN_ERR(9:12)).^2);
   TEMP_P_DIF(1,4) = sum((TEMP_P(1,4)-P_FIN_ERR(13:16)).^2);
   TEMP_P_DIF(1,5) = sum((TEMP_P(1,5)-P_FIN_ERR(17:20)).^2);
   TEMP_P_DIF(1,6) = sum((TEMP_P(1,6)-P_FIN_ERR(21:24)).^2);
   len = length(A_ERRS(:,1));
   
   % Combining standard deviations from different channels
   P_FIN_STD = sqrt((4*(len-1).*TEMP_PS+4*len*TEMP_P_DIF)/(4*len-1));
   % Relative standard deviation of Position
   P_FIN_RSTD = (P_FIN_STD./TEMP_P)*100;
      
   % Converet to LaTeX table 
   string_arr1 = [ "",  "$m_{error}(\%)$", "$s^2$", "RSTD (\%)" ];
   string_arr2 = ["GP (PG1)"; "GP (PG2)"; "GP (PG3)"; "Spline"; "PLSQ"; "LM"];
    
   fin2 = [string_arr2 TEMP_A' A_FIN_STD' A_FIN_RSTD'];
   fin2 = [string_arr1; fin2];
   
   
   path = "C:\Users\jil\Desktop\Kaivos\PARAINEN\lets try something\amplitude.tex";
   matrix2latex(fin2, path, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny');
   
   %fin3 = [string_arr2 TEMP_P' P_FIN_STD' P_FIN_RSTD'];
   %fin3 = [string_arr1; fin3];
   
   %path = "C:\Users\jil\Desktop\Kaivos\PARAINEN\lets try something\position.tex";
   %matrix2latex(fin3, path, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny');
  
   %} 
   %% Calculation of errors for waveform fit between measured and fitted data
   % Will have a mean value in each column of a 1D vector for each channel
   % and each fitting method. Ie there will be a total of 24 different
   % values.
   
   %Gaussian:   4 channels * 3 PG   = 12
   %Spline:     4 channels          = 4
   %Poly:       4 channels          = 4
   %LM:         4 channels          = 4
   %-------------------------------------
   %Total:                          = 24
   
   % The order these values are arranged in the array is the same as above
   % strating from PG1 PG2 PG3 then Spline and so on...

   %% Concatenate WF Means
   WF_A_ERRS = [GAME SAME PAME LAME];
   WF_A_STD = [GASD SASD PASD LASD];
   % Remove all rows skipped in the previous loop. The skipped row numbers
   % have been stored in a vector calles skpd_r in the for loop above
 
   WF_A_ERRS(skpd_r,:)  = [];
   WF_A_STD(skpd_r,:) = [];
   n_dat(skpd_r) = [];
   A_ERRS(skpd_r,:) = [];
   
   % Values needed to form point cloud
   trigs2(skpd_r,:) = [];
   peak_pos(skpd_r,:) = [];
   peak_val(skpd_r,:) = [];
   dpx(skpd_r,:) = [];
   dpy(skpd_r,:) = [];
   
   %disp('This should be l_rpp by 24 vector')
   %[wf_ar, wf_ac] = size(WF_A_ERRS)
   disp('length of n_dat and WF_A_ERRS before NaN removal')
   length(n_dat)
   length(WF_A_ERRS(:,1))
   length(A_ERRS)
   %find number of rows with nans in both arrays before and after
   %the removal
   
   %find row numbers with all zeros before removing
   disp('number of all zero rows n_dat and WF_A_ERRS before NaN removal')
   [azr, azc] = find(~all(WF_A_ERRS==0,2));
   [nzr, nzc] = find(~all(n_dat ==0,2));
   len_azr = length(azr);
   len_nzr = length(nzr);
   
   %% Find NANs and remove them from all arrays
   disp('n_dat and WF_A_ERRS before NaN removal')
   length(n_dat)
   length(WF_A_ERRS(:,1))
   
   re = 1-any(isnan(WF_A_ERRS),2);
   
   WF_A_ERRS = WF_A_ERRS.*re;
   WF_A_ERRS = WF_A_ERRS(any(WF_A_ERRS,2),:);
   
   A_ERRS = A_ERRS.*re;
   A_ERRS = A_ERRS(any(A_ERRS,2),:);
   
   WF_A_STD = WF_A_STD.*re;
   WF_A_STD = WF_A_STD(any(WF_A_STD,2),:);
   
   n_dat = n_dat.*re;
   n_dat = n_dat(any(n_dat,2),:);
   
   trigs2 = trigs2.*re;
   trigs2 =trigs2(any(trigs2,2),:);
   
   peak_pos = peak_pos.*re;
   peak_pos = peak_pos(any(peak_pos,2),:);
   
   peak_val = peak_val.*re;
   peak_val = peak_val(any(peak_val,2),:);
   
   % find all zero rows in "re" and use them to remove the values from 
   % data_pos_x and data_pos_y.
   %{
   rr = zeros(length(re),1);
   tracker = 1;
   for i=1:length(re)
        if re(i) == 0
           rr(tracker,1) = i;
           tracker = tracker + 1;
        end
   end
   rr = rr(1:tracker,1);
   dpx(rr(:,1)) = [];
   dpy(rr(:,1)) = [];
   %} 
   disp('n_dat and WF_A_ERRS after NaN removal')
   length(n_dat)
   length(WF_A_ERRS(:,1))
   
   %% Remove inf.s 
   
   disp('n_dat, WF_A_ERRS and A_ERRS before Inf removal')
   length(n_dat)
   length(WF_A_ERRS(:,1))
   length(A_ERRS(:,1))
   
   remz = 1-any(isinf(WF_A_ERRS),2);
   WF_A_ERRS = WF_A_ERRS.*remz;
   WF_A_ERRS = WF_A_ERRS(any(WF_A_ERRS,2),:);
   
   WF_A_STD = WF_A_STD.*remz;
   WF_A_STD = WF_A_STD(any(WF_A_STD,2),:);
   
   A_ERRS = A_ERRS.*remz;
   A_ERRS = A_ERRS(any(A_ERRS,2),:);
   
   n_dat = n_dat.*remz;
   n_dat = n_dat(any(n_dat,2),:);
   
   trigs2 = trigs2.*remz;
   trigs2 =trigs2(any(trigs2,2),:);
   
   peak_pos = peak_pos.*remz;
   peak_pos = peak_pos(any(peak_pos,2),:);
   
   peak_val = peak_val.*remz;
   peak_val = peak_val(any(peak_val,2),:);
   
   rremz = zeros(length(remz),1);
   tracker = 0;
   for i=1:length(remz)
        if remz(i) == 0
           tracker = tracker + 1;
           rremz(tracker,1) = i;
        end
   end
   if tracker > 1
        rremz = rremz(1:tracker,1);
        dpx(rremz(:,1)) = [];
        dpy(rremz(:,1)) = [];
   end
   
   disp('n_dat, WF_A_ERRS and A_ERRS after Inf removal')
   length(n_dat)
   length(WF_A_ERRS(:,1))
   length(A_ERRS(:,1))
   ldpx = length(dpx)
   
   %% Reomve all errors larger than 50% These pulses were not fitted properly
   [re_r,re_c] = find(WF_A_ERRS(:,1)>=50);
   for i =1:length(re_r)
       WF_A_ERRS(i,:) = [];
       WF_A_STD(i,:) = [];
       A_ERRS(i,:) = [];
       n_dat(i,:) = [];
       trigs2(i,:) = [];
       peak_pos(i,:) = [];
       peak_val(i,:) = [];
       dpx(i) = [];
       dpy(i) = [];
       
   end
   
   disp('n_dat, WF_A_ERRS and A_ERRS after 50% error removal')
   length(n_dat)
   length(WF_A_ERRS(:,1))
   length(A_ERRS);
   ldpx = length(dpx)
   
   %% Remove nan.s again using WF_A_STD
   re3 = 1-any(isnan(WF_A_STD),2);
   WF_A_ERRS = WF_A_ERRS.*re3;
   WF_A_ERRS = WF_A_ERRS(any(WF_A_ERRS,2),:);
   
   WF_A_STD = WF_A_STD.*re3;
   WF_A_STD = WF_A_STD(any(WF_A_STD,2),:);
   
   A_ERRS = A_ERRS.*re3;
   A_ERRS = A_ERRS(any(A_ERRS,2),:);
   
   n_dat = n_dat.*re3;
   n_dat = n_dat(any(n_dat,2),:);
   
   trigs2 = trigs2.*re3;
   trigs2 = trigs2(any(trigs2,2),:);
   
   peak_pos = peak_pos.*re3;
   peak_pos = peak_pos(any(peak_pos,2),:);
   
   peak_val = peak_val.*re3;
   peak_val = peak_val(any(peak_val,2),:);
   
   rr3 = zeros(length(re3),1);
   tracker = 1;
   for i=1:length(re3)
        if re3(i) == 0
           rr3(tracker,1) = i;
           tracker = tracker + 1;
        end
   end
   if length(rr3) > 1
    rr3 = rr3(1:tracker-1,1);
    dpx(rr3(:,1)) = [];
    dpy(rr3(:,1)) = [];
   end
   disp('n_dat, WF_A_ERRS and A_ERRS after 50% error removal')
   length(n_dat)
   length(WF_A_ERRS(:,1))
   length(WF_A_STD(:,1))
   length(A_ERRS);
   
   %% Channel and algorithm specific mean errors for WF MaxAmplitude
   WF_FIN_ERR = mean(WF_A_ERRS);
   %disp('This should be 1 by 24 vector')
   %[wf_finr, wf_finc] = size(WF_FIN_ERR) 
   CH_FIN_AERR = mean(A_ERRS);
   CH_FIN_STD = sqrt(sum((CH_FIN_AERR - A_ERRS).^2)/(length(A_ERRS(:,1))-1));
   CH_FIN_RSTD = (CH_FIN_STD./CH_FIN_AERR)*100;
   
   % Calculation of combined standard deviation for WF amplitude
   % Each column of the resulting vector will hold a CSD value expressing 
   % the std of a particular fitting method or PG at a particular channel
   WF_CSD1 = sqrt((sum((n_dat-1).*(WF_A_STD.^2))+ ...
       sum(n_dat.*((WF_A_ERRS-WF_FIN_ERR).^2)))...
       /(sum(n_dat)-1));
   
   % CSD from all channels. WF_CSD2 holds the final standard deviation
   % values for each fitting method with the cahnnel error data combined
   WF_CSD_FIN = zeros(1,6);
   
   %Mean errors for each algorihtm calculated form all channels
   FIN_A_ARR = zeros(1,6);
   FIN_A_ARR(1,1) = mean(WF_FIN_ERR(1:4));      % Gauss PG1
   FIN_A_ARR(1,2) = mean(WF_FIN_ERR(5:8));      % Gauss PG2
   FIN_A_ARR(1,3) = mean(WF_FIN_ERR(9:12));     % Gauss PG3
   FIN_A_ARR(1,4) = mean(WF_FIN_ERR(13:16));    % Spline
   FIN_A_ARR(1,5) = mean(WF_FIN_ERR(17:20));    % Polynomial
   FIN_A_ARR(1,6) = mean(WF_FIN_ERR(21:24));    % Levenberg-Marquardt
   
   % Total number of data points on each channel
   tot_dp = sum(n_dat);
   
   WF_CSD_FIN(1,1) = sqrt(( sum((tot_dp-1)*(WF_CSD1(1:4).^2)) + ...
       sum(tot_dp*((WF_FIN_ERR(1:4)-FIN_A_ARR(1,1)).^2)))...
       /(4*tot_dp-1));
   
   WF_CSD_FIN(1,2) = sqrt(( sum((tot_dp-1)*(WF_CSD1(5:8).^2)) + ...
       sum(tot_dp*((WF_FIN_ERR(5:8)-FIN_A_ARR(1,2)).^2)))...
       /(4*tot_dp-1));
   
   WF_CSD_FIN(1,3) = sqrt(( sum((tot_dp-1)*(WF_CSD1(9:12).^2)) + ...
       sum(tot_dp*((WF_FIN_ERR(9:12)-FIN_A_ARR(1,3)).^2)))...
       /(4*tot_dp-1));
   
   WF_CSD_FIN(1,4) = sqrt(( sum((tot_dp-1)*(WF_CSD1(13:16).^2)) + ...
       sum(tot_dp*((WF_FIN_ERR(13:16)-FIN_A_ARR(1,4)).^2)))...
       /(4*tot_dp-1));
   
   WF_CSD_FIN(1,5) = sqrt(( sum((tot_dp-1)*(WF_CSD1(17:20).^2)) + ...
       sum(tot_dp*((WF_FIN_ERR(17:20)-FIN_A_ARR(1,5)).^2)))...
       /(4*tot_dp-1));
   
   WF_CSD_FIN(1,6) = sqrt(( sum((tot_dp-1)*(WF_CSD1(21:24).^2)) + ...
       sum(tot_dp*((WF_FIN_ERR(21:24)-FIN_A_ARR(1,6)).^2)))...
       /(4*tot_dp-1));
   
   % Implement calculation for relative standard deviation (WF)
   RSTD = (WF_CSD_FIN./FIN_A_ARR)*100;
   
   % Converet to LaTeX table 
   string_arr1 = [ "",  "$m_{error}(\%)$", "$s^2$", "RSTD (\%)" ];
   string_arr2 = ["GP (PG1)"; "GP (PG2)"; "GP (PG3)"; "Spline"; "PLSQ"; "LM"];
    
   fin = [string_arr2 FIN_A_ARR' WF_CSD_FIN' RSTD'];
   fin = [string_arr1; fin];

   %path = "C:\Users\jil\Desktop\Kaivos\PARAINEN\lets try something\out.tex";
   path = "C:\Users\julia\Desktop\work\PARAINEN\output\out.tex";
   matrix2latex(fin, path, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny');

   
   
   %% Convert Single channel data to LaTeX format
   
   % Reorginze data to be channel wise
   CH_FIN_A = zeros(1,24); CH_FIN_S = CH_FIN_A; CH_FIN_RS = CH_FIN_A;
   % tracker and increment variablea
   trk = 1; inc = 0;
   for i=1:4
       for j=1:6
            CH_FIN_A(1,trk) = CH_FIN_AERR(i+inc); 
            CH_FIN_S(1,trk) = CH_FIN_STD(i+inc);
            CH_FIN_RS(1,trk) = CH_FIN_RSTD(i+inc);
            trk = trk +1;
            inc = inc +4;
       end
       inc = 0;
   end    
   
   string_arr3 = ["GP (PG1)"; "GP (PG2)"; "GP (PG3)"; "Spline"; "PLSQ"; "LM";...
       "GP (PG1)"; "GP (PG2)"; "GP (PG3)"; "Spline"; "PLSQ"; "LM"; ...
       "GP (PG1)"; "GP (PG2)"; "GP (PG3)"; "Spline"; "PLSQ"; "LM"; ...
       "GP (PG1)"; "GP (PG2)"; "GP (PG3)"; "Spline"; "PLSQ"; "LM"];
    fin2 = [string_arr3 CH_FIN_A' CH_FIN_S' CH_FIN_RS'];
    fin2 = [string_arr1; fin2];
    
    %path2 = "C:\Users\jil\Desktop\Kaivos\PARAINEN\lets try something\channels.tex";
    path2 = "C:\Users\julia\Desktop\work\PARAINEN\output\channels.tex";
    matrix2latex(fin2, path2, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny');
          
end    

%%%% THIS NEEDS TO BE DONE TO GET THE RESULTS CORRECT FOR
% PG2 of channel 1.

%{
a = find(A_ERRS(:,5)>10);
length(a)
AER = A_ERRS(:,5);
AER(a(:,1)) = [];
l = length(AER)
maer = mean(AER)
AERSTD = sqrt(sum((maer - AER(:,1)).^2)/(l-1))
AERRSTD = (AERSTD./maer)*100
%}    

%end