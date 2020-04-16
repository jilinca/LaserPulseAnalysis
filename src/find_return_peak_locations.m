%% Parameters and variables
% Finds the position (coordinates) of return peaks in each waveform
% Allocating space for return locations

%% A MINOR BUG EXISTS IN THIS FUNCTION, SEE NOTES ON PAGE 9.

function [finds, skipped, skipped1, ret_pos_r_nums, rows_to_be_removed, temp2,i, r_mult, r2_mult] = find_return_peak_locations(mult, ch1_t)
finds = zeros(300000,25); 
i =1;   % index for finds array (also indicates the number of returns recorded)
m=1;    % Waveform number
skipped = zeros(length(mult(:,1)),1);
skipped1 = ones(length(mult(:,1)),1);

% all rows in mult that contain only zeros
r_mult = 1-all(mult==0,2);
r2_mult = find(r_mult == 0);

% Calculate slopes
% This operation is done to make the smaller peaks distictive
mult_a = circshift(mult, [0 1]);
time = circshift(ch1_t, [0,1]);
slp = (mult - mult_a)./(ch1_t-time); 
% the first slp will be erronous due to the circshift operation
%slp(:,1) = []; (this line is not needed)
% the second column will represent the slope between the second and the
% first value.
% Thus the position index at which the peak occurs is 1 + the postion index
% of where the slp still passes the threshold below.
%% Peak position finder
while m <= length(mult(:,1)) 
    [~, pos] = find(slp(m,:) > 1e11);   % THIS THRESHOLD NEEDS TO BE CHANGED DEPENDING ON THE DATA SET USED (works for channels 1-4)
    %check if there are multiple returns
    check = pos-circshift(pos,[0,1]); % check will be empty if no peaks are found
   
    % NO PEAKS ARE FOUND
    if isempty(check) %==1 % condition for when no peaks are found
        %store row number of skipped pulses
        skipped(m,1) = m; 
        skipped1(m,1) =0;
        %finds(i,1) =m;
        %i=i+1;
        m=m+1;
        
    % AT LEAST ONE PEAK IS FOUND    
    else % else if isempty(check) ==0
        check(1) = [];
        %holds the positions in which the postion value of the peak location
        %changes in the pos vector
        %check will hold only 1.s if there is only one return
        check_p = find(check ~= 1); %check if there are multiple returns
        % Thus check_p will be empty if there is only one return
        % Check will have the length of 1 if there are two returns
        
        % STORING MULTIPLE RETURN PEAKS
        if isempty(check_p) == 0   
            % number of returns in waveform
            n_of_rets = length(check_p)+1; % +1 was added here as a correction to account for the first return
            % "a" Holds the index values for the index values at which
            % "pos" indicates at which index values of "slp" a new pulse
            % begins
            a = zeros(1,n_of_rets);
            a(1,1)=1;
            % The length of check_p is equal to the number of returns -1
            % The last value of check_p indicates from where the last
            % return begins -1
            for k=1:length(check_p)
                a(1,k+1) = check_p(1,k)+1;
            end
            n =0;
            % The values in check_p indicates the locations in the 'finds'
            % array (which holds postion values for the peaks in slp)
            for j=1:length(a)
                % This is mean to store the last peak
                if j == length(a) && length(pos)-(a(j)-1) > 3 %(try changing from 2 to three)
                    % the second conditon checks that there at least three
                    % consecutive positive slopes
                    finds(i,1) = m;
                    finds(i,2:1+length(pos(a(j):length(pos)))) = pos(a(j):length(pos)); %pos(a(j):length(pos)-1);

                    i=i+1;
                    n=n+1;

                % This is for all peaks except for the last    
                elseif j < length(a) && a(j+1)-a(j) > 3   % --> changing this number to a larger value will define how many consecutive slopes are needed to store a peak. The larger the number the more peaks will be left out 
                    % the second conditon checks that there at least three
                    % consecutive positive slopes per pulse
                    finds(i,1) = m;
                    % a(j) the position index for where a new peak begins
                    % check_p(j) holds the position index for where a peak
                    % ends
                    finds(i,2:length(pos(a(j):check_p(j)))+1) = pos(a(j):check_p(j)); 
                    i=i+1;
                    n=n+1;
                end
               
                % Add a skipped waveform if all pulses found had less than
                % three consecutive positive slopes
                if j == length(a) && n == 0
                    skipped(m,1) = m; 
                    skipped1(m,1) =0;
                end
            end
            
        % STORING THE ONE AND ONLY RETURN PEAK
        % used only if there is only one return
        else  
            % Check that the found peak is longer than 3 consecutive
            % positive slopes
            if (length(check)+1) > 3  
                finds(i,1) =m;
                finds(i,2:length(pos)+1) = pos;
                i=i+1;
            else
                skipped(m,1) = m; 
                skipped1(m,1) =0;
            end
        end
        m=m+1;
    end
end

%% Remove all rows containing only zeros
% The last rows of 'finds' are padded with zeros, all other rows contain 
% both the waveform number 'm' and some position information

temp = finds(:,2:length(finds(1,:)-1));
ret_pos_r_nums = finds(:,1);
%find rows containing all zeros 
% temp2 will contain a zero on every row where all values finds are zeros
temp2 = 1-all(temp==0,2);
rows_to_be_removed = find(temp2 == 0);
%Remove rows containing all zeros
finds = finds(any(finds,2),:);
 
end
