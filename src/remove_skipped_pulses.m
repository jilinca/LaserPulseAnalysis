%% Preprocess Data and removes extra angles
function [ch1_data_, ch1_time_, ch2_data_, ch2_time_, data_pos_x, data_pos_y]...
    = remove_skipped_pulses(data_time, timeStamp_1, timeStamp_2,...
    ch1_data_, ch1_time_, ch2_data_,ch2_time_, data_pos_x, data_pos_y)
%#ok<*AGROW>
%% REMOVAL OF CORRUPTED DATA
% for some reason the 42 first recoded pulses of the Parainen data set 
% have corrupted information and thus are removed. Feel free to comment
% out this section if your data integrity is fine.
    %%{
    if length(timeStamp_1(:,1)) == 178858 % This is the number of waveforms before processing is started. Change accordingly
        ch1_data_(1:42,:)=[];
        %{
    ch2_data_(1:42,:) = [];
    ch3_data_(1:42,:) = [];
    ch4_data_(1:42,:) = [];
    ch5_data_(1:42,:) = [];
    ch6_data_(1:42,:) = [];
    ch7_data_(1:42,:) = [];
    ch8_data_(1:42,:) = [];
    %}
        ch1_time_(1:42,:)=[];
        %{
    ch2_time_(1:42,:) = [];
    ch3_time_(1:42,:) = [];
    ch4_time_(1:42,:) = [];
    ch5_time_(1:42,:) = [];
    ch6_time_(1:42,:) = [];
    ch7_time_(1:42,:) = [];
    ch8_time_(1:42,:) = [];
    %}
        timeStamp_1(1:42)=[];
        timeStamp_2(1:42)=[];
        data_time(1:42)=[];
        data_pos_x(1:42)=[];
        data_pos_y(1:42)=[];
    end
    %%}
    %%{
%% COMPARISON AND REMOVAL OF EXTRA ANGLES AND WAVEFORMS
% NOTE: When dealing with all 8 channels the array timeStamp_1 for
% channels 1-4 and the array timeStamp_2 for channels 5-8. For this to 
% go correctly a if clause hase to be implemented below

    j=1; i=1;
    i_del =[];
    limit = length(data_time(:,1))-length(timeStamp_1(:,1));
    while i <= length(timeStamp_1(:,1))    
        if timeStamp_1(i,1) - data_time(j) > 1
            i_del=[i_del,j];
            i=i-1;
        end
        if j == length(data_time)
            rts = length(i_del)-limit;
            break
        end     
        j=j+1;
        i=i+1;
    end

    for i=1:length(i_del)
        data_pos_x(i_del(i))=1000;
        data_pos_y(i_del(i))=1000;
    end 

    data_pos_x(data_pos_x==1000)=[];
    data_pos_y(data_pos_y==1000)=[];

    if rts > 0
        timeStamp_1(length(timeStamp_1)-rts+1:length(timeStamp_1))=[];
        %data_pos_x(length(data_pos_x)-rts:length(data_pos_x))=[];
        %data_pos_y(length(data_pos_y)-rts:length(data_pos_y))=[];
        ch1_data_(length(ch1_data_(:,1))-rts+1:length(ch1_data_(:,1)),:)=[];
        ch1_time_(length(ch1_time_(:,1))-rts+1:length(ch1_time_(:,1)),:)=[];
        ch2_data_(length(ch2_data_(:,1))-rts+1:length(ch2_data_(:,1)),:)=[];
        ch2_time_(length(ch2_time_(:,1))-rts+1:length(ch2_time_(:,1)),:)=[];        
    end
    %%}
    if length(data_pos_x(:,1)) ~= length(timeStamp_1)
        warning('Something went wrong in function "remove_skipped_pulses.m". The number of angles does not match the number of waveforms')
    else
        disp('Measurement angles have been set to correspond to their relative wavefroms.') 
        disp('Done')
    end

%% REMOVE IRRELEVANT DATA FROM BEGGINGING AND END OF EACH WAVEFORM
% this can be commented out if all data points in each waveform need to be
% analyzed
    %%{
    len = length(ch1_data_(1,:));
    ch1_data_(:,len-200:len) = [];
    %{
ch2_data_(:,len-200:len) = [];
ch3_data_(:,len-200:len) = [];
ch4_data_(:,len-200:len) = [];
ch5_data_(:,len-200:len) = [];
ch6_data_(:,len-200:len) = [];
ch7_data_(:,len-200:len) = [];
ch8_data_(:,len-200:len) = [];
%}
    ch1_time_(:,len-200:len) = [];
    %{
ch2_time_(:,len-200:len) = [];
ch3_time_(:,len-200:len) = [];
ch4_time_(:,len-200:len) = [];
ch5_time_(:,len-200:len) = [];
ch6_time_(:,len-200:len) = [];
ch7_time_(:,len-200:len) = [];
ch8_time_(:,len-200:len) = [];
%}
    ch1_data_(:,1:40) = [];
    %{
ch2_data_(:,1:40) = [];
ch3_data_(:,1:40) = [];
ch4_data_(:,1:40) = [];
ch5_data_(:,1:40) = [];
ch6_data_(:,1:40) = [];
ch7_data_(:,1:40) = [];
ch8_data_(:,1:40) = [];
%}
    ch1_time_(:,1:40) = [];
    %{
ch2_time_(:,1:40) = [];
ch3_time_(:,1:40) = [];
ch4_time_(:,1:40) = [];
ch5_time_(:,1:40) = [];
ch6_time_(:,1:40) = [];
ch7_time_(:,1:40) = [];
ch8_time_(:,1:40) = [];
%}
    %%}
%%END OF FUNCTION
end
