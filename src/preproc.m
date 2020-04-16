% Preprocessing
% Remove the first 42 rows from the channel datasets
%function preproc()

ch1_data_(1:42,:) = [];
%{
ch2_data_(1:42,:) = [];
ch3_data_(1:42,:) = [];
ch4_data_(1:42,:) = [];
ch5_data_(1:42,:) = [];
ch6_data_(1:42,:) = [];
ch7_data_(1:42,:) = [];
ch8_data_(1:42,:) = [];
%}
ch1_time_(1:42,:) = [];
%{
ch2_time_(1:42,:) = [];
ch3_time_(1:42,:) = [];
ch4_time_(1:42,:) = [];
ch5_time_(1:42,:) = [];
ch6_time_(1:42,:) = [];
ch7_time_(1:42,:) = [];
ch8_time_(1:42,:) = [];
%}

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
%end