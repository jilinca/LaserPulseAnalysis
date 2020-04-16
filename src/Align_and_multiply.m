%% ALIGN AND MULTIPLY PEAKS
% aligns corresponding waveforms from different channels. Then multiplies 
% corresponding channels to make peaks more distinct. Multplied waveforms 
% are stored in the 'mult' array. find_returns.m is called to find the 
% locations of the returns


%% Initialization of arrays
rows = zeros(length(ch1_data_(:,1)),4);
cols = rows; 

%% Check if trigger is missing from waveform and remove such rows
% Find trigger peak position for each waveform. Cols holds location of
% value in the array.
[~,cols(:,1)]= max(ch1_data_(:,1:90)');
[~,cols(:,2)]= max(ch2_data_(:,1:90)');
[~,cols(:,3)]= max(ch3_data_(:,1:90)');
[~,cols(:,4)]= max(ch4_data_(:,1:90)');

% Remove waveforms (and corresponding angular data) that don't have trigger
% Threshold for minumvalue of peak location in the array
th = 45;
temp_r = (cols(:,:) < 45);
[r,c] = find(temp_r(:,:) == 1);
I= unique(r);

% Remove rows with no trigger
ch1_data_(I,:)= [];
ch2_data_(I,:)= [];
ch3_data_(I,:)= [];
ch4_data_(I,:)= [];
data_pos_x(I) = []; 
data_pos_y(I) = []; 
cols(I,:)= [];
ch1_time_(I,:)= [];

%{
%ch1_data_= ch1_data_(any(ch1_data_,2),:);
%ch2_data_= ch2_data_(any(ch2_data_,2),:);
%ch3_data_= ch3_data_(any(ch3_data_,2),:);
%ch4_data_= ch4_data_(any(ch4_data_,2),:);
%cols = cols(any(cols,2),:);
%ch1_time_= ch1_time_(any(ch1_time_,2),:);
%}

% Remove rows with only partial trigger
part_trig_r = find(ch1_data_(:,1)>15);
ch1_data_(part_trig_r,:)=[];
ch2_data_(part_trig_r,:)=[];
ch3_data_(part_trig_r,:)=[];
ch4_data_(part_trig_r,:)=[];
ch1_time_(part_trig_r,:)=[];
data_pos_x(part_trig_r)=[];
data_pos_y(part_trig_r)=[];
cols(part_trig_r,:)=[];

%% Align corresponding waveforms
% matrix holding the value telling how much a waveform has to be shifted
% relative to the channel 1 waveform
[r,c] = size(ch1_data_);
shift = zeros(r,3);
shift(:,1) = cols(:,1)-cols(:,2);
shift(:,2) = cols(:,1)-cols(:,3);
shift(:,3) = cols(:,1)-cols(:,4);
% Remove rows with exceedingly large shifts
sh_r1 = find(shift(:,1) < -5);
part_trig_r = find(ch1_data_(:,1)>15); % what is the purpose of this line (an attempt to remove partial triggers)
% remove waveforms with a large shift between columns 1 and 2
if isempty(sh_r1) == 1
    shift(sh_r1,:) = [];
    ch1_data_(sh_r1,:)=[];
    ch2_data_(sh_r1,:)=[];
    ch3_data_(sh_r1,:)=[];
    ch4_data_(sh_r1,:)=[];
    ch1_time_(sh_r1,:)=[];  
    data_pos_x(sh_r1)=[];
    data_pos_y(sh_r1)=[];
    cols(sh_r1,:)=[];
    r = r-length(sh_r1);
end
% remove waveforms with a large shift between columns 1 and 3
sh_r2 = find(shift(:,2) < -5);
if isempty(sh_r2) == 1
    shift(sh_r2,:) = [];
    ch1_data_(sh_r2,:)=[];
    ch2_data_(sh_r2,:)=[];
    ch3_data_(sh_r2,:)=[];
    ch4_data_(sh_r2,:)=[];
    ch1_time_(sh_r2,:)=[];
    data_pos_x(sh_r2)=[];
    data_pos_y(sh_r2)=[];
    cols(sh_r2,:)=[];
    r = r-length(sh_r2);
end
sh_r3 = find(shift(:,3) < -5);
% remove waveforms with a large shift between columns 1 and 4
if isempty(sh_r3) == 1
    shift(sh_r3,:) = [];
    ch1_data_(sh_r3,:)=[];
    ch2_data_(sh_r3,:)=[];
    ch3_data_(sh_r3,:)=[];
    ch4_data_(sh_r3,:)=[];
    ch1_time_(sh_r3,:)=[];
    data_pos_x(sh_r3)=[];
    data_pos_y(sh_r3)=[];
    cols(sh_r3,:)=[];
    r = r-length(sh_r3);
end    

% Assign Arrays for aligned waveforms
ch1_d = zeros(r,c); ch2_d = ch1_d; ch3_d =ch1_d; ch4_d = ch1_d; 
% Align the waveforms
for rs =1:r
    ch2_d(rs,:) = circshift(ch2_data_(rs,:), [0 shift(rs,1)]);
    ch3_d(rs,:) = circshift(ch3_data_(rs,:), [0 shift(rs,2)]);
    ch4_d(rs,:) = circshift(ch4_data_(rs,:), [0 shift(rs,3)]);
end
% ch1_t acts a time refernce for channels 1-4, since now all the triggers 
% have been shifted to correspond to the trigger location on channel 1
ch1_d = ch1_data_;
ch1_t = ch1_time_;

%% Remvoe all rows where all data values are less than 20
[~,m1] = max(ch1_d,[],2);

rm1 = m1(m1<20);

ch1_d(rm1,:) =[];
ch2_d(rm1,:) =[];
ch3_d(rm1,:) =[];
ch4_d(rm1,:) =[];
ch1_t(rm1,:) =[];
data_pos_x(rm1) =[];
data_pos_y(rm1) =[];
cols(rm1,:) =[];
r = r-length(rm1);
%% Multiply corresponding waveform arrays
% and all data points by them selves to make peaks distinct
% helps in locating far away peaks
mult = ch1_d.*(ch2_d.^3);
mult = mult.*(ch3_d.^3);
mult = mult.*(ch4_d.^3);

%% PLot a non-fitted waveform that still contains both the trigger and
%{
% echoes
%save path for figures
s_path = '\\tsclient\C\Users\jil\Desktop\Projects and Byrocracy\FGI Research Seminar\';

figure(204626)
hold on
plot(ch1_t(17568,:),ch1_d(17568,:),'LineWidth',2)
plot(ch1_t(17568,:),ch2_d(17568,:),'LineWidth',2)
plot(ch1_t(17568,:),ch3_d(17568,:),'LineWidth',2)
plot(ch1_t(17568,:),ch4_d(17568,:),'LineWidth',2)
xlabel('time (ns)')
ylabel('Amplitude')
hold off


figure(204627)
plot(ch1_t(17568,:),ch1_d(17568,:),'LineWidth',2)
xlabel('time (ns)')
ylabel('Amplitude')
%}
%% Fit all triggers and extract max amplitudes from channels 1-4
% The 'trigs' array is arranged so that the first column holds the 
% spatio-temporal trigger peak location and the columns after that hold the
% channel wise peak amplitude value of the trigger.
trigs = zeros(r,5);
%pp = zeros(4,1);
pp = zeros(4,r);
for i =1:r
    tx = ch1_t(i,cols(i,1)-6:cols(i,1)+6);
    
    %Check that triggers are actually aligned
    %{
    if i == 3
        figure(7165)
        hold on
        plot(tx,ch1_d(i,cols(i,1)-6:cols(i,1)+6))
        plot(tx,ch2_d(i,cols(i,1)-6:cols(i,1)+6))
        plot(tx,ch3_d(i,cols(i,1)-6:cols(i,1)+6))
        plot(tx,ch4_d(i,cols(i,1)-6:cols(i,1)+6))
        hold off
        figure(7175)
        plot(ch1_t(3,:),ch1_d(3,:))
    end
    
    if i == 54653
        figure(7166)
        hold on
        plot(tx,ch1_d(i,cols(i,1)-6:cols(i,1)+6))
        plot(tx,ch2_d(i,cols(i,1)-6:cols(i,1)+6))
        plot(tx,ch3_d(i,cols(i,1)-6:cols(i,1)+6))
        plot(tx,ch4_d(i,cols(i,1)-6:cols(i,1)+6))
        hold off
        figure(7176)
        plot(ch1_t(54653,:),ch1_d(54653,:))
    end
    
    if i == 105426
        figure(7167)
        hold on
        plot(tx,ch1_d(i,cols(i,1)-6:cols(i,1)+6))
        plot(tx,ch2_d(i,cols(i,1)-6:cols(i,1)+6))
        plot(tx,ch3_d(i,cols(i,1)-6:cols(i,1)+6))
        plot(tx,ch4_d(i,cols(i,1)-6:cols(i,1)+6))
        hold off
        figure(7177)
        plot(ch1_t(105426,:),ch1_d(105426,:))
    end
    
    if i == 135842
        figure(7168)
        hold on
        plot(tx,ch1_d(i,cols(i,1)-6:cols(i,1)+6))
        plot(tx,ch2_d(i,cols(i,1)-6:cols(i,1)+6))
        plot(tx,ch3_d(i,cols(i,1)-6:cols(i,1)+6))
        plot(tx,ch4_d(i,cols(i,1)-6:cols(i,1)+6))
        hold off
        figure(7178)
        plot(ch1_t(135842,:),ch1_d(135842,:))
    end
    
    if i == 169999
        figure(7169)
        hold on
        plot(tx,ch1_d(i,cols(i,1)-6:cols(i,1)+6))
        plot(tx,ch2_d(i,cols(i,1)-6:cols(i,1)+6))
        plot(tx,ch3_d(i,cols(i,1)-6:cols(i,1)+6))
        plot(tx,ch4_d(i,cols(i,1)-6:cols(i,1)+6))
        hold off
        figure(7179)
        plot(ch1_t(169999,:),ch1_d(169999,:))
    end
    %}
    
    % Use a spline fit to retreive the maximum trigger value and 
    % temporal location
    xx = linspace(tx(1),tx(length(tx)),100);
    [trigs(i,2),p1]=max(spline(tx,ch1_d(i,cols(i,1)-6:cols(i,1)+6),xx));
    [trigs(i,3),p2]=max(spline(tx,ch2_d(i,cols(i,1)-6:cols(i,1)+6),xx));
    [trigs(i,4),p3]=max(spline(tx,ch3_d(i,cols(i,1)-6:cols(i,1)+6),xx));
    [trigs(i,5),p4]=max(spline(tx,ch4_d(i,cols(i,1)-6:cols(i,1)+6),xx));
    %calculate mean temporal location for trigger peaks
    
    % Check that fitted peaks match the actual triggers
    %{
    if i == 3
        figure(7185)
        hold on
        plot(xx,spline(tx,ch1_d(i,cols(i,1)-6:cols(i,1)+6),xx))
        plot(xx,spline(tx,ch2_d(i,cols(i,1)-6:cols(i,1)+6),xx))
        plot(xx,spline(tx,ch3_d(i,cols(i,1)-6:cols(i,1)+6),xx))
        plot(xx,spline(tx,ch4_d(i,cols(i,1)-6:cols(i,1)+6),xx))
        hold off
    end
    
    if i == 54653
        figure(7186)
        hold on
        plot(xx,spline(tx,ch1_d(i,cols(i,1)-6:cols(i,1)+6),xx))
        plot(xx,spline(tx,ch2_d(i,cols(i,1)-6:cols(i,1)+6),xx))
        plot(xx,spline(tx,ch3_d(i,cols(i,1)-6:cols(i,1)+6),xx))
        plot(xx,spline(tx,ch4_d(i,cols(i,1)-6:cols(i,1)+6),xx))
        hold off
    end
    if i == 105426
        figure(7187)
        hold on
        plot(xx,spline(tx,ch1_d(i,cols(i,1)-6:cols(i,1)+6),xx))
        plot(xx,spline(tx,ch2_d(i,cols(i,1)-6:cols(i,1)+6),xx))
        plot(xx,spline(tx,ch3_d(i,cols(i,1)-6:cols(i,1)+6),xx))
        plot(xx,spline(tx,ch4_d(i,cols(i,1)-6:cols(i,1)+6),xx))
        hold off
    end
    
    if i == 135842
        figure(7188)
        hold on
        plot(xx,spline(tx,ch1_d(i,cols(i,1)-6:cols(i,1)+6),xx))
        plot(xx,spline(tx,ch2_d(i,cols(i,1)-6:cols(i,1)+6),xx))
        plot(xx,spline(tx,ch3_d(i,cols(i,1)-6:cols(i,1)+6),xx))
        plot(xx,spline(tx,ch4_d(i,cols(i,1)-6:cols(i,1)+6),xx))
        hold off
    end

    if i == 169999
        figure(7189)
        hold on
        plot(xx,spline(tx,ch1_d(i,cols(i,1)-6:cols(i,1)+6),xx))
        plot(xx,spline(tx,ch2_d(i,cols(i,1)-6:cols(i,1)+6),xx))
        plot(xx,spline(tx,ch3_d(i,cols(i,1)-6:cols(i,1)+6),xx))
        plot(xx,spline(tx,ch4_d(i,cols(i,1)-6:cols(i,1)+6),xx))
        hold off
    end
    %}
    
    pp(1,i) = xx(p1);
    pp(2,i) = xx(p2);
    pp(3,i) = xx(p3);
    pp(4,i) = xx(p4);
    trigs(i,1)=mean(pp(:,i));
end
%{
trigs = zeros(r,5);
for i =1:r
    trigs(i,1) = ch1_t(i,cols(i,1));
    trigs(i,2) = ch1_d(i,cols(i,1));
    trigs(i,3) = ch2_d(i,cols(i,1));
    trigs(i,4) = ch3_d(i,cols(i,1));
    trigs(i,5) = ch4_d(i,cols(i,1));
end
%}

%% Remove triggers (so that finding returns becomes easier
ch1_d(:,1:90)=[];
ch2_d(:,1:90)=[];
ch3_d(:,1:90)=[];
ch4_d(:,1:90)=[];
ch1_t(:,1:90)=[];
mult(:,1:90)=[];
% Remove useless data from end (so that processing becomes faster)
% if there are far-away returns comment this out 
%(this is not strictly needed in any case)
%%{
[~,len] =size(ch1_d);
ch1_d(:,len-10:len)=[];
ch2_d(:,len-10:len)=[];
ch3_d(:,len-10:len)=[];
ch4_d(:,len-10:len)=[];
ch1_t(:,len-10:len)=[];
mult(:,len-10:len)=[]; 
%%}

% Remove rows with no or extremly small returns far away returns
% that are uninteresting

%{
ch1_dtemp = ch1_d;
ch1_dtemp(ch1_dtemp < 10) = 0;
ch1_sum = sum(ch1_dtemp,2);
[deler,~] = find(ch1_sum == 0);
%}


%% Extract return postions
% Find peaks by calculating slopes: Returns postions of peaks. The
% amplitude and time values can be then extracted from the arrays above
[ret_pos, skipped, skipped1, ret_pos_r_nums, rows_to_be_removed, temp2, i2,r_mult, r2_mult]...
    = find_return_peak_locations(mult,ch1_t);


% Remove the pulses from proccessing for which no returns were found
% The pulse numbers are stored in the array 'skipped' and the array 'skipped1'
% holds a logical vector for the rows that can be removed from the data i.e
% the skipped rows.
% The pulses for which returns were not found are not stored in ret_pos
% by default

% Remove all zeros from the "skipped" array to only have the row numbers
% of the waveforms that were skipped. Use this array in a for loop to get
% rid of data in "data_pos_x" and "data_pos_y"

skipped3 = zeros(length(skipped),1);
trk = 1;
for i=1:length(skipped)
    if skipped(i) ~= 0
        skipped3(trk,1) = skipped(i);
        trk = trk +1;
    end
end    

% Remove all zeros froom the end of "skipped3"
skipped3 = skipped3(1:trk-1,1);


% Conctenate the channels
con = horzcat(ch1_t,ch1_d, ch2_d, ch3_d, ch4_d); 
% remove skipped from all neccessay arrays
con = con.*skipped1;
con = con(any(con,2),:);
mult = mult.*skipped1;
mult = mult(any(mult,2),:);
trigs = trigs.*skipped1;
trigs = trigs(any(trigs,2),:);

% Remove data from data_pos_x and data_pos_y using "skipped3"

trk2 = 0;
for i=1:length(skipped3(:,1))
    data_pos_x(skipped3(i,1)-trk2) = [];
    data_pos_y(skipped3(i,1)-trk2) = [];
    trk2 = trk2+1;
end

%data_pos_x = data_pos_x.*skipped1;
%data_pos_x = data_pos_x(any(data_pos_x,2),:);
%data_pos_y = data_pos_y.*skipped1;
%data_pos_y = data_pos_y(any(data_pos_y,2),:);
% Bring data and time arrays back to original form and place into struct.
% reshape(array_to_be_reshaped, row_size_of_new_arrays, column_size_of_new_arrays, n_of_new_arrays)
B=reshape(con,length(con(:,1)),length(ch1_t(1,:)),5); 
% 'darrs' will hold all the arrays with the pulses to be processed
str = ["td1","ad1","ad2","ad3","ad4"];
for k=1:5
    darrs.(sprintf('%s',str(k)))=B(:,:,k);
end

% CHECK THAT THE NUMBER OF DIFFERET WAVEFORMS STORED IN RET POS IS THE 
% SAME AS THE NUMBER OF WAVEFORMS IN EACH CAHNNEL OF CONS.

arr = find(skipped1==0);
chk = length(ch1_d(:,1)) - length(arr);
if (chk == length(darrs.ad1(:,1)))
    disp("YEY")
end    
%% Extract actual return amplitude and spatio-temporal postion
%[] = process_for_fitting(darrs,ret_pos);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%