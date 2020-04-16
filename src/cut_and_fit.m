
%% PARMETERS FOR FINDING AMPLITUDE AND LOCATION

% ret_pos holds the positional indexes of the mult array for where 
% the return peak begins to rise until the peak point is reached. The first
% column also holds the waveform number/index.
l_rp = length(ret_pos(:,1)); % finds = ret_pos
%normalized amplitude
norms = zeros(4,l_rp);
range = zeros(1,l_rp);
%mean spatio-temporal location calculated from 4 channels
%m_time = zeros(l_rp,1);
p_amp = zeros(1,4);
pa_loc = zeros(4,1);
p_time = zeros(4,1);
rl = length(ret_pos(1,:));

%% PLOTTING PARAMETERS
% transform angular to cartesian
len = length(data_pos_x(:,1));
carts = zeros(l_rp,3);
transz = sin((pi/180)*data_pos_y);
transx = cos((pi/180)*data_pos_x).*cos((pi/180)*data_pos_y);
transy = sin((pi/180)*data_pos_x).*cos((pi/180)*data_pos_y);

%% CALCULATION OF NORMALIZERD AMPLITUDE, RANGE AND COORDINATE TRANSFORMS

for i=1:l_rp
   % positions in waveform
   
   pos_w = ret_pos(ret_pos(i,2:rl)~=0);
   lpw = length(pos_w);
   pos_w = ret_pos(i,2:length(pos_w)+1);

   cut_d1 = ch1_d(ret_pos(i,1),pos_w(1):pos_w(lpw)+lpw);
   cut_d2 = ch2_d(ret_pos(i,1),pos_w(1):pos_w(lpw)+lpw);
   cut_d3 = ch3_d(ret_pos(i,1),pos_w(1):pos_w(lpw)+lpw);
   cut_d4 = ch4_d(ret_pos(i,1),pos_w(1):pos_w(lpw)+lpw);
   % The time values from channel one are only needed.
   cut_t = ch1_t(ret_pos(i,1),pos_w(1):pos_w(lpw)+lpw);
   % fit using cubic spline
   xx = linspace(cut_t(1),cut_t(length(cut_t)),100);
   s1 = spline(cut_t,cut_d1,xx);
   s2 = spline(cut_t,cut_d2,xx);
   s3 = spline(cut_t,cut_d3,xx);
   s4 = spline(cut_t,cut_d4,xx);
   if i == 178661
    figure(11)
    plot(xx,s1,'LineWidth',2)
    figure(12)
    plot(xx,s2,'LineWidth',2)
    figure(13)
    plot(xx,s3,'LineWidth',2)
    figure(14)
    plot(xx,s4,'LineWidth',2)
    figure(15)
    plot(ch1_t(i,:),ch1_d(ret_pos(i,1),:),'LineWidth',2)
    figure(16)
    plot(ch1_t(i,:),ch2_d(ret_pos(i,1),:),'LineWidth',2)
    figure(17)
    plot(ch1_t(i,:),ch3_d(ret_pos(i,1),:),'LineWidth',2)
    figure(18)
    plot(ch1_t(i,:),ch4_d(ret_pos(i,1),:),'LineWidth',2)
    figure(19)
    plot(cut_t,cut_d1,'Linewidth',2)
   end
   % find peak maximum (p_amp) and its temporal location (pa_loc)
   [p_amp(1,1),pa_loc(1,1)] = max(s1);
   [p_amp(1,2),pa_loc(2,1)] = max(s2);
   [p_amp(1,3),pa_loc(3,1)] = max(s3);
   [p_amp(1,4),pa_loc(4,1)] = max(s4);
   p_time(1,1) = xx(pa_loc(1,1));
   p_time(2,1) = xx(pa_loc(2,1));
   p_time(3,1) = xx(pa_loc(3,1));
   p_time(4,1) = xx(pa_loc(4,1));
   % Mean of the spatio-temporal location
   m_time= mean(p_time);
   % normalize with corresponding trigger 
   % all channels are normalized at once
   norms(:,i) = p_amp./trigs(ret_pos(i,1),2:5);
   range(i) = ((m_time-trigs(ret_pos(i,1),1))*2.99792e-1)/2;
   carts(i,3) = range(i)*transz(ret_pos(i,1));
   carts(i,2) = range(i)*transy(ret_pos(i,1));
   carts(i,1) = range(i)*transx(ret_pos(i,1));
end
