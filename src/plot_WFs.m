%% WAVEFORM PLOT
s_path = "C:\Users\jil\Desktop\Kaivos\PAPERI\FINAL_OSA_ARTICLE";
hf=figure(1)
plot(ch1_time_(5354,1:400), ch1_data_(5354,1:400), 'LineWidth',2)
hold on
plot(ch1_time_(5354,1:400), ch2_data_(5354,1:400), 'LineWidth',2)
plot(ch1_time_(5354,1:400), ch3_data_(5354,1:400), 'LineWidth',2)
plot(ch1_time_(5354,1:400), ch4_data_(5354,1:400), 'LineWidth',2)

xlim([10 80])
ylim([-10 100])
legend('Channel 1','Channel 2','Channel 3','Channel 4');
xlabel('time (ns)','FontSize',28)
ylabel('Amplitude','FontSize',28)
grid on
set(gca,'FontSize',28)
          
% Saving Figure
set(hf,'PaperSize',fliplr(get(hf,'PaperSize'))) 
print(hf,[s_path 'WaveformsCh1-4'],'-dpdf')