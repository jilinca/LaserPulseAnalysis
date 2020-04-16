
s_path = 'C:\Users\jil\Desktop\Kaivos\PARAINEN\Codes\Plots\';

% Extract some noise data
thresh = 10;
nd  = ch1_data_(:,400:500);
nd = ch1_data_<thresh

minh =min(nd);
maxh =max(nd);
edges = linspace(minh,maxh+(maxh/800),800);
figure(1)
histogram(nd,edges);
%title('Error distribution for Max Amplitude','FontSize',16)
xlabel('Amplitude [mV]','FontSize',16);
ylabel('N','FontSize',16)
% Saving Figure
fig = gcf;
fig.PaperPosition = [1 1 20 20];
print([s_path 'NoiseDistCh1'],'-dpng', '-r0')

figure(2)
histfit(nd,800,'Normal')