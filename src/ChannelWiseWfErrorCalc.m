
   %% Convert Single channel WF data to LaTeX format
   % For the code to work use the saved data files called
   % "ChannelWiseFwfAmpMeanErr" and "ChannelWiseFwfAmpStd"
   
   
   % Reorginze data to be channel wise
   CH_FIN_WF = zeros(1,24); CH_FIN_WFS = CH_FIN_WF; CH_FIN_WFRS = CH_FIN_WF;
   CH_WF_RSTD = (WF_CSD1./WF_FIN_ERR)*100;
   
   % tracker and increment variable
   trk = 1; inc = 0;
   for i=1:4
       for j=1:6
            CH_FIN_WF(1,trk) = WF_FIN_ERR(i+inc); 
            CH_FIN_WFS(1,trk) = WF_CSD1(i+inc);
            CH_FIN_WFRS(1,trk) = CH_WF_RSTD(i+inc);
            trk = trk +1;
            inc = inc +4;
       end
       inc = 0;
   end   
   
   string_arr1 = [ "",  "$m_{error}(\%)$", "$s^2$", "RSTD (\%)" ];
   string_arr3 = ["GP (PG1)"; "GP (PG2)"; "GP (PG3)"; "Spline"; "PLSQ"; "LM";...
       "GP (PG1)"; "GP (PG2)"; "GP (PG3)"; "Spline"; "PLSQ"; "LM"; ...
       "GP (PG1)"; "GP (PG2)"; "GP (PG3)"; "Spline"; "PLSQ"; "LM"; ...
       "GP (PG1)"; "GP (PG2)"; "GP (PG3)"; "Spline"; "PLSQ"; "LM"];
    fin2 = [string_arr3 CH_FIN_WF' CH_FIN_WFS' CH_FIN_WFRS'];
    fin2 = [string_arr1; fin2];
    
    path2 = "C:\Users\jil\Desktop\Kaivos\PARAINEN\lets try something\channelsWF.tex";
    %path2 = "C:\Users\julia\Desktop\work\PARAINEN\output\channels.tex";
    matrix2latex(fin2, path2, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny');
      

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