%====================== savesignals.m ====================================
% save data of a day pair into directorysave2daysignals
%==================== Warning ===============
% we have observed huge outliers in the following files:
% ihc==1, date==2015/10/07 ,from sample index 2.4e6
% ihc==1, date==2015/10/09 ,from sample index 2.4e6
% ihc==2, date==2015/08/07 ,from sample index 2.4e6
% ihc==2, date==2015/10/05 ,from sample index 2.4e6
% ihc==5, date==2015/10/05
% ihc==6, date==2015/10/07
% ihc==8, date==2015/10/07 , from sample index 2.5e6

% This program does not remove them for examination purposes
%
%============================================
setimesE_ihc                 = zeros(1,2);
setimesN_ihc                 = zeros(1,2);
setimesZ_ihc                 = zeros(1,2);
problemHC                    = zeros(1,2);

idSE1 = 1;
idSZ1 = 1;
idSN1 = 1;
idSE2 = 1;
idSZ2 = 1;
idSN2 = 1;

idT  = 1;
signals     = zeros(34560000,6);
Lrecords    = length(records);

for ir = 1:Lrecords
    switch records{ir}.station
        case 'PD31'
            switch records{ir}.channel(3)
                case 'E'
                    auxE1 = records{ir}.etime;
                    LLE1 = length(records{ir}.data);
                    signalsE1 = [records{ir}.data];
                    signals(idSE1:idSE1+LLE1-1,1)=signalsE1;
                    idSE1 = idSE1+LLE1;
                    
                case 'N'
                    auxN1 = records{ir}.etime;
                    LLN1 = length(records{ir}.data);
                    signalsN1 = [records{ir}.data];
                    signals(idSN1:idSN1+LLN1-1,2)=signalsN1;
                    idSN1 = idSN1+LLN1;
                    
                case 'Z'
                    auxZ1 = records{ir}.etime;
                    LLZ1 = length(records{ir}.data);
                    signalsZ1 = [records{ir}.data];
                    signals(idSZ1:idSZ1+LLZ1-1,3)=signalsZ1;
                    idSZ1 = idSZ1+LLZ1;
            end
        case 'PD32'
            switch records{ir}.channel(3)
                case 'E'
                    auxE2 = records{ir}.etime;
                    LLE2 = length(records{ir}.data);
                    signalsE2 = [records{ir}.data];
                    signals(idSE2:idSE2+LLE2-1,4)=signalsE2;
                    idSE2 = idSE2+LLE2;
                    
                case 'N'
                    auxN2 = records{ir}.etime;
                    LLN2 = length(records{ir}.data);
                    signalsN2 = [records{ir}.data];
                    signals(idSN2:idSN2+LLN2-1,5)=signalsN2;
                    idSN2 = idSN2+LLN2;
                    
                case 'Z'
                    auxZ2 = records{ir}.etime;
                    LLZ2 = length(records{ir}.data);
                    signalsZ2 = [records{ir}.data];
                    signals(idSZ2:idSZ2+LLZ2-1,6)=signalsZ2;
                    idSZ2 = idSZ2+LLZ2;
                    
                    
            end
    end
end

idScMin = min([idSE1, idSN1, idSZ1, idSE2, idSN2, idSZ2]);

signals     = signals(1:idScMin-1,:);
Fs_Hz       = records{1}.Fs_Hz;
Ts_sec      = 1/Fs_Hz;
signals_centered = signals-ones(size(signals,1),1)*mean(signals);

commandsave = sprintf...
    ('save %syear%smonth%sday%s signals_centered Fs_Hz', ...
    directorysave2daysignals,yearstart,monthstart,daystart);
eval(commandsave)




