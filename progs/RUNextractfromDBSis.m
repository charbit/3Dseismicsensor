%=============== RUNextractfromDB.m ======================================
clear

addpath   00pierrick/

directorysave2daysignals  = './';

%=== temporary files
temporary_gparse_dir = '00pierrick/tempfiles/';
if not(exist(temporary_gparse_dir,'dir'))
    rmdir(temporary_gparse_dir,'s')
    mkdir(temporary_gparse_dir)
else
    commandrmparse = sprintf('!rm %s.*',temporary_gparse_dir);
    eval(commandrmparse);
end
if exist('gparse.wfdisc','file')
    !rm gparse*.*;
end
%========= source of data
data_source = 'operations_archive';
user        = 'charbit';
password    = 'sqlmomo';
channel     = '(''BHE'',''BHN'',''BHZ'',''SHE'',''SHN'',''SHZ'')';

yearstart   =  '2015';
monthstart  =  '10';
HMSstart    = '00:00:10';
yearend     =  '2015';
monthend    =  '10';
HMSend      = '05:00:10';

stations    = '(''PD31'',''PD32'')';%'(''FIA1'',''FIA1'')';%(''ZAA0'',''ZAA0B'')';%
daystart_num    =  22;
if daystart_num<10
    daystart    = ['0' num2str(daystart_num)];
else
    daystart    = num2str(daystart_num);
end
%=== clean temporary files
commandclean = sprintf('!rm %s/*.*',temporary_gparse_dir);
eval(commandclean)
%=== extract data from the database
h_starttime   = sprintf('%s/%s/%s %s',yearstart,monthstart,daystart, HMSstart);
h_endtime     = sprintf('%s/%s/%s %s',yearend,monthend,daystart, HMSend);
[~,starttime] = unix(['h2e ',h_starttime,' ofmt="%#"']);
[~,endtime]   = unix(['h2e ',h_endtime,' ofmt="%#"']);
starttime     = str2double(starttime);
endtime       = str2double(endtime);
wlength       = endtime-starttime;
extractfromDBSis
if nowfdiscflag
    display('***** .wfdisc does not exist');
else
    savesignalsSis
end
if exist(sprintf('%sgparse.wfdisc',temporary_gparse_dir),'file')
    commandrmparse = sprintf('!rm %s*.*',temporary_gparse_dir);
    eval(commandrmparse);
end
%=================================================================

