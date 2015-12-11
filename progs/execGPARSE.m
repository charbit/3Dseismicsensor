%====================================================
% this program executes GPARSE, with
% gparse_temp.par file (considered as temporary)
%
% to extract the signal informations from
% the database. Results are saved into 
%      - gparse.wfdisc 
%      - gparse.w 
% written in this directory (considered as temporary)
% usually after we run CONVERTCSS2matlab.m
% 
%====================================================
clear all
close all
h_starttime = '2015/04/22 05:00';
h_endtime   = '2015/04/22 20:00';
% sta         = '(''ZAA0'',''ZAA0B'')';
% sta         = '(''PD31'',''PD32'')';
sta         = '(''FIA1'',''FIA1'')';
channel     = '(''BHZ'',''HHZ'')';

% Time transformation performing
[~,starttime] = unix(['h2e ',h_starttime,' ofmt="%#"']);
[~,endtime]   = unix(['h2e ',h_endtime,' ofmt="%#"']);
starttime     = str2double(starttime);
endtime       = str2double(endtime);
wlength       = endtime-starttime;

% Gparse script edition:
fid=fopen('gparse_temp.par','w');
%Write the query
fprintf(fid, '%s\n', 'open data_source=testbed_archive user=charbit password=sqlmomo');
fprintf(fid, '%s\n',['query wfdisc select * from sel3.wfdisc where sta in ', ...
              sta, 'and chan in ', channel,' and time between ',num2str(starttime),...
              ' and ',num2str(endtime+wlength),' order by sta,chan,time']);
fprintf(fid, '%s\n','read waveforms');
fprintf(fid, '%s\n','write waveforms');
fclose(fid);

% Gparse run:
unix('setenv ORACLE_HOME /cots/oracle/oracle-10.2;');
unix('setenv D_LIBRARY_PATH $ORACLE_HOME/lib:$ORACLE_HOME/lib32;');
unix('/ctbto/ims/sm/local/linux/Geotool++/2.3.10/bin/gparse<gparse_temp.par;');
%====================== END ============================