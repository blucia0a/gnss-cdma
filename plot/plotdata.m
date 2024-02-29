%% function [] = plotdata()
% ==============================================================================
%> @brief Plot tracking results of the operation
%> @details A simple routine to plot the 6 raw array values which are generated
%> by the executables from the `install` directory. This should provide insight
%> into the approximate scales and ranges of the values computed, which may aid
%> in optimization).
%> @date      January 30th, 2018
%> @author    Damian Miralles
% ==============================================================================
function [] = plotdata()

%% Cleaning workspace
close all;
clear;

%% File configuration settings
%--- File configuration settings
conf.sw_version = 'Nominal implementation';
conf.comp_flags = 'O3';
conf.cpu_time = 'N/A';
conf.data_proc = '37 s';
conf.path = 'data_reg/';

%% Read data from file
%open the data file - this assumes that both this m-file and the
%output file from the lab4 program are in the current working
%directory for matlab
%and
%read in all the data - assuming it is saved as a double if any
%optimization attempts results in saving the data in a different
%format then that format must be specified below
% (see 'help fread' from within Matlab)
%--- Add path to find files

fid=fopen([conf.path,'/Q_E_output.bin'],'rb');
data_Q_E=fread(fid,'double');

fid=fopen([conf.path,'/Q_P_output.bin'],'rb');
data_Q_P=fread(fid,'double');

fid=fopen([conf.path,'/Q_L_output.bin'],'rb');
data_Q_L=fread(fid,'double');

fid=fopen([conf.path,'/I_E_output.bin'],'rb');
data_I_E=fread(fid,'double');

fid=fopen([conf.path,'/I_P_output.bin'],'rb');
data_I_P=fread(fid,'double');

fid=fopen([conf.path,'/I_L_output.bin'],'rb');
data_I_L=fread(fid,'double');

%assume all has been read in correct and close the data file
fclose(fid);

%% Plotting procedure
%--- Read time in seconds
timeAxisInSeconds = (1:length(data_I_P))/1000;

%create a figure and plot the raw data values (this should
%provide some idea about the range necessary to save the data)
figure(1)
set(gca,'FontSize',10)
suptitle({['Version: ', conf.sw_version,', Flags: ',conf.comp_flags,...
       ', CPU Time: ',conf.cpu_time,', Data Processed: ',conf.data_proc],''});
subplot(2,3,1),plot(timeAxisInSeconds, data_I_E)
ylabel('IE');xlabel('Time (s)');
grid on;
subplot(2,3,2),plot(timeAxisInSeconds, data_I_P)
ylabel('IP');xlabel('Time (s)');
grid on;
subplot(2,3,3),plot(timeAxisInSeconds, data_I_L)
ylabel('IL');xlabel('Time (s)');
grid on;
subplot(2,3,4),plot(timeAxisInSeconds, data_Q_E)
ylabel('QE');xlabel('Time (s)');
grid on;
subplot(2,3,5),plot(timeAxisInSeconds, data_Q_P)
ylabel('QP');xlabel('Time (s)');
grid on;
subplot(2,3,6),plot(timeAxisInSeconds, data_Q_L)
ylabel('QL');xlabel('Time (s)');
grid on;

%create a 2nd figure and plot the 'more interesting' values
%any optimization should provide similiar 'interesting' values
% (try zooming on on the top plot to see actual data bits if
% interested - dma will be happy to explain these values if interested)
figure(2)
set(gca,'FontSize',10)
suptitle({['Version: ', conf.sw_version,', Flags: ',conf.comp_flags,...
       ', CPU Time: ',conf.cpu_time,', Data Processed: ',conf.data_proc],''});
subplot(211);
plot(timeAxisInSeconds, data_I_P,'c');
hold on
plot(timeAxisInSeconds, data_Q_P,'m');
grid on;
title('Decoded Data Bits');
xlabel('Time (s)');
legend('IP', 'QP');

subplot(212);
plot(timeAxisInSeconds, (data_I_E .* data_I_E) + (data_Q_E .* data_Q_E),'b');
hold on;
plot(timeAxisInSeconds, (data_I_P .* data_I_P) + (data_Q_P .* data_Q_P),'g');
plot(timeAxisInSeconds, (data_I_L .* data_I_L) + (data_Q_L .* data_Q_L),'r');
title('Correlation Results');
grid on;
xlabel('Time (s)');
cr_leg = legend( '$\sqrt{I_{E}^2 + Q_{E}^2}$', ...
        '$\sqrt{I_{P}^2 + Q_{P}^2}$', ...
        '$\sqrt{I_{L}^2 + Q_{L}^2}$');
set(cr_leg, 'Interpreter', 'Latex');
