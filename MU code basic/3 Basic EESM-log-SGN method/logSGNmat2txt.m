clear all
load('snr_LogSGNParam_Config192_Model-B_1-by-1_MCS5.mat')
len = length(snrs);
fid = fopen('snr_LogSGNParam_Config192_Model-B_1-by-1_MCS5.txt', 'w');
for rowIdx = 1 : len
    fprintf(fid, '{%4.2f, {%6.4f, %6.4f, %6.4f, %6.4f}},\n', snrs(rowIdx), logSGNParam(rowIdx,:));
end
fclose(fid);