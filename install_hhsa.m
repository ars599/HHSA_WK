function [] = install_hhsa
% before run this script, please decompress the split spm12.zip files by winzip
nss=which('holo_eeg_tf_par');
[filepath,name,ext] = fileparts(nss);
sep=filesep;
addpath(filepath);
unzip('spm12',filepath);
fi=fullfile(filepath,'IFS');
ft=fullfile(filepath,'statistics');
fs=fullfile(filepath,'spm12');
addpath(fs,fi,ft);
spm eeg
spm('Quit');
savepath

