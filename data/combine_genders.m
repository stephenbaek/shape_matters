clear all
close all
clc

load('AE_male.mat');
AE_male = (AE - mean(AE))./std(AE);
subjectid_male = subjectid;
clear AE subjectid
load('AE_female2.mat');
AE_female = (AE - mean(AE))./std(AE);
AE_female(:,2) = -AE_female(:,2);
subjectid_female = subjectid;
clear AE subjectid

AE = zeros(length(subjectid_male)+length(subjectid_female), 3);
AE(subjectid_male+1, [2,1]) = AE_male;
AE(subjectid_female+1, :) = AE_female;

save('bodyshape.mat', 'AE');

%%
demographics = readtable('Demographics.xls');
measurements = readtable('Measurements.xls');
caesar = join(demographics, measurements);
clear demographics measurements

%%
