%% The main function of this program is to calculate errors of charateristic pressures and RMSE of waveforms.
clc
clear all
close all
load('T.mat');
load('P_est.mat');
load('reAP.mat');
fss=1000;
t_est=0:1/fss:T;
AP=reAP(t_est);
SBP_origin=roundn(max(AP),-1);
SBP_est=roundn(max(P_est),-1);
DBP_origin=roundn(min(AP),-1);
DBP_est=roundn(min(P_est),-1);
MBP_origin=roundn(mean(AP),-1);
MBP_est=roundn(mean(P_est),-1);
PP_origin=roundn(SBP_origin-DBP_origin,-1);
PP_est=roundn(SBP_est-DBP_est,-1);
Error_SBP=roundn(SBP_est-SBP_origin,-1);
Error_DBP=roundn(DBP_est-DBP_origin,-1);
Error_MBP=roundn(MBP_est-MBP_origin,-1);
Error_PP=roundn(PP_est-PP_origin,-1);
RMSE= roundn(sqrt(mean((P_est-AP').^2)),-1);
Origin=[SBP_origin DBP_origin MBP_origin PP_origin ];
Est=[SBP_est DBP_est MBP_est PP_est RMSE Error_SBP Error_DBP Error_MBP Error_PP];
Error=[Origin Est];
