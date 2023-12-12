clc
clear

p_range = [350,0,0,0;550,90,180,3];  %设计变量的取值范围   对应参数分别为 轨道高度h 轨道倾角i 升交点赤经Raan 以及walker星座相位因子F
p_discrete = [1,0.1,0.1,1]; %设计变量的离散程度，0表示是连续值
Ini_gene_method = 0; %采用第一种初始化方法，随机产生初值
p0_ini_para = 0; %第一种方法无需这个参数
NO = 1; %目标数量

% 启动STK
uiap = actxserver('STK11.application');
root = uiap.Personality2;
root.NewScenario('walker');

pop = PSO_MO_main(p_range,p_discrete,Ini_gene_method,p0_ini_para,NO);