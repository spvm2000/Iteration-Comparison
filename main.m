clc;clear all;
warning off;
close all force;                    
addpath(genpath(pwd));                                  
F_name='F18';                                           
[Data_i,Data_o] = init_param(F_name);                
h=waitbar(0,'Ready progress......');
for i = 1:Data_i.test_iter
    Data_i.now_test_iter=i;                             
    for j=1:size(Data_i.Ai_name,2)
        Data_i.now_fun_index=j;                         
        I_name=str2func(Data_i.Ai_name(j));              
        str=['当前实验进度:',num2str(floor((((i-1)*size(Data_i.Ai_name,2)+j)/(Data_i.test_iter*size(Data_i.Ai_name,2)))*100)),'%[',num2str(j),'.',func2str(I_name),'/',F_name,']'];
        waitbar(i/Data_i.test_iter,h,str);
        Data_o = I_name(Data_i,Data_o);
    end
end
pause(0.1);close(h);                                    
fid=fopen('.\输出日志\ans_science.txt','w'); fclose(fid);
for j=1:size(Data_i.Ai_name,2)
    Data_o.Ai_mean(j)=mean(Data_o.Best_fitness(j,:));                           
    Data_o.Ai_std(j)=std(Data_o.Best_fitness(j,:));         
    Data_o.Ai_best(j)=min(Data_o.Best_fitness(j,:));        
    Data_o.Ai_worst(j)=max(Data_o.Best_fitness(j,:));       
    Data_o.Ai_result(j,:)=[Data_o.Ai_best(j),Data_o.Ai_worst(j),Data_o.Ai_mean(j),Data_o.Ai_std(j)];
    Data_o.Time_mean(j)=mean(Data_o.exe_time(j,:));             
    log_out(Data_i.Ai_name(j),Data_o.Ai_result(j,:));
    Data_i.ALLR=[Data_i.ALLR;Data_o.Ai_result(j,:)];
end                                      
cost_input(Data_i.now_fun,Data_i.Ai_name,Data_o.Time_mean,Data_o.Item_mean,Data_i.file_no);
input_ans(Data_i,Data_o);
close all;
rmpath(genpath(strcat(pwd,'\实验数据\实验结果(',num2str(Data_i.save_num),')\',F_name)));      
disp('实验结束!');