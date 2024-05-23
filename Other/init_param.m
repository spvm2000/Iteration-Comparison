function [Data_i,Data_o] = init_param(F_name)
%% 所有参数初始化
Data_i=struct('param',struct(),'pop',30,'maxIter',40,'test_error',1E-06,'dim',30,'test_iter',30,'fobj','', ...
    'ub',[],'lb',[],'X',[],'now_test_iter',0,'now_fun',"",'Ai_name',[],'ALLR',[],'F_value',[] ...
    ,'now_fun_index',1,'file_no',2,'auto_flag',0,'Color_ary',[],'total_test_file',["输入数据\total.txt"]);
Data_i.file_no=Data_i.maxIter;
data=readlines(Data_i.total_test_file);
if length(data)~=length(unique(data))
    error('所选的测试文件列表中存在重复函数名!');
else
    Data_i.Ai_name=data';
end
Data_i.now_fun=F_name;
Data_i.Color_ary=['r','g','b','y'];                             
Data_i.save_num=1;                                                      
Data_o=struct('fun_best_fit',0,'Best_Pos',[],'Best_fitness',[],'exe_time',[],'min_t',[],'IterCurve',[],'Ai_mean',[], ...
    'Ai_std',[],'Ai_best',[],'Ai_worst',[],'Time_mean',[],'Item_mean',[],'meanCurve',[],'Ai_result',[], ...
    'F',[],'F_value_x',[],'F_value_y',[],'F_value_z',[],'history_position',[], ...
    'test_history_best',struct(),'history_best',struct(),'massey_param',struct());
Data_o.test_history_best=struct('fit',[],'best_x',[]);
Data_o.history_best=struct('fit',[],'x',[],'item_no',[],'test_no',[],'best_x',[]);
[Data_i.lb,Data_i.ub,Data_i.dim,Data_i.fobj]=Get_F(F_name);                
Data_i.ub = Data_i.ub.*ones(1,Data_i.dim);                   
Data_i.lb = Data_i.lb.*ones(1,Data_i.dim);
Data_o.Best_fitness=ones(length(Data_i.Ai_name),Data_i.test_iter).*100000;    
Data_i=rand_init(Data_i);                                              
for i=1:Data_i.pop              
    Data_i.F_value(i)=Data_i.fobj(Data_i.X(i,:));
end   
end