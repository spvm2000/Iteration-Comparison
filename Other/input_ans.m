function  input_ans(Data_i,Data_o)
% 结果输出函数
path=['.\实验数据\实验结果(',num2str(Data_i.file_no),')\'];
if ~isfolder(path)
    mkdir(['.\实验数据\实验结果(',num2str(Data_i.file_no),')']);
end    

path1=['.\实验数据\实验结果(',num2str(Data_i.file_no),')\',num2str(Data_i.now_fun)];
if ~isfolder(path1)
    mkdir(['.\实验数据\实验结果(',num2str(Data_i.file_no),')\',num2str(Data_i.now_fun)]);
end

%2.保存日志文件
copyfile('.\输出日志\ans_science.txt',['.\实验数据\实验结果(',num2str(Data_i.file_no),')\',num2str(Data_i.now_fun)]);
close all;
end