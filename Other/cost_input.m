function cost_input(F_name,Ai_name,Time_mean,Item_mean,num)
path=['.\实验数据\实验结果(',num2str(num),')\实验代价\',num2str(F_name)];
if ~isfolder(path)
    mkdir(['.\实验数据\实验结果(',num2str(num),')\实验代价\',num2str(F_name)]);   
end
fid1=fopen(['.\实验数据\实验结果(',num2str(num),')\实验代价\',num2str(F_name),'\time_cost.txt'],'w+');
for i=1:size(Ai_name,2)
    fprintf(fid1,'%s平均时间花费:%f\n',Ai_name(i),Time_mean(i));
end
end