function ret = log_out(name,data);
fid=fopen('.\输出日志\ans_science.txt','a+');
fprintf(fid,'-------%s实验结果-------\n',name);
formatSpec = '%.2e';
fprintf(fid,'最优值:%s\n',strrep(num2str(data(1),formatSpec),'e','E'));
fprintf(fid,'最差值:%s\n',strrep(num2str(data(2),formatSpec),'e','E'));
fprintf(fid,'均值:%s\n',strrep(num2str(data(3),formatSpec),'e','E'));
fprintf(fid,'标准差:%s\n',strrep(num2str(data(4),formatSpec),'e','E'));
fprintf(fid,'\n');
fclose(fid);
end