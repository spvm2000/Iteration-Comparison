function Data_o = GTO(Data_i,Data_o)
rand_num=[];                         
ti=clock;                              
X = Data_i.X;                        
Ffun=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Ffun);
IterCurve=zeros(1,Data_i.maxIter);
Pop_Fit=Ffun;                      
Silverback_Score=Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);    
p=0.03;
Beta=3;
w=0.8;
It=0;
GX=X(:,:);

variables_no=Data_i.dim;
lb=Data_i.lb; 
ub=Data_i.ub;
max_iter=Data_i.maxIter;
pop_size=Data_i.pop;
while It<max_iter 
    a=(cos(2*rand)+1)*(1-It/max_iter);
    C=a*(2*rand-1);
    for i=1:pop_size
        if rand<p    
            GX(i,:) =(ub-lb)*rand+lb;
        else  
            if rand>=0.5
                Z = unifrnd(-a,a,1,variables_no);
                H=Z.*X(i,:);   
                GX(i,:)=(rand-a)*X(randi([1,pop_size]),:)+C.*H; 
            else   
                GX(i,:)=X(i,:)-C.*(C*(X(i,:)- GX(randi([1,pop_size]),:))+rand*(X(i,:)-GX(randi([1,pop_size]),:))); %ok ok 
            end
        end
        GX(i,:) = BoundaryCheck(GX(i,:), ub, lb,Data_i.dim);
    end    
    
    for i=1:pop_size
         New_Fit=Data_i.fobj(GX(i,:));
         if New_Fit<Pop_Fit(i)
            Pop_Fit(i)=New_Fit;
            X(i,:)=GX(i,:);
         end

         if New_Fit<Silverback_Score                
            Silverback_Score=New_Fit; 
            Silverback=GX(i,:);
            Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=Silverback_Score;
            Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
         end
    end
    It=It+1;
     if It>max_iter
         break;
     else
        IterCurve(It)=Silverback_Score;
     end

end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=It;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end