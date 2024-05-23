function Data_o = GJO(Data_i,Data_o)
rand_num=[];                        
ti=clock;                            
Positions = Data_i.X;                        
Cost=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);   
dim=Data_i.dim;
Male_Jackal_pos=zeros(1,dim);
Male_Jackal_score=inf; 

Female_Jackal_pos=zeros(1,dim);  
Female_Jackal_score=inf; 
l=1;
while l<=Data_i.maxIter
    E1=1.5*(1-(l/Data_i.maxIter));
    RL=0.05*levy(Data_i.pop,dim,1.5);
    for i=1:size(Positions,1)
        for j=1:size(Positions,2)     
            r1=rand(); % r1 is a random number in [0,1]
            E0=2*r1-1;            
            E=E1*E0;
            if abs(E)<1
                D_male_jackal=abs((RL(i,j)*Male_Jackal_pos(j)-Positions(i,j))); 
                Male_Positions(i,j)=Male_Jackal_pos(j)-E*D_male_jackal;
                D_female_jackal=abs((RL(i,j)*Female_Jackal_pos(j)-Positions(i,j))); 
                Female_Positions(i,j)=Female_Jackal_pos(j)-E*D_female_jackal;
            else
                D_male_jackal=abs( (Male_Jackal_pos(j)- RL(i,j)*Positions(i,j)));
                Male_Positions(i,j)=Male_Jackal_pos(j)-E*D_male_jackal;
                D_female_jackal=abs( (Female_Jackal_pos(j)- RL(i,j)*Positions(i,j)));
                Female_Positions(i,j)=Female_Jackal_pos(j)-E*D_female_jackal;
            end
            Positions(i,j)=(Male_Positions(i,j)+Female_Positions(i,j))/2;
        end  
    end
    for i=1:Data_i.pop
        Positions(i,:)=BoundaryCheck(Positions(i,:),Data_i.ub,Data_i.lb,Data_i.dim);
        fitness=Data_i.fobj(Positions(i,:));
        if fitness<Male_Jackal_score 
            Male_Jackal_score=fitness; 
            Male_Jackal_pos=Positions(i,:);
        end  
        if fitness>Male_Jackal_score && fitness<Female_Jackal_score 
            Female_Jackal_score=fitness; 
            Female_Jackal_pos=Positions(i,:);
        end
        if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>Male_Jackal_score
            Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=Male_Jackal_score;
            Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
        end
    end
    IterCurve(l)=Male_Jackal_score;
    l=l+1;
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=l;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function [z] = levy(n,m,beta)

    num = gamma(1+beta)*sin(pi*beta/2); % used for Numerator 
    
    den = gamma((1+beta)/2)*beta*2^((beta-1)/2); % used for Denominator

    sigma_u = (num/den)^(1/beta);% Standard deviation

    u = random('Normal',0,sigma_u,n,m); 
    
    v = random('Normal',0,1,n,m);

    z =u./(abs(v).^(1/beta));
end
