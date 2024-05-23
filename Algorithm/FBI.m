function Data_o = FBI(Data_i,Data_o)
rand_num=[];                         
ti=clock;                               
Pop = Data_i.X;                        
ObjVal=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
LB=Data_i.lb;       UB=Data_i.ub;   NP=Data_i.pop;  GEN = Data_i.maxIter;      
iBest = find(ObjVal==min(ObjVal));    D=Data_i.dim;
iBest = iBest(end); 
GlobalMin = ObjVal(iBest); 
Xbest = Pop(iBest,:); 
g = 1;
Fit_store = zeros(1,GEN);
PopA = Pop; ObjValA = ObjVal;
PopB = Pop; ObjValB = ObjVal;
while g <= GEN
    for i=1:NP 
        Change=fix(rand*D)+1;
        nb1=floor(rand*NP)+1;      
            while(nb1==i)
                nb1=floor(rand*NP)+1;
            end
        nb2=floor(rand*NP)+1;      
            while(nb1==nb2 || nb2==i)
                nb2=floor(rand*NP)+1;
            end   
       solA=PopA(i,:);
       solA(Change)=PopA(i,Change)+(PopA(i,Change)-(PopA(nb1,Change)+PopA(nb2,Change))/2)*(rand-0.5)*2; %Eq.(2) in FBI Inspired Meta-Optimization
       for Change = 1:D
            if (solA(Change)<=LB(Change))||(solA(Change)>=UB(Change))
                solA(Change) = LB(Change) + (UB(Change)-LB(Change))*rand();
            end
       end
        f_a = Data_i.fobj(solA);
        if f_a <= ObjValA(i)
            PopA(i,:) = solA; 
            ObjValA(i) = f_a ;
            if f_a <= GlobalMin
                Xbest = solA ; 
                GlobalMin = f_a; 
            end  
        end
    end

    if min(ObjValA) < max(ObjValA)
    prob = probability(ObjValA);
    for i = 1:NP
        if (rand>prob(i))
            r(1) = floor(rand()* NP) + 1;
            while r(1)==i
                r(1) = floor(rand()* NP) + 1;
            end
            r(2) = floor(rand()* NP) + 1;
            while (r(2)==r(1))||(r(2)==i)
                r(2) = floor(rand()* NP) + 1;
            end
             r(3) = floor(rand()* NP) + 1; 
             while (r(3)==r(2))||(r(3)==r(1))||(r(3)==i)
                 r(3) = floor(rand()* NP) + 1; 
             end 
            solA = PopA(i,1:D);
            Rnd = floor(rand()*D) + 1;
            for j = 1:D
                if (rand()< rand()) || ( Rnd==j)
                    solA(j) = Xbest(j) + PopA(r(1),j) + rand() * (PopA(r(2),j) - PopA(r(3),j)); %Eq.(5) in FBI Inspired Meta-Optimization
                else
                    solA(j) = PopA(i,j);
                end
            end
            for j = 1:D
                if (solA(j)<=LB(j))||(solA(j)>=UB(j))
                    solA(j) = LB(j) + (UB(j)-LB(j))*rand();
                end
            end
            f_a = Data_i.fobj(solA);
            if f_a <= ObjValA(i)
                PopA(i,:) = solA; 
                ObjValA(i) = f_a;
                if f_a <= GlobalMin
                    Xbest = solA ; 
                end  
            end
        end 
    end
    end

    for i = 1:NP 
        SolB = PopB(i,1:D);
        for j = 1:D
            SolB(j) = rand()*PopB(i,j) + rand()*(Xbest(j) - PopB(i,j)); %Eq.(6) in FBI Inspired Meta-Optimization
            if (SolB(j)<LB(j))||(SolB(j)>UB(j))
                SolB(j) = LB(j) + (UB(j)-LB(j))*rand();
            end
        end
        f_b = Data_i.fobj(SolB);
        if f_b <= ObjValB(i)
            PopB(i,1:D) = SolB; 
            ObjValB(i) = f_b;
            if f_b <= GlobalMin
                Xbest = SolB ; 
                GlobalMin = f_b; 
            end
        end
    end 

    for i = 1:NP 
        rr = floor(rand()* NP) + 1;
        while rr ==i
            rr = floor(rand()* NP) + 1;
        end
        if ObjValB(i) > ObjValB(rr)
            SolB = PopB(rr,:) + rand(1,D).* (PopB(rr,:) - PopB(i,:))+ rand(1,D).* (Xbest - PopB(rr,:)); %Eq.(7) in FBI Inspired Meta-Optimization
            for j = 1:D
                if (SolB(j)<LB(j))||(SolB(j)>UB(j))
                    SolB(j) = LB(j) + (UB(j)-LB(j))*rand();
                end
            end
        else
            SolB = PopB(i,:) + rand(1,D).* (PopB(i,:)-PopB(rr,:))+ rand(1,D).* (Xbest - PopB(i,:)); %Eq.(8) in FBI Inspired Meta-Optimization
            for j = 1:D
                if (SolB(j)<LB(j))||(SolB(j)>UB(j))
                    SolB(j) = LB(j) + (UB(j)-LB(j))*rand();
                end
            end
        end
        f_b = Data_i.fobj(SolB);
        if f_b <= ObjValB(i)
            PopB(i,:) = SolB; 
            ObjValB(i) = f_b;
            if f_b <= GlobalMin
                Xbest = SolB ; 
                GlobalMin = f_b; 
            end
            if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>f_b
                Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=f_b;
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
            end
        end
    end 
    IterCurve(g) = GlobalMin;
    g = g + 1;
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=g;                              
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function f = rescale_matrix(X, LB, UB)
[NP,D] = size(X);
f = zeros(NP,D);
for i = 1:D
    f(:, i) = LB(i)*ones(NP,1) + (UB(i) - LB(i))*X(:,i);
end
end

function prob = probability(fObjV)
prob = (max(fObjV)-fObjV)/((max(fObjV)-min(fObjV))); %Eq.(3) in FBI Inspired Meta-Optimization
end
