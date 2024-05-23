function Data_o = HPO(Data_i,Data_o)
t_c=clock;                             
cnt=1;                                  
%% HPO Parameters
MaxIt = Data_i.maxIter;     % Maximum Nomber of Iterations
nPop = Data_i.pop;         % Population Size
dim = Data_i.dim;             
HPpos = Data_i.X;               
HPposFitness = Data_i.F_value;         
lb = Data_i.lb;
ub = Data_i.ub;
Convergence_curve = zeros(1,MaxIt);

% Constriction Coefeicient
B = 0.1;
% NFE = nPop;
[TargetScore,indx] = min(HPposFitness);
Target = HPpos(indx,:);   % Target HPO
Convergence_curve(1)=TargetScore;
Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=TargetScore;
Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=indx;

%% HPO Main Loop
for it = 2:MaxIt

   c = 1 - it*((0.98)/MaxIt);   % Update C Parameter
    kbest=round(nPop*c);        % Update kbest
     for i = 1:nPop
            r1=rand(1,dim)<c;
            r2=rand;
            r3=rand(1,dim);
            idx=(r1==0);
            z=r2.*idx+r3.*~idx;
%             r11=rand(1,dim)<c;
%             r22=rand;
%             r33=rand(1,dim);
%             idx=(r11==0);
%             z2=r22.*idx+r33.*~idx;
        if rand<B
        xi=mean(HPpos);
        dist = pdist2(xi,HPpos);
        [~,idxsortdist]=sort(dist);
        SI=HPpos(idxsortdist(kbest),:);
        HPpos(i,:) =HPpos(i,:)+0.5*((2*(c)*z.*SI-HPpos(i,:))+(2*(1-c)*z.*xi-HPpos(i,:)));
        else
          for j=1:dim
            rr=-1+2*z(j);
          HPpos(i,j)= 2*z(j)*cos(2*pi*rr)*(Target(j)-HPpos(i,j))+Target(j);

          end
        end  
        HPpos(i,:) = min(max(HPpos(i,:),lb),ub);
        % Evaluation
        HPposFitness(i) = Data_i.fobj(HPpos(i,:));
        % Update Target
        if HPposFitness(i)<TargetScore 
            Target = HPpos(i,:);
            TargetScore = HPposFitness(i);
            if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) >TargetScore
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
                Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) =TargetScore;
            end
        end
     end 
  Convergence_curve(it)=TargetScore;
 end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=it;                                
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=Convergence_curve;
end

