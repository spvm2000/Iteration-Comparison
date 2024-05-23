function Data_o = CAPSA(Data_i,Data_o)
t_c=clock;                             
cnt=1;                                  
%% Problem Definition
maxite = Data_i.maxIter;      
noP = Data_i.pop;            
dim=Data_i.dim;
CapPos = Data_i.X;              
CapFit = Data_i.F_value;         
LB = Data_i.lb;
UB = Data_i.ub;
cg_curve=zeros(1,maxite);
v=0.1*CapPos;% initial velocity
v0=zeros(noP,dim); 

% Initial Fit of the random Pos

Fit = CapFit;

[fitCapSA,index]=min(CapFit);
Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) = fitCapSA;
Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=index;
CapBestPos = CapPos; 
Pos= CapPos; 
gFoodPos = CapPos(index,:); 
bf=0.70;
cr=11.0;  
g=9.81;

% CapSA velocity updates
a1=1.250; a2=1.5;   

beta=[2 11 2];
wmax=0.8;
wmin=0.1;

  
for t = 1 : maxite

    % Life time convergence
        tau = beta(1) * exp(-beta(2) * t/maxite)^beta(3);
        w   = wmax-(wmax-wmin)*(t/maxite); 
        fol=ceil((noP-1).*rand(noP,1))'; %  

        % CapSA velocity update
for i=1:noP 
    for j=1:dim
        v(i,j)=  w* v(i,j) + ...
                           a1*(CapBestPos(i,j)- CapPos(i,j))*rand + ...
                           a2*(gFoodPos(j)  - CapPos(i,j))*rand;      
    end
end
     
% CapSA Pos update

for i=1:noP
   if i<noP/2
          if (rand()>=0.1)
               r=rand;
              if r<=0.15
                 CapPos(i,:) =  gFoodPos +    bf*((v(i,:).^2)*sin(2*rand()*1.5))/g;  % Jumping (Projection)
              elseif   r>0.15 && r<=0.30  
                  CapPos(i,:) =  gFoodPos +   cr*bf*((v(i,:).^2)*sin(2*rand()*1.5))/g;  % Jumping (Land)  
              elseif   r>0.30 && r<=0.9      
                  CapPos(i,:) =    CapPos(i,:) +  v(i,:); % Movement on the ground    
              elseif  r>0.9 && r<=0.95 
                 CapPos(i,:) =      gFoodPos   +  bf*sin(rand()*1.5);   % Swing % Local search  
              elseif   r>0.95 
                CapPos(i,:) =       gFoodPos   +  bf*(v(i,:)- v0(i,:));    % Climbing   % Local search
              end
           else
               CapPos(i,:) =           tau*(LB  + rand *(UB- LB));   
           end
    
% Let the followers follow the leaders (update their Pos)
elseif i>=noP/2 && i<=noP 
            
           eps=((rand()+2*rand())-(3*rand()))/(1+rand()); 
     
           Pos(i,:)=gFoodPos+2*(CapBestPos(fol(i),:)-CapPos(i,:))*eps +...
                                 2*(CapPos(i,:)-CapBestPos(i,:))*eps;      
 
          CapPos(i,:)=(Pos(i,:)+ CapPos(i-1,:))/(2); 
    
   end
end 
v0 = v;   
 

for i=1:noP % relocation (Update, exploration)
        u=UB-CapPos(i,:)<0;
        l=LB-CapPos(i,:)>0;
     
         CapPos(i,:)= LB.*l+UB.*u+CapPos(i,:).*~xor(u,l);
    
         CapFit(i)=Data_i.fobj (CapPos(i,:)) ;
         
            if CapFit(i)<Fit(i)
                CapBestPos(i,:)=CapPos(i,:);
                Fit(i)=CapFit(i);
            end 
end
[fmin,index]=min(Fit); % finding out the best Pos  

% Updating gPos and best Fit
if fmin < fitCapSA
    gFoodPos = CapBestPos(index,:); % Update the global best Pos
    fitCapSA = fmin;
    Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=index;
end

   cg_curve(t)=fitCapSA;
   Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) = fitCapSA;
   
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=t;                                
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=cg_curve;
end