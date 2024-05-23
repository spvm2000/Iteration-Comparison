function Data_o = BDO(Data_i,Data_o)
t_c=clock;                             
cnt=1;                                  
%% Problem Definition
M = Data_i.maxIter;      
pop = Data_i.pop;            
dim=Data_i.dim;
lb= Data_i.lb;    % Lower limit/bounds/     a vector
ub= Data_i.ub;    % Upper limit/bounds/     a vector
x = Data_i.X;      
fit = Data_i.F_value;                
P_percent = 0.2;    % The population size of producers accounts for "P_percent" percent of the total population size       
pNum = round( pop *  P_percent );    % The population size of the producers   
Convergence_curve=zeros(1,Data_i.maxIter);
pFit = fit;                       
pX = x; 
XX=pX;    
[ fMin, bestI ] = min( fit );      % fMin denotes the global optimum fitness value
bestX = x( bestI, : );             % bestX denotes the global optimum position corresponding to fMin
Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=fMin;
Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=bestI;
% Start updating the solutions.
for t = 1 : M    
       
       [~,B]=max(fit);
       worse= x(B,:);   
       r2=rand(1);  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1 : pNum    
        if(r2<0.9)
            r1=rand(1);
          a=rand(1,1);
          if (a>0.1)
           a=1;
          else
           a=-1;
          end
       x( i , : ) =  pX(  i , :)+0.3*abs(pX(i , : )-worse)+a*0.1*(XX( i , :)); % Equation (1)
       else
            
           aaa= randperm(180,1);
           if ( aaa==0 ||aaa==90 ||aaa==180 )
            x(  i , : ) = pX(  i , :);   
           end
           theta= aaa*pi/180;   
           x(  i , : ) = pX(  i , :)+tan(theta).*abs(pX(i , : )-XX( i , :));    % Equation (2)      

        end
      
        x(  i , : ) = Bounds( x(i , : ), lb, ub );    
        fit(  i  ) = Data_i.fobj( x(i , : ) );
    end 
  [ ~, bestII ] = min( fit );      % fMin denotes the current optimum fitness value
  bestXX = x( bestII, : );             % bestXX denotes the current optimum position 

  R=1-t/M;                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Xnew1 = bestXX.*(1-R); 
   Xnew2 =bestXX.*(1+R);                    %%% Equation (3)
   Xnew1= Bounds( Xnew1, lb, ub );
   Xnew2 = Bounds( Xnew2, lb, ub );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Xnew11 = bestX.*(1-R); 
   Xnew22 =bestX.*(1+R);                     %%% Equation (5)
   Xnew11= Bounds( Xnew11, lb, ub );
   Xnew22 = Bounds( Xnew22, lb, ub );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    for i = ( pNum + 1 ) :12                  % Equation (4)
     x( i, : )=bestXX+((rand(1,dim)).*(pX( i , : )-Xnew1)+(rand(1,dim)).*(pX( i , : )-Xnew2));
   x(i, : ) = Bounds( x(i, : ), Xnew1, Xnew2 );
  fit(i ) = Data_i.fobj(  x(i,:) ) ;
   end
   
  for i = 13: 19                  % Equation (6)

   
        x( i, : )=pX( i , : )+((randn(1)).*(pX( i , : )-Xnew11)+((rand(1,dim)).*(pX( i , : )-Xnew22)));
       x(i, : ) = Bounds( x(i, : ),lb, ub);
       fit(i ) = Data_i.fobj(  x(i,:) ) ;
  
  end
  
  for j = 20 : pop                 % Equation (7)
       x( j,: )=bestX+randn(1,dim).*((abs(( pX(j,:  )-bestXX)))+(abs(( pX(j,:  )-bestX))))./2;
      x(j, : ) = Bounds( x(j, : ), lb, ub );
      fit(j ) = Data_i.fobj(  x(j,:) ) ;
  end
   % Update the individual's best fitness vlaue and the global best fitness value
     XX=pX;
    for i = 1 : pop 
        if ( fit( i ) < pFit( i ) )
            pFit( i ) = fit( i );
            pX( i, : ) = x( i, : );
        end
        
        if( pFit( i ) < fMin )
           % fMin= pFit( i );
           fMin= pFit( i );
           bestX = pX( i, : );
          %  a(i)=fMin;
           Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i; 
        end
    end
    Convergence_curve(t)=fMin;   
    Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) = fMin;
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);           
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=t;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=Convergence_curve;
end

% Application of simple limits/bounds
function s = Bounds( s, Lb, Ub)
  % Apply the lower bound vector
  temp = s;
  I = temp < Lb;
  temp(I) = Lb(I);
  
  % Apply the upper bound vector 
  J = temp > Ub;
  temp(J) = Ub(J);
  % Update this new move 
  s = temp;
end

function S = Boundss( SS, LLb, UUb)
  % Apply the lower bound vector
  temp = SS;
  I = temp < LLb;
  temp(I) = LLb(I);
  
  % Apply the upper bound vector 
  J = temp > UUb;
  temp(J) = UUb(J);
  % Update this new move 
  S = temp;
end
%---------------------------------------------------------------------------------------------------------------------------
