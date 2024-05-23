
% D=10, 9, 16, 18,F1 is defined on D=9，
% F2 is defined on D=16,
% F3 is defined on D=18
%if (n==10)
    
        %mexPrintf ("usage: f = cec19_func(x, func_num);\n");
		%mexErrMsgTxt ("n=10;");
    
% F4-F10 are defined on D=10%

for i=4:10
    x=rand(1,10);
    
    o(i-3,:)=ZHX2019(x',i);
end
o