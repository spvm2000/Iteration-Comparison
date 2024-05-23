%format long;
fhd=@cec22_test_func;

D=20;
for fun=1:12
    %eval(['load input_data/shift_data_' num2str(fun) '.txt']);
    %eval(['O=shift_data_' num2str(fun) '(1,1:D);']);
    O = rand(1,22);
    f(fun,:)=F38(O);
    %f(fun,:)= feval(fhd, O', fun);
end
disp(f)
    


