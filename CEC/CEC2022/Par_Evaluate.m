function o = Par_Evaluate(x,fun)
    fhd=@CEC2022;
    o=feval(fhd, x', fun);
end