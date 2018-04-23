function y=CalculateDamp(x,p)

   if x>0
       y=0;
   elseif x<0
       y=-p.d*x;
   end
    
end