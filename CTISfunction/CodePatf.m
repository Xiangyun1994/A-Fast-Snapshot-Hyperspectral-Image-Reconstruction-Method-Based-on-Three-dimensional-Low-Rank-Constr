function [fct]=CodePatf(par)

par.dn = 2;
g = par.data;
iter = par.iter;

[ f_iter1 ] = codeg2f1( g,par );
f_iter1(f_iter1<1e-10)=1e-10; 

for iii=1:iter
        
   [ ig ] = codef2g1( f_iter1,par );  
   for i = 1:2
        ig{i}(ig{i}<1e-10)=1e-10;
        g_div{i} = g{i}./ig{i}; 
        g_div{i}(isnan(g_div{i}) | isinf(g_div{i})) = 1e-10;        
        g_div{i}(g_div{i}>3e0)=3e0;
   end
    [ f_div ] = codeg2f1( g_div,par );
    f_iter1 = f_iter1.*f_div;
    f_iter1(f_iter1<1e-10)=1e-10; 
    tic
    [f_iter1,par,Srt]=LowrankSpectra02(f_iter1,par);     
    toc
  
end

fct=f_iter1;


end