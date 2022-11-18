function [pl, pl_fact, tmp_pl] = Compute_pl(bp,dot_d3,dot_d3_intr)

warning('off')

x = [ zeros(numel(bp),1) , bp ] ;
y = log(dot_d3) ;

tmp_pl = x \ y ;
pl = -1/tmp_pl(2) ; 

y = log(dot_d3) - log(dot_d3_intr) ;

tmp_pl = x \ y ;
pl_fact = -1/tmp_pl(2) ; 


end

