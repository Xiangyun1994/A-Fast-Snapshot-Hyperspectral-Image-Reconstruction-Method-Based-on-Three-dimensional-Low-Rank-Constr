function [fct]=CodeLRf(par)
%para

g = par.data;
err_type = 'AUC';
options.maxIter = 1;
options.shrinkIter = 1;
options.maxK = 200 * ones(3, 1);
% 
tau =2.2e4;
options.ex_dims = [];
options.eps = 5e-5;
options.rho = 0.4;
options.l_search = 1;
options.compute_test_per = 1;
options.accelerated = 0;
options.verbose = 1;
options.err_type = err_type;
options.AUC_size = 1e6;
mk=3.0;

[ f_iter1 ] = codeg2f1( g,par );
f_iter1(f_iter1<1e-10)=1e-10; 

data0 = f_iter1;
data = sptensor(data0);
A.vals = data.vals;
A.nmodes = ndims(data);
A.size = size(data);
A.subs = data.subs;
nnzA = length(A.subs);
Xt = zeros(nnzA,1);

[ fct,Xt,ig,g_div,minq] = acctenfwEMCoded(f_iter1,options,tau,mk,par,Xt);     


end