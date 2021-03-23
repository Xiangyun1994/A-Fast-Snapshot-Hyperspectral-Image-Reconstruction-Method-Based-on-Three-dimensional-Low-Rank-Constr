function [ f_iter1, Xt, ig, g_div, minq] = acctenfwEMCoded(  f_iter1, opts ,tau,mk,par, Xt )

g = par.data;
iter=par.iter;

data = sptensor(f_iter1);   
maxIter = opts.maxIter;
maxK = opts.maxK;
ex_dims =opts.ex_dims;
l_search = opts.l_search;

A.vals = data.vals;
A.nmodes = ndims(data);
A.size = size(data);
A.subs = data.subs;

D = A.nmodes;
Asize = A.size;
nnzA = length(A.subs);

fac_size = sqrt(Asize)*mk;
dims = 1:D;
dims = setdiff(dims, ex_dims);

nmodes = A.nmodes;
Sigma = cell(nmodes,1);
YSold = cell(nmodes, 1);

U = cell(nmodes, 1);
V = cell(nmodes, 1);

Pi = cell(nmodes, 1);
Pj = cell(nmodes, 1);
sppattern = cell(nmodes, 1);
sppattern_index = cell(nmodes, 1);

for d = dims

    Sigma{d} = zeros(maxK(d), 1);
    U{d} = zeros(Asize(d), maxK(d));
    V{d} = zeros(prod(Asize)/Asize(d), maxK(d));
    
    YSold{d} = zeros(maxK(d), 1);
    
    rdim = d;
    cdims = 1:nmodes;
    cdims = cdims(cdims ~= rdim);
    
    transP = spfold(A.subs, rdim, cdims, Asize);
    Pi{d} = transP(:,1);
    Pj{d} = transP(:,2);
        
    sppattern{d} = sparse(Pi{d}, Pj{d}, ones(size(Pi{d})), Asize(d), prod(Asize)/Asize(d));
    sppattern_index{d} = sparse_order(Pi{d}, Pj{d}, [Asize(d), prod(Asize)/Asize(d)]);
end
    X2.vals = zeros(nnzA,1);
    X2.subs = A.subs;

    X.vals = Xt;
    X.subs = A.subs;
    minq = 10;

for iii=1:iter
           
    data = sptensor(f_iter1);          
    A.vals = data.vals;
    X.vals = zeros(nnzA,1);
    
    
    for t = 1:maxIter

        X2.vals = A.vals - X.vals;    
        X2.nmodes = nmodes;
        X2.size = Asize;
        [j, ~, u, v] = fwsubproblem(X2, dims, sppattern, sppattern_index, fac_size);
        newcomp = tau * fac_size(j)*  spmultic(Pi{j}, Pj{j}, u, v');    
        gamma = linesearch(X, newcomp, X2);
        if gamma == 0
            gamma = 1e-6;
        end
        X.vals  = (1-gamma)*X.vals + gamma * (newcomp);        
    end
    
    data_LR = sptensor(X.subs, X.vals, size(data));
    data_LR = double(data_LR);
    f_iter1 = data_LR;
    Xt = X.vals;
    tic 
    for ij = 1:2
        [ ig ] = codef2g1( f_iter1,par );
        for i = 1:2
            ig{i}(ig{i}<1e-10)=1e-10;
            g_div{i} = g{i}./ig{i}; 
            g_div{i}(isnan(g_div{i}) | isinf(g_div{i})) = 1e-10; 
            g_div{i}(g_div{i}<1e-25)=1e-25;
            g_div{i}(g_div{i}>5e1)=5e1;
        end
        [ f_div ] = codeg2f1( g_div,par );
        f_iter1 = f_iter1.*f_div;
        f_iter1(f_iter1<1e-10)=1e-10; 
    end
       
    fct2 = zeros(210,201,139);
    toc
end
    

   
end
