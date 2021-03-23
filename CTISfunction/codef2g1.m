function [ g ] = codef2g1( f,par )
    
    n = par.n;
    m = par.m;
    k = par.k;
    Mu1 = par.Mu1;
    Mu2 = par.Mu2;
    
    cc1 = zeros(n,k-1,k);
    f1 = [f,cc1];   
     
    for i = 1:k
        f1(:,:,i) = circshift(f1(:,:,i),i-1,2);
    end    

    G1 = f1.*Mu1;
    g{1} = sum(f.*Mu2,3);
    g{2} = sum(G1,3);


end

