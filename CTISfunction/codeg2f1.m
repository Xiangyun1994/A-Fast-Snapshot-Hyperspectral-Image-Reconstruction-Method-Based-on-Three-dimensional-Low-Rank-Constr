function [ f ] = codeg2f1( g,par )

    n = par.n;
    m = par.m;
    k = par.k;
    Mu1 = par.Mu1;
    Mu2 = par.Mu2;

    G1 = repmat(g{1},[1,1,k]); 
    f1 = G1.*Mu2;
    
    G2 = repmat(g{2},[1,1,k]); 
    f2 = G2.*Mu1;
    for ii = 1:k
        f2(:,:,ii) = circshift(f2(:,:,ii),1-ii,2);
    end
    f2 = f2(:,1:m,:);

    
    dd = 0.5;
    f = dd*f1+(1-dd)*f2;

        
end
    
