function [f_iter1,par,Sr]=LowrankSpectra02(f_iter1,par)

lambda = 0.3;

n = par.n;
m = par.m;
kf = par.k;
sizeA = par.sizeA;
jumpA = par.jumpA ;
dd = par.dd;
f_weight1 = par.weight;

f_rank1 = zeros(n,m,kf);
    
for ix = 1:(n - dd)/jumpA
    for iy = 1:(m - dd)/jumpA

        A1{ix,iy} = f_iter1(jumpA*(ix-1)+1:jumpA*(ix-1)+sizeA,jumpA*(iy-1)+1:jumpA*(iy-1)+sizeA,:); 
        A1{ix,iy} = reshape(A1{ix,iy},sizeA.^2,kf);
        [U,S,V] =svd( A1{ix,iy},'econ');
        sig=diag(S)';
        tmp = (sig-lambda.*max(sig))./lambda;
        sigm = max([tmp; zeros(1,length(sig))]);
        tmp = [sig; sigm];
        Sr = diag(min(tmp).*sign(sig));
        a=U*Sr*V';   
        Ar = reshape(a,sizeA,sizeA,kf);
        f_rank1(jumpA*(ix-1)+1:jumpA*(ix-1)+sizeA,jumpA*(iy-1)+1:jumpA*(iy-1)+sizeA,:) = f_rank1(jumpA*(ix-1)+1:jumpA*(ix-1)+sizeA,jumpA*(iy-1)+1:jumpA*(iy-1)+sizeA,:) + Ar;
    end
end
f_rank1 = f_rank1./f_weight1;
f_iter1 =  f_rank1;
f_iter1(f_iter1<1e-10)=0;
f_iter1(isnan(f_iter1) | isinf(f_iter1)) = 0;  




end