function [ s] = disp0( f )

[m,n,l] = size(f);
lt = (l - mod(l,5))./5 + 1;
d = 0;
while lt>2
    figure,
    for i = 1+d:10+d
        subplot(2,5,i-d);imagesc(f(:,:,i));
        title(i)
    end
    lt = lt - 2;
    d = d + 10;
end
if i~=l
    figure,
    for i = 1+d:l
        subplot(2,5,i-d);imagesc(f(:,:,i));
        title(i)
    end
end

end

