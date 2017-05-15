function sync = generateSync2(ifftLen,carriers)
    carriers = sort([carriers]);
    seq = randsrc(1,length(carriers),[-1 1]);
    
    sync = zeros(1,ifftLen);
    sync(carriers) = seq;
    sync(1) = 0+0i;
    sync = ifftshift(sync);
end