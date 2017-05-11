function sync = generateSync1(ifftLen,occupiedCarriers,pilotCarriers)
    carriers = [];
    for k = 1:length(occupiedCarriers)
        if (rem(occupiedCarriers(k),2) == 1)
            carriers = [carriers occupiedCarriers(k)];
        end
    end
    for k = 1:length(pilotCarriers)
        if (rem(pilotCarriers(k),2) == 1)
            carriers = [carriers pilotCarriers(k)];
        end
    end
    carriers = sort([carriers]);
    seq = randsrc(1,length(carriers));
    
    sync = zeros(1,ifftLen);
    sync(carriers) = seq*sqrt(2);
    sync(1) = 0+0i;
    sync = ifftshift(sync);
end