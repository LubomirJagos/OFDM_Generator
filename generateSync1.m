<<<<<<< HEAD
function sync = generateSync1(ifftLen,carriers)
    carriers = sort([carriers]);
    seq = randsrc(1,length(carriers),[-1 1]);
=======
function sync = generateSync1(ifftLen,occupiedCarriers,pilotCarriers)
    carriers = [];
    for k = 1:length(occupiedCarriers)
        if (occupiedCarriers(k) < 0)
            occupiedCarriers(k) = occupiedCarriers(k)+ifftLen+1;
        end
        if (rem(occupiedCarriers(k),2) == 0)
            carriers = [carriers occupiedCarriers(k)];
        end
    end
    for k = 1:length(pilotCarriers)
        if (rem(pilotCarriers(k),2) == 0)
            carriers = [carriers pilotCarriers(k)];
        end
    end
    carriers = sort([carriers]);
    seq = randsrc(1,length(carriers));
>>>>>>> caa2e9fcfafb0a4fa2651ef2c22ef7d53706222e
    
    sync = zeros(1,ifftLen);
    sync(carriers) = seq*sqrt(2);
    sync(1) = 0+0i;
    sync = ifftshift(sync);
end