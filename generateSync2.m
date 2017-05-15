<<<<<<< HEAD
function sync = generateSync2(ifftLen,carriers)
    carriers = sort([carriers]);
    seq = randsrc(1,length(carriers),[-1 1]);
=======
function sync = generateSync2(ifftLen,occupiedCarriers,pilotCarriers)
    carriers = [];
    for k = 1:length(occupiedCarriers)
        if (occupiedCarriers(k) < 0)
            carriers = [carriers occupiedCarriers(k)+ifftLen+1];
        else
            carriers = [carriers occupiedCarriers(k)];
        end
    end
    for k = 1:length(pilotCarriers)
        if (pilotCarriers(k) < 0)
            carriers = [carriers pilotCarriers(k)+ifftLen+1];
        else
            carriers = [carriers pilotCarriers(k)];
        end
    end
    carriers = sort([carriers]);
    seq = randsrc(1,length(carriers));
>>>>>>> caa2e9fcfafb0a4fa2651ef2c22ef7d53706222e
    
    sync = zeros(1,ifftLen);
    sync(carriers) = seq;
    sync(1) = 0+0i;
    sync = ifftshift(sync);
end