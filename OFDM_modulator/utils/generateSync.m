function [sync1 sync2] = generateSync(ifftLen,occupiedCarriers,pilotCarriers)
    carriers1 = [];
    carriers2 = [];
    for k = 1:length(occupiedCarriers)
        if (occupiedCarriers(k) ~= 0)
            if (abs(rem(occupiedCarriers(k),2)) == 1)
                carriers1 = [carriers1 occupiedCarriers(k)];
            end
            carriers2 = [carriers2 occupiedCarriers(k)];        
        end
    end
    for k = 1:length(pilotCarriers)
        if (pilotCarriers(k) ~= 0)        
            if (abs(rem(pilotCarriers(k),2)) == 1)
                carriers1 = [carriers1 pilotCarriers(k)];
            end
            carriers2 = [carriers2 pilotCarriers(k)];
        end
    end
    
    for k = 1:length(carriers1)
        if (carriers1(k) < 0)
            carriers1(k) = carriers1(k) + ifftLen + 1;
        else
            carriers1(k) = carriers1(k) + 1;
        end
    end
    for k = 1:length(carriers2)
        if (carriers2(k) < 0)
            carriers2(k) = carriers2(k) + ifftLen + 1;
        else
            carriers2(k) = carriers2(k) + 1;
        end
    end
    
    sync1 = generateSync1(ifftLen,carriers1);
    sync2 = generateSync2(ifftLen,carriers2);

    dlmwrite('sync1.txt',sync1,',');
    f = fopen('sync1.txt','r');
    sync1 = fread(f,inf)';
    fclose(f);
    f = fopen('sync1.txt','w');
    sync1 = fwrite(f,horzcat('(',sync1(1:end-1),')'),'uint8');
    fclose(f);

    dlmwrite('sync2.txt',sync2,',');
    f = fopen('sync2.txt','r');
    sync2 = fread(f,inf)';
    fclose(f);
    f = fopen('sync2.txt','w');
    sync1 = fwrite(f,horzcat('(',sync2(1:end-1),')'),'uint8');
    fclose(f);
    
end