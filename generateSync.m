function sync = generateSync(ifftLen,occupiedCarriers,pilotCarriers)
    sync1 = generateSync1(ifftLen, occupiedCarriers, pilotCarriers);
    dlmwrite('sync1.txt',sync1,',');
    sync2 = generateSync2(ifftLen, occupiedCarriers, pilotCarriers);
    dlmwrite('sync2.txt',sync2,',');

    f = fopen('sync1.txt','r');
    content1 = fread(f, inf,'uint8')';
    fclose(f);
    f = fopen('sync1.txt','w');
    fwrite(f,horzcat('[',content1(1:end-1),']'));
    fclose(f);

    f = fopen('sync2.txt','r');
    content2 = fread(f, inf,'uint8')';
    fclose(f);
    f = fopen('sync2.txt','w');
    fwrite(f,horzcat('[',content2(1:end-1),']'));
    fclose(f);
    
    % for gnuradio
    
    
end