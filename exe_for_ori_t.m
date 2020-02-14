G_flag_c = sum(G_flag(:,j));
if G_flag_c > 0
    k = rem(k+1, K);
    if k==0
        k = K;
    else
    end
    bklb_temp = zeros(G_flag_c, 1);
    bklm_temp = zeros(G_flag_c, 1);
    for p = 1 : Nc
        bklb = ceil((cos(phib(k, p)) * sin(thetab(k, p))+1)*Nbx*0.5) + (ceil((sin(phib(k, p))+1)*Nby*0.5)-1)*Nbx;
        bklm = ceil((cos(phim(k, p)) * sin(thetam(k, p))+1)*Nmx*0.5) + (ceil((sin(phim(k, p))+1)*Nmy*0.5)-1)*Nmx;
        flag_temp = 0;
        for gg_temp = 1 : G_flag_c
            if bklb == indexbklbs(gg_temp, j) || G_forbindden(bklb, 1) > 0
                flag_temp = 1;
            else
            end
        end
        if flag_temp<1
            indexbklbs(G_flag_c+1, j) = bklb;
            indexbklus(G_flag_c+1, j) = bklm;
            G_flag(k, j) = 1;
            G_in(G_flag_c+1, j) = k;
            I(k, 1) = 0;
            for pp = 1 : Nc
                if pp == p
                else
                    bklb = ceil((cos(phib(k, pp)) * sin(thetab(k, pp))+1)*Nbx*0.5) + (ceil((sin(phib(k, pp))+1)*Nby*0.5)-1)*Nbx;
                    bklm = ceil((cos(phim(k, pp)) * sin(thetam(k, pp))+1)*Nmx*0.5) + (ceil((sin(phim(k, pp))+1)*Nmy*0.5)-1)*Nmx;
                    if bklm == indexbklus(G_flag_c+1, j) && bklb ~= indexbklus(G_flag_c+1, j)
                        G_forbindden(bklb, 1) = 1;
                    else
                    end
                end
            end
            break;
        else
        end
    end
else
    p = 1;
    bklb = ceil((cos(phib(k, p)) * sin(thetab(k, p))+1)*Nbx*0.5) + (ceil((sin(phib(k, p))+1)*Nby*0.5)-1)*Nbx;
    bklm = ceil((cos(phim(k, p)) * sin(thetam(k, p))+1)*Nmx*0.5) + (ceil((sin(phim(k, p))+1)*Nmy*0.5)-1)*Nmx;
    indexbklbs(G_flag_c+1, j) = bklb;
    indexbklus(G_flag_c+1, j) = bklm;
    G_flag(k, j) = 1;
    G_in(G_flag_c+1, j) = k;
    I(k, 1) = 0;
    for pp = 1 : Nc
        if pp == p
        else
            bklb = ceil((cos(phib(k, pp)) * sin(thetab(k, pp))+1)*Nbx*0.5) + (ceil((sin(phib(k, pp))+1)*Nby*0.5)-1)*Nbx;
            bklm = ceil((cos(phim(k, pp)) * sin(thetam(k, pp))+1)*Nmx*0.5) + (ceil((sin(phim(k, pp))+1)*Nmy*0.5)-1)*Nmx;
            if bklm == indexbklus(G_flag_c+1, j) && bklb ~= indexbklus(G_flag_c+1, j)
                G_forbindden(bklb, 1) = 1;
            else
            end
        end
    end
end