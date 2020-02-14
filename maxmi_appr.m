meth = 2;
capacity_cumu = zeros(K, 1);
indexs = zeros(K, T);
indexus = zeros(K, T);
for tfr = 1 : T
    indextemp = zeros(Nb, 1);
    [~, Csel] = sort(capacity_cumu, 'ascend');
    k = Csel(1, 1);
    p = 1;
    indexs(k, tfr) = ceil((cos(phib(k, p)) * sin(thetab(k, p))+1)*Nbx*0.5) + (ceil((sin(phib(k, p))+1)*Nby*0.5)-1)*Nbx;
    indexus(k, tfr) = ceil((cos(phim(k, p)) * sin(thetam(k, p))+1)*Nmx*0.5) + (ceil((sin(phim(k, p))+1)*Nmy*0.5)-1)*Nmx;
    k_xx = 1;
    indextemp(indexs(k, tfr), 1) = 1;
    for k_x = 2 : K
        k = Csel(k_x, 1);
        for p = 1 : Nc
            indexk = ceil((cos(phib(k, p)) * sin(thetab(k, p))+1)*Nbx*0.5) + (ceil((sin(phib(k, p))+1)*Nby*0.5)-1) * Nbx;
            if indextemp(indexk, 1) > 0
            else
                k_xx = k_xx + 1;
                indexs(k, tfr) = indexk;
                indextemp(indexs(k, tfr), 1) = 1;
                indexus(k, tfr) = ceil((cos(phim(k, p)) * sin(thetam(k, p))+1)*Nmx*0.5) + (ceil((sin(phim(k, p))+1)*Nmy*0.5)-1)*Nmx;
                break;
            end
        end
        if Nr == k_xx
            break;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    WHF = zeros(Nr, Nr);
    Hqtn = zeros(Nr, Nb);
    Ftr = zeros(Nb, Nr);
    kx = 1;
    for kk = 1 : K
        if indexs(kk, tfr) > 0
            Ftr(:,kx) = Ub(:, indexs(kk, tfr));
            kx = kx + 1;
        else
        end
    end
    for n = 1 : N
        H = H_s(:, Nb*K*(n-1)+1:Nb*K*n);
        kx = 1;
        for kk = 1 : K
            if indexs(kk, tfr) > 0
                wtk = Um(:, indexus(kk, tfr));
                Hkn = H(:, Nb*(kk-1)+1:Nb*kk);
                WHF(kx, :) = wtk' * Hkn * Ftr;
                Hqtn(kx, :) = wtk' * Hkn;
                kx = kx + 1;
            else
            end
        end
        X = Hqtn * Ftr * (Ftr' * Ftr)^(-0.5);
        [~,Dtn,Utn] = svd(X);
        Gamatn = zeros(Nr, Nr);
        miu = P / Nr;
        Gamatnf2 = 0;
        while abs(Gamatnf2-P)>1e-2*P
            miu = miu + 1e-3*P;
            Gamatnf2 = 0;
            for jgama = 1 : Nr
                Gamatn(jgama, jgama) = sqrt(max(miu-1/Dtn(jgama, jgama)^2, 0));
                Gamatnf2 = Gamatnf2 + abs(Gamatn(jgama, jgama))^2;
            end
        end
        Ftn = (Ftr' * Ftr)^(-0.5) * Utn * Gamatn;
        kx = 1;
        for kk = 1 : K
            if indexs(kk, tfr) > 0
                SNRt = abs(WHF(kx, :) * Ftn(:, kx))^2;
                interf = 0;
                for ktt = 1 : Nr
                    if ktt == kx
                    else
                        interf = interf + abs(WHF(kx, :) * Ftn(:, ktt))^2;
                    end
                end
                SNRt = SNRt / (1 + interf);
                capacity_cumu(kk, 1) = capacity_cumu(kk, 1) + log2(1+SNRt);
                capa_s(kk, (variable_n-1)*5+meth) = capa_s(kk, (variable_n-1)*5+meth) + log2(1+SNRt);
                kx = kx + 1;
            else
            end
        end
    end
end
capacity(variable_n, meth) = capacity(variable_n, meth) + sum(capa_s(:, (variable_n-1)*5+meth));
jain(variable_n, meth) = jain(variable_n, meth) + (sum(capa_s(:, (variable_n-1)*5+meth)))^2 / K / (capa_s(:, (variable_n-1)*5+meth)'*capa_s(:, (variable_n-1)*5+meth));


