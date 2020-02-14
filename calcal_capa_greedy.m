indexs_cc(k, 1) = indexk;
indexus_cc(k, 1) = ceil((cos(phim(k, p)) * sin(thetam(k, p))+1)*Nmx*0.5) + (ceil((sin(phim(k, p))+1)*Nmy*0.5)-1)*Nmx;
WHF = zeros(k_xx+1, k_xx+1);
Hqtn = zeros(k_xx+1, Nb);
Ftr = zeros(Nb, k_xx+1);
kxcc = 1;
tfr = 1;
for kk = 1 : K
    if indexs(kk, tfr) > 0
        Ftr(:,kxcc) = Ub(:, indexs(kk, tfr));
        kxcc = kxcc + 1;
    else
    end
end
Ftr(:,kxcc) = Ub(:, indexs_cc(k, 1));
for n = 1 : N
    H = H_s(:, Nb*K*(n-1)+1:Nb*K*n);
    kxcc = 1;
    for kk = 1 : K
        if indexs(kk, tfr) > 0
            wtk = Um(:, indexus(kk, tfr));
            Hkn = H(:, Nb*(kk-1)+1:Nb*kk);
            WHF(kxcc, :) = wtk' * Hkn * Ftr;
            Hqtn(kxcc, :) = wtk' * Hkn;
            kxcc = kxcc + 1;
        else
        end
    end
    wtk = Um(:, indexus_cc(k, 1));
    Hkn = H(:, Nb*(k-1)+1:Nb*k);
    WHF(kxcc, :) = wtk' * Hkn * Ftr;
    Hqtn(kxcc, :) = wtk' * Hkn;
    X = Hqtn * Ftr * (Ftr' * Ftr)^(-0.5);
    [~,Dtn,Utn] = svd(X);
    Gamatn = zeros(k_xx+1, k_xx+1);
    miu = P / (k_xx+1);
    Gamatnf2 = 0;
    while abs(Gamatnf2-P)>1e-2*P
        miu = miu + 1e-3*P;
        Gamatnf2 = 0;
        for jgama = 1 : k_xx+1
            Gamatn(jgama, jgama) = sqrt(max(miu-1/Dtn(jgama, jgama)^2, 0));
            Gamatnf2 = Gamatnf2 + abs(Gamatn(jgama, jgama))^2;
        end
    end
    kkcc = 1;
    for kxcc = 1 : K
        if indexs(kxcc, tfr) > 0 || kxcc == k
            Ftn = (Ftr' * Ftr)^(-0.5) * Utn * Gamatn;
            SNRt = abs(WHF(kkcc, :) * Ftn(:, kkcc))^2;
            interf = 0;
            for ktt = 1 : k_xx+1
                if ktt == kkcc
                else
                    interf = interf + abs(WHF(kkcc, :) * Ftn(:, ktt))^2;
                end
            end
            SNRt = SNRt / (1 + interf);
            capa_greedy(k, 1) = capa_greedy(k, 1)  + log2(1+SNRt);
            kkcc = kkcc + 1;
        else
        end
    end
end