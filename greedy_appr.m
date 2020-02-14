meth = 1;
capacity_cumu = zeros(K, 1);
capa_cu_s = zeros(K, 1);
indexs = zeros(K, 1);
indexus = zeros(K, 1);

for k = 1 : K
    p = 1;
    tfr = 1;
    indexs(k, tfr) = ceil((cos(phib(k, p)) * sin(thetab(k, p))+1)*Nbx*0.5) + (ceil((sin(phib(k, p))+1)*Nby*0.5)-1)*Nbx;
    indexus(k, tfr) = ceil((cos(phim(k, p)) * sin(thetam(k, p))+1)*Nmx*0.5) + (ceil((sin(phim(k, p))+1)*Nmy*0.5)-1)*Nmx;
    wtk = Um(:, indexus(k, tfr));
    ftr = Ub(:, indexs(k, tfr));
    for n = 1 : N
        H = H_s(:, Nb*K*(n-1)+1:Nb*K*n);
        Hkn = H(:, Nb*(k-1)+1:Nb*k);
        whf = wtk' * Hkn * ftr;
        capa_cu_s(k, 1) = capa_cu_s(k, 1) + log2(1 + P*abs(whf)^2);
    end
end
indextemp = zeros(Nb, 1);
[~, Csel] = sort(capa_cu_s, 'descend');
k = Csel(1, 1);
p = 1;
indexs = zeros(K, 1);
indexus = zeros(K, 1);
indexs(k, 1) = ceil((cos(phib(k, p)) * sin(thetab(k, p))+1)*Nbx*0.5) + (ceil((sin(phib(k, p))+1)*Nby*0.5)-1)*Nbx;
indexus(k, 1) = ceil((cos(phim(k, p)) * sin(thetam(k, p))+1)*Nmx*0.5) + (ceil((sin(phim(k, p))+1)*Nmy*0.5)-1)*Nmx;
k_xx = 1;
indextemp(indexs(k, 1), 1) = 1;
while k_xx < Nr
    capa_greedy = zeros(K, 1);
    indexs_cc = zeros(K, 1);
    indexus_cc = zeros(K, 1);
    for k_x = 2 : K
        k = Csel(k_x, 1);
        if indexs(k, 1) > 0
        else
            for p = 1 : Nc
                indexk = ceil((cos(phib(k, p)) * sin(thetab(k, p))+1)*Nbx*0.5) + (ceil((sin(phib(k, p))+1)*Nby*0.5)-1) * Nbx;
                if indextemp(indexk, 1) > 0
                else
                    calcal_capa_greedy;%%%%%%%%%%%%%%%%%%%%%%%
                    break;
                end
            end
        end
    end
    [~, capa_gree_ind] = max(capa_greedy);
    k_xx = k_xx + 1;
    indexs(capa_gree_ind, 1) = indexs_cc(capa_gree_ind, 1);
    indextemp(indexs(capa_gree_ind, 1), 1) = 1;
    indexus(capa_gree_ind, 1) = indexus_cc(capa_gree_ind, 1);
end
%%%%%%%%%%%%%%%%%%%%%%
% if variable_n > 1
%     capa_greedy_compa = zeros(K, 1);
%     WHF = zeros(Nr, Nr);
%     Hqtn = zeros(Nr, Nb);
%     Ftr = zeros(Nb, Nr);
%     kx = 1;
%     tfr = 1;
%     for kk = 1 : K
%         if indexs_store(kk, tfr) > 0
%             Ftr(:,kx) = Ub(:, indexs_store(kk, tfr));
%             kx = kx + 1;
%         else
%         end
%     end
%     for n = 1 : N
%         H = H_s(:, Nb*K*(n-1)+1:Nb*K*n);
%         kx = 1;
%         for kk = 1 : K
%             if indexs_store(kk, tfr) > 0
%                 wtk = Um(:, indexus_store(kk, tfr));
%                 Hkn = H(:, Nb*(kk-1)+1:Nb*kk);
%                 WHF(kx, :) = wtk' * Hkn * Ftr;
%                 Hqtn(kx, :) = wtk' * Hkn;
%                 kx = kx + 1;
%             else
%             end
%         end
%         X = Hqtn * Ftr * (Ftr' * Ftr)^(-0.5);
%         [~,Dtn,Utn] = svd(X);
%         Gamatn = zeros(Nr, Nr);
%         miu = P / Nr;
%         Gamatnf2 = 0;
%         while abs(Gamatnf2-P)>1e-2*P
%             miu = miu + 1e-3*P;
%             Gamatnf2 = 0;
%             for jgama = 1 : Nr
%                 Gamatn(jgama, jgama) = sqrt(max(miu-1/Dtn(jgama, jgama)^2, 0));
%                 Gamatnf2 = Gamatnf2 + abs(Gamatn(jgama, jgama))^2;
%             end
%         end
%         Ftn = (Ftr' * Ftr)^(-0.5) * Utn * Gamatn;
%         kx = 1;
%         for kk = 1 : K
%             if indexs_store(kk, tfr) > 0
%                 SNRt = abs(WHF(kx, :) * Ftn(:, kx))^2;
%                 interf = 0;
%                 for ktt = 1 : Nr
%                     if ktt == kx
%                     else
%                         interf = interf + abs(WHF(kx, :) * Ftn(:, ktt))^2;
%                     end
%                 end
%                 SNRt = SNRt / (1 + interf);
%                 capa_greedy_compa(kk,1) = capa_greedy_compa(kk,1) + log2(1+SNRt);
%                 kx = kx + 1;
%             else
%             end
%         end
%     end
% else
% end
%%%%%%%%%%%%%%%%%%%%%%%%
WHF = zeros(Nr, Nr);
Hqtn = zeros(Nr, Nb);
Ftr = zeros(Nb, Nr);
kx = 1;
tfr = 1;
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
            capa_s(kk, (variable_n-1)*5+meth) = capa_s(kk, (variable_n-1)*5+meth) + log2(1+SNRt);
            kx = kx + 1;
        else
        end
    end
end
% if variable_n > 1
%     if sum(capa_greedy_compa) > sum(capa_s(:, (variable_n-1)*5+meth))
%         capa_s(:, (variable_n-1)*5+meth) = capa_greedy_compa;
%     else
%         indexs_store = indexs;
%         indexus_store = indexus;
%     end
% else
%     indexs_store = indexs;
%     indexus_store = indexus;
% end
capa_s(:, (variable_n-1)*5+meth) = capa_s(:, (variable_n-1)*5+meth) * T;
capacity(variable_n, meth) = capacity(variable_n, meth) + sum(capa_s(:, (variable_n-1)*5+meth));
jain(variable_n, meth) = jain(variable_n, meth) + (sum(capa_s(:, (variable_n-1)*5+meth)))^2 / K / (capa_s(:, (variable_n-1)*5+meth)'*capa_s(:, (variable_n-1)*5+meth));

