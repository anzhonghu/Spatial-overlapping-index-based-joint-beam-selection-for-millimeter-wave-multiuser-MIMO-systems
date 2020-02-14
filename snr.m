clear all;
close all;
sigma_s = [1e-2;1e-1;1;1e1;1e2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
Nb = 256;
Fb = zeros(Nb, Nb);
for n1 = 1 : Nb
    phin = -pi/2 + pi/2/Nb + (n1-1)*pi/Nb;
    for n2 = 1 : Nb
        Fb(n2, n1) = 1/sqrt(Nb)*exp(1i*pi*(n2-1)*sin(phin));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nm = 16;
sum_rate = zeros(length(sigma_s), 7);
K = 6;
Ncl = 8;
Nray = 10;
sigma_angl = 5 / 180 * pi;
Niten = 1e3;
W = zeros(Nm, K);
am = zeros(Nm, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fm = zeros(Nm, Nm);
for n1 = 1 : Nm
    phin = -pi/2 + pi/2/Nm + (n1-1)*pi/Nm;
    for n2 = 1 : Nm
        Fm(n2, n1) = 1/sqrt(Nm)*exp(1i*pi*(n2-1)*sin(phin));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%
thetams = zeros(Ncl, K);
thetabs = zeros(Ncl, K);
for ite = 1 : Niten
    Hyy = zeros(Nm, K);
    H = zeros(Nb, Nm*K);
    Hxx = zeros(Nb, K);
    ab = zeros(Nb, 1);
    for k = 1 : K
        for lcl = 1 : Ncl
            thetabs(lcl, k) = (rand - 0.5) * pi;
            thetams(lcl, k) = (rand - 0.5) * pi;
            for lray = 1 : Nray
                xxx_temp = rand - 0.5;
                thetabstt = -sigma_angl / sqrt(2) *sign(xxx_temp) * log(1-2*abs(xxx_temp));
                thetabstt = thetabs(lcl, k) + thetabstt;
                xxx_temp = rand - 0.5;
                thetamstt = -sigma_angl / sqrt(2) *sign(xxx_temp) * log(1-2*abs(xxx_temp));
                thetamstt = thetams(lcl, k) + thetamstt;
                for n = 1 : Nm
                    am(n, 1) = 1/sqrt(Nm)*exp(1i*pi*(n-1)*sin(thetamstt));
                end
                alpha_tt = (randn + 1i*randn) / sqrt(2);
                Hyy(:, k) = Hyy(:, k) + alpha_tt * am;
                for n = 1 : Nb
                    ab(n, 1) = 1/sqrt(Nb)*exp(1i*pi*(n-1)*sin(thetabstt));
                end
                H(:, (k-1)*Nm+1:k*Nm) = H(:, (k-1)*Nm+1:k*Nm) + alpha_tt * ab * am';
                Hxx(:, k) = Hxx(:, k) + alpha_tt * ab;
            end
        end
    end
    H = H * sqrt(Nb*Nm/Ncl/Nray);
    FA = zeros(Nb, K);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%Beam mask [6]
    beam_ind = zeros(Nb, 1);
    Hb = Fb' * Hxx;
    for k = 1 : K
        beam_ind = beam_ind + abs(abs(Hb(:, k)).^2 / (Hb(:, k)'*Hb(:, k)));
    end
    beam_ind = sqrt(beam_ind);
    [~, bb_ind] = sort(beam_ind, 'descend');
    bb_ind = bb_ind(1:K, 1);
    Hyyq = abs(Fm' * Hyy);
    mm_ind = zeros(K, 1);
    for k = 1 : K
        [~, mm_ind(k, 1)] = max(Hyyq(:, k));
    end
    for k = 1 : K
        FA(:, k) = Fb(:, bb_ind(k, 1));
    end
    Hq = zeros(K, K);
    for k = 1 : K
        Hq(:, k) = FA' * H(:, (k-1)*Nm+1:k*Nm) * Fm(:, mm_ind(k, 1));
    end
    for nn = 1 : length(sigma_s)
        sigma2 = sigma_s(nn, 1);
        Hqinv = inv(Hq);
        for k = 1 : K
            sum_rate(nn, 1) = sum_rate(nn, 1) + log2(1 + abs(sigma2 / (Hqinv(k, :) * Hqinv(k, :)')));
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%'M-SINR [7]'
    Hbtt = Hb;
    for stepcc1 = 1 : Nb - K
        amphbttt = zeros(Nb - stepcc1 + 1, 1);
        flag_cond_mapp = 0;
        for stepcc2 = 1 : Nb - stepcc1 + 1
            Hbttt = [Hbtt(1:stepcc2-1, :); Hbtt(stepcc2+1:Nb - stepcc1 + 1, :)];
            if cond(Hbttt' * Hbttt) > 1e5
                flag_cond_mapp = 1;
                break;
            else
                amphbttt(stepcc2, 1) = trace(inv(Hbttt' * Hbttt));
            end
        end
        if flag_cond_mapp < 0.5
            [~, stepcc2opt] = min(abs(amphbttt));
        else
            stepcc2opt = stepcc2;
        end
        Hbtt = [Hbtt(1:stepcc2opt-1, :); Hbtt(stepcc2opt+1:Nb - stepcc1 + 1, :)];
    end
    for k = 1 : K
        diff = zeros(Nb, 1);
        for step1 = 1 : Nb
            diff(step1, 1) = norm(Hbtt(k, :) - Hb(step1, :));
        end
        [~, bb_ind(k, 1)] = min(diff);
        FA(:, k) = Fb(:, bb_ind(k, 1));
    end
    Hq = zeros(K, K);
    for k = 1 : K
        Hq(:, k) = FA' * H(:, (k-1)*Nm+1:k*Nm) * Fm(:, mm_ind(k, 1));
    end
    Hqinv = inv(Hq);
    for nn = 1 : length(sigma_s)
        sigma2 = sigma_s(nn, 1);
        for k = 1 : K
            sum_rate(nn, 2) = sum_rate(nn, 2) + log2(1 + abs(sigma2 / (Hqinv(k, :) * Hqinv(k, :)')));
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%'M-capacity [7]'
    Hbtt = Hb;
    for stepcc1 = 1 : Nb - K
        flag_cond_mapp = 0;
        amphbttt = zeros(Nb - stepcc1 + 1, 1);
        for stepcc2 = 1 : Nb - stepcc1 + 1
            Hbttt = [Hbtt(1:stepcc2-1, :); Hbtt(stepcc2+1:Nb - stepcc1 + 1, :)];
            if cond(Hbttt' * Hbttt) > 1e5
                flag_cond_mapp = 1;
                break;
            else
                capa_sl_t = inv(Hbttt' * Hbttt);
                for k = 1 : K
                    amphbttt(stepcc2, 1) = amphbttt(stepcc2, 1) + log2(1 + abs(sigma2 / capa_sl_t(k, k)));
                end
            end
        end
        if flag_cond_mapp < 0.5
            [~, stepcc2opt] = max(abs(amphbttt));
        else
            stepcc2opt = stepcc2;
        end
        Hbtt = [Hbtt(1:stepcc2opt-1, :); Hbtt(stepcc2opt+1:Nb - stepcc1 + 1, :)];
    end
    for k = 1 : K
        diff = zeros(Nb, 1);
        for step1 = 1 : Nb
            diff(step1, 1) = norm(Hbtt(k,:) - Hb(step1, :));
        end
        [~, bb_ind(k, 1)] = min(diff);
        FA(:, k) = Fb(:, bb_ind(k, 1));
    end
    Hq = zeros(K, K);
    for k = 1 : K
        Hq(:, k) = FA' * H(:, (k-1)*Nm+1:k*Nm) * Fm(:, mm_ind(k, 1));
    end
    Hqinv = inv(Hq);
    for nn = 1 : length(sigma_s)
        sigma2 = sigma_s(nn, 1);
        for k = 1 : K
            sum_rate(nn, 3) = sum_rate(nn, 3) + log2(1 + abs(sigma2 / (Hqinv(k, :) * Hqinv(k, :)')));
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%Group [8]
    Hbgg = abs(Hb);
    str_b = zeros(K, 1);
    for k = 1 : K
        [~, str_b(k, 1)] = max(Hbgg(:, k));
    end
    IU_NIU_in = zeros(K, 1);
    for k = 1 : K
        for kk = 1 : K
            if kk ~= k && str_b(k, 1)==str_b(kk, 1)
                IU_NIU_in(k, 1) = 1;
            else
            end
        end
    end
    A = zeros(K-sum(IU_NIU_in), K);
    A_cc = 0;
    for k = 1 : K
        if IU_NIU_in(k, 1)==0
            A_cc = A_cc + 1;
            bb_ind(A_cc, 1) = str_b(k, 1);
            A(A_cc, :) = Hb(bb_ind(A_cc, 1), :);
        else
        end
    end
    Axx = A' * A + 1e-3*eye(K);
    for step0 = A_cc+1 : K
        capa_inc_t = zeros(Nb, 1);
        for step1 = 1 : Nb
            if min(abs(step1-bb_ind)) < 0.5
            else
                Ayy = Axx + Hb(step1, :)' * Hb(step1, :);
                Ayy = inv(Ayy);
                for k = 1 : K
                    capa_inc_t(step1, 1) = capa_inc_t(step1, 1) + log2(1 + abs(sigma2 / Ayy(k, k)));
                end
            end
        end
        [~, m1star] = max(capa_inc_t);
        bb_ind(step0, 1) = m1star;
        Axx = Axx + Hb(m1star, :)' * Hb(m1star, :);
    end
    for k = 1 : K
        FA(:, k) = Fb(:, bb_ind(k, 1));
    end
    Hq = zeros(K, K);
    for k = 1 : K
        Hq(:, k) = FA' * H(:, (k-1)*Nm+1:k*Nm) * Fm(:, mm_ind(k, 1));
    end
    Hqinv = inv(Hq);
    for nn = 1 : length(sigma_s)
        sigma2 = sigma_s(nn, 1);
        for k = 1 : K
            sum_rate(nn, 4) = sum_rate(nn, 4) + log2(1 + abs(sigma2 / (Hqinv(k, :) * Hqinv(k, :)')));
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%Fully-digital
    Hqq = zeros(Nb, K);
    for k = 1 : K
        Hqq(:, k) =  H(:, (k-1)*Nm+1:k*Nm) * Fm(:, mm_ind(k, 1));
    end
    Hq = Hqq'*Hqq;
    Hqinv = inv(Hq);
    for nn = 1 : length(sigma_s)
        sigma2 = sigma_s(nn, 1);
        for k = 1 : K
            sum_rate(nn, 6) = sum_rate(nn, 6) + log2(1 + abs(sigma2 / Hqinv(k,k)));
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%SOI
    bb_ind = zeros(K, 1);
    sklclm = zeros(K*Ncl, Nb);
    beamindic = zeros(K*Ncl, 2);
    Beam_se_indi = zeros(Nb, 1);
    for k = 1 : K
        for lcl = 1 : Ncl
            stepb = floor((thetabs(lcl, k) + pi/2)/pi*Nb) + 1;
            if stepb > Nb
                stepb = Nb;
            else
            end
            sklclm((k-1)*Ncl+lcl, stepb) = 1;
            beamindic((k-1)*Ncl+lcl, 1) = stepb;
            stepb = floor((thetams(lcl, k) + pi/2)/pi*Nm) + 1;
            if stepb > Nm
                stepb = Nm;
            else
            end
            beamindic((k-1)*Ncl+lcl, 2) = stepb;
        end
    end
    for k = 1 : K
        Fa_tt = zeros(Nb, k);
        for kk = 1 : k-1
            Fa_tt(:, kk) = Fb(:, bb_ind(kk, 1));
        end
        Hq_tt = zeros(k, k);
        sumrate_cc = zeros(Ncl, 1);
        for lcl = 1 : Ncl
            over_indik = zeros(K, 1);
            if Beam_se_indi(beamindic((k-1)*Ncl+lcl, 1), 1) > 0
                continue;
            else
            end
            oversi = 0;
            for kk = k+1 : K
                oversi_tt = 0;
                for lllcl = 1 : Ncl
                    oversi_tt = oversi_tt + sklclm((kk-1)*Ncl+lllcl, beamindic((k-1)*Ncl+lcl, 1));
                end
                if oversi_tt > 0
                    over_indik(kk, 1) = 1;
                else
                end
                oversi = oversi + oversi_tt;
            end
            Fa_tt(:, k) = Fb(:, beamindic((k-1)*Ncl+lcl, 1));
            if oversi > 0
                add_user = sum(over_indik);
                Fa_ttxx = zeros(Nb, k+add_user);
                Hq_ttxx = zeros(k+add_user, k+add_user);
                Fa_ttxx(:, 1:k) = Fa_tt;
                add_c = 0;
                mm_indxx = zeros(K, 1);
                for kk = k+1 : K
                    if over_indik(kk, 1) > 0
                        add_c = add_c + 1;
                        for llcl = 1 : Ncl
                            if beamindic((kk-1)*Ncl+llcl, 1) ~= beamindic((k-1)*Ncl+lcl, 1)
                                Fa_ttxx(:, k+add_c) = Fb(:, beamindic((kk-1)*Ncl+llcl, 1));
                                mm_indxx(k+add_c, 1) = beamindic((kk-1)*Ncl+llcl, 2);
                                break;
                            else
                            end
                        end
                    else
                    end
                end
                for kk = 1 : k-1
                    Hq_ttxx(:, kk) = Fa_ttxx' * H(:, (kk-1)*Nm+1:kk*Nm) * Fm(:, mm_ind(kk, 1));
                end
                kk = k;
                Hq_ttxx(:, kk) = Fa_ttxx' * H(:, (kk-1)*Nm+1:kk*Nm) * Fm(:, beamindic((k-1)*Ncl+lcl, 2));
                add_c = 0;
                for kk = k+1 : K
                    if over_indik(kk, 1) > 0
                        add_c = add_c + 1;
                        Hq_ttxx(:, k+add_c) = Fa_ttxx' * H(:, (kk-1)*Nm+1:kk*Nm) * Fm(:, mm_indxx(k+add_c, 1));
                    else
                    end
                end
                Hqinv_tt = inv(Hq_ttxx);
                for kk = 1 : k
                    sumrate_cc(lcl, 1) = sumrate_cc(lcl, 1) + log2(1 + abs(sigma2 / (Hqinv_tt(kk, :) * Hqinv_tt(kk, :)')));
                end
            else
                for kk = 1 : k-1
                    Hq_tt(:, kk) = Fa_tt' * H(:, (kk-1)*Nm+1:kk*Nm) * Fm(:, mm_ind(kk, 1));
                end
                kk = k;
                Hq_tt(:, kk) = Fa_tt' * H(:, (kk-1)*Nm+1:kk*Nm) * Fm(:, beamindic((k-1)*Ncl+lcl, 2));
                Hqinv_tt = inv(Hq_tt);
                for kk = 1 : k
                    sumrate_cc(lcl, 1) = sumrate_cc(lcl, 1) + log2(1 + abs(sigma2 / (Hqinv_tt(kk, :) * Hqinv_tt(kk, :)')));
                end
            end
        end
        [~, sum_cc_inc] = max(sumrate_cc);
        bb_ind(k, 1) = beamindic((k-1)*Ncl+sum_cc_inc, 1);
        mm_ind(k, 1) = beamindic((k-1)*Ncl+sum_cc_inc, 2);
        Beam_se_indi(bb_ind(k, 1), 1) = 1;
    end
    for k = 1 : K
        FA(:, k) = Fb(:, bb_ind(k, 1));
    end
    Hq = zeros(K, K);
    for k = 1 : K
        Hq(:, k) = FA' * H(:, (k-1)*Nm+1:k*Nm) * Fm(:, mm_ind(k, 1));
    end
    Hqinv = inv(Hq);
    for nn = 1 : length(sigma_s)
        sigma2 = sigma_s(nn, 1);
        for k = 1 : K
            sum_rate(nn, 5) = sum_rate(nn, 5) + log2(1 + abs(sigma2 / (Hqinv(k, :) * Hqinv(k, :)')));
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Qupp, Rupp] = qr(Hq);
    for nn = 1 : length(sigma_s)
        sigma2 = sigma_s(nn, 1);
        for k = 1 : K
            sum_rate(nn, 7) = sum_rate(nn, 7) + log2(1 + abs(sigma2 * abs(Rupp(k,k))^2));
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(ite)
end
sum_rate = sum_rate / Niten;
sigma_s = 10 * log10(sigma_s);
%%%%%%%%%%%%%%%%%%%%%%%%%
plot(sigma_s, sum_rate(:, 1), 'k--^','LineWidth',1,'MarkerSize',10)
hold on
plot(sigma_s, sum_rate(:, 2), 'k--s','LineWidth',1,'MarkerSize',10)
plot(sigma_s, sum_rate(:, 3), 'k--x','LineWidth',1,'MarkerSize',10)
plot(sigma_s, sum_rate(:, 4), 'k--d','LineWidth',1,'MarkerSize',10)
plot(sigma_s, sum_rate(:, 5), 'k-*','LineWidth',1,'MarkerSize',10)
plot(sigma_s, sum_rate(:, 7), 'k-o','LineWidth',1,'MarkerSize',10)
plot(sigma_s, sum_rate(:, 6), 'k-','LineWidth',1,'MarkerSize',10)
xlim([min(sigma_s), max(sigma_s)])
%ylim([0, max(max(sum_rate))+10])
le = legend('Beam mask [6]','M-SINR [7]','M-capacity [7]', 'Group [8]', 'SOI','Upper bound', 'Fully-digital', 'Location', 'northwest');
set(le,'Fontname','Times')
set(gca,'XTick',sigma_s)
xlabel('SNR(dB)','Fontname','Times')
ylabel('Sum rate (bps/Hz)','Fontname','Times')
grid on%%%%%%%%%%%%%%%%%%%%%%%%
