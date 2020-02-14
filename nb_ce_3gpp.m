clear all;
close all;
sigma2 = 10.^(60*0.1);%%%%%%%%%%%%40 dbm
rhob = sigma2;%%%%
rhom = sigma2;%%%%%%%%%%%%
Nbs = [32;64;128;128+64;256;];
K = 6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
Nb = Nbs(1, 1);
Fb1 = zeros(Nb, Nb);
FqA11 = zeros(Nb, K);
FqA12 = zeros(Nb, 2*K);
FqA16 = zeros(Nb, 6*K);
FqA130 = zeros(Nb, 30*K);
for n1 = 1 : Nb
    phin = -pi/2 + pi/2/Nb + (n1-1)*pi/Nb;
    for n2 = 1 : Nb
        Fb1(n2, n1) = 1/sqrt(Nb)*exp(1i*pi*(n2-1)*sin(phin));
    end
end
for n1 = 1 : K
    phin = -pi/2 + pi/2/Nb + floor((n1-1)/K*Nb)*pi/Nb;
    for n2 = 1 : Nb
        FqA11(n2, n1) = 1/sqrt(Nb)*exp(1i*pi*(n2-1)*sin(phin));
    end
end
for n1 = 1 : 2*K
    phin = -pi/2 + pi/2/Nb + floor((n1-1)/2/K*Nb)*pi/Nb;
    for n2 = 1 : Nb
        FqA12(n2, n1) = 1/sqrt(Nb)*exp(1i*pi*(n2-1)*sin(phin));
    end
end
for n1 = 1 : 6*K
    phin = -pi/2 + pi/2/Nb + floor((n1-1)/6/K*Nb)*pi/Nb;
    for n2 = 1 : Nb
        FqA16(n2, n1) = 1/sqrt(Nb)*exp(1i*pi*(n2-1)*sin(phin));
    end
end
for n1 = 1 : 30*K
    phin = -pi/2 + pi/2/Nb + floor((n1-1)/30/K*Nb)*pi/Nb;
    for n2 = 1 : Nb
        FqA130(n2, n1) = 1/sqrt(Nb)*exp(1i*pi*(n2-1)*sin(phin));
    end
end
Nb = Nbs(2, 1);
Fb2 = zeros(Nb, Nb);
FqA21 = zeros(Nb, K);
FqA22 = zeros(Nb, 2*K);
FqA26 = zeros(Nb, 6*K);
FqA230 = zeros(Nb, 30*K);
for n1 = 1 : Nb
    phin = -pi/2 + pi/2/Nb + (n1-1)*pi/Nb;
    for n2 = 1 : Nb
        Fb2(n2, n1) = 1/sqrt(Nb)*exp(1i*pi*(n2-1)*sin(phin));
    end
end
for n1 = 1 : K
    phin = -pi/2 + pi/2/Nb + floor((n1-1)/K*Nb)*pi/Nb;
    for n2 = 1 : Nb
        FqA21(n2, n1) = 1/sqrt(Nb)*exp(1i*pi*(n2-1)*sin(phin));
    end
end
for n1 = 1 : 2*K
    phin = -pi/2 + pi/2/Nb + floor((n1-1)/2/K*Nb)*pi/Nb;
    for n2 = 1 : Nb
        FqA22(n2, n1) = 1/sqrt(Nb)*exp(1i*pi*(n2-1)*sin(phin));
    end
end
for n1 = 1 : 6*K
    phin = -pi/2 + pi/2/Nb + floor((n1-1)/6/K*Nb)*pi/Nb;
    for n2 = 1 : Nb
        FqA26(n2, n1) = 1/sqrt(Nb)*exp(1i*pi*(n2-1)*sin(phin));
    end
end
for n1 = 1 : 30*K
    phin = -pi/2 + pi/2/Nb + floor((n1-1)/30/K*Nb)*pi/Nb;
    for n2 = 1 : Nb
        FqA230(n2, n1) = 1/sqrt(Nb)*exp(1i*pi*(n2-1)*sin(phin));
    end
end
Nb = Nbs(3, 1);
Fb3 = zeros(Nb, Nb);
FqA31 = zeros(Nb, K);
FqA32 = zeros(Nb, 2*K);
FqA36 = zeros(Nb, 6*K);
FqA330 = zeros(Nb, 30*K);
for n1 = 1 : Nb
    phin = -pi/2 + pi/2/Nb + (n1-1)*pi/Nb;
    for n2 = 1 : Nb
        Fb3(n2, n1) = 1/sqrt(Nb)*exp(1i*pi*(n2-1)*sin(phin));
    end
end
for n1 = 1 : K
    phin = -pi/2 + pi/2/Nb + floor((n1-1)/K*Nb)*pi/Nb;
    for n2 = 1 : Nb
        FqA31(n2, n1) = 1/sqrt(Nb)*exp(1i*pi*(n2-1)*sin(phin));
    end
end
for n1 = 1 : 2*K
    phin = -pi/2 + pi/2/Nb + floor((n1-1)/2/K*Nb)*pi/Nb;
    for n2 = 1 : Nb
        FqA32(n2, n1) = 1/sqrt(Nb)*exp(1i*pi*(n2-1)*sin(phin));
    end
end
for n1 = 1 : 6*K
    phin = -pi/2 + pi/2/Nb + floor((n1-1)/6/K*Nb)*pi/Nb;
    for n2 = 1 : Nb
        FqA36(n2, n1) = 1/sqrt(Nb)*exp(1i*pi*(n2-1)*sin(phin));
    end
end
for n1 = 1 : 30*K
    phin = -pi/2 + pi/2/Nb + floor((n1-1)/30/K*Nb)*pi/Nb;
    for n2 = 1 : Nb
        FqA330(n2, n1) = 1/sqrt(Nb)*exp(1i*pi*(n2-1)*sin(phin));
    end
end
Nb = Nbs(4, 1);
Fb4 = zeros(Nb, Nb);
FqA41 = zeros(Nb, K);
FqA42 = zeros(Nb, 2*K);
FqA46 = zeros(Nb, 6*K);
FqA430 = zeros(Nb, 30*K);
for n1 = 1 : Nb
    phin = -pi/2 + pi/2/Nb + (n1-1)*pi/Nb;
    for n2 = 1 : Nb
        Fb4(n2, n1) = 1/sqrt(Nb)*exp(1i*pi*(n2-1)*sin(phin));
    end
end
for n1 = 1 : K
    phin = -pi/2 + pi/2/Nb + floor((n1-1)/K*Nb)*pi/Nb;
    for n2 = 1 : Nb
        FqA41(n2, n1) = 1/sqrt(Nb)*exp(1i*pi*(n2-1)*sin(phin));
    end
end
for n1 = 1 : 2*K
    phin = -pi/2 + pi/2/Nb + floor((n1-1)/2/K*Nb)*pi/Nb;
    for n2 = 1 : Nb
        FqA42(n2, n1) = 1/sqrt(Nb)*exp(1i*pi*(n2-1)*sin(phin));
    end
end
for n1 = 1 : 6*K
    phin = -pi/2 + pi/2/Nb + floor((n1-1)/6/K*Nb)*pi/Nb;
    for n2 = 1 : Nb
        FqA46(n2, n1) = 1/sqrt(Nb)*exp(1i*pi*(n2-1)*sin(phin));
    end
end
for n1 = 1 : 30*K
    phin = -pi/2 + pi/2/Nb + floor((n1-1)/30/K*Nb)*pi/Nb;
    for n2 = 1 : Nb
        FqA430(n2, n1) = 1/sqrt(Nb)*exp(1i*pi*(n2-1)*sin(phin));
    end
end
Nb = Nbs(5, 1);
Fb5 = zeros(Nb, Nb);
FqA51 = zeros(Nb, K);
FqA52 = zeros(Nb, 2*K);
FqA56 = zeros(Nb, 6*K);
FqA530 = zeros(Nb, 30*K);
for n1 = 1 : Nb
    phin = -pi/2 + pi/2/Nb + (n1-1)*pi/Nb;
    for n2 = 1 : Nb
        Fb5(n2, n1) = 1/sqrt(Nb)*exp(1i*pi*(n2-1)*sin(phin));
    end
end
for n1 = 1 : K
    phin = -pi/2 + pi/2/Nb + floor((n1-1)/K*Nb)*pi/Nb;
    for n2 = 1 : Nb
        FqA51(n2, n1) = 1/sqrt(Nb)*exp(1i*pi*(n2-1)*sin(phin));
    end
end
for n1 = 1 : 2*K
    phin = -pi/2 + pi/2/Nb + floor((n1-1)/2/K*Nb)*pi/Nb;
    for n2 = 1 : Nb
        FqA52(n2, n1) = 1/sqrt(Nb)*exp(1i*pi*(n2-1)*sin(phin));
    end
end
for n1 = 1 : 6*K
    phin = -pi/2 + pi/2/Nb + floor((n1-1)/6/K*Nb)*pi/Nb;
    for n2 = 1 : Nb
        FqA56(n2, n1) = 1/sqrt(Nb)*exp(1i*pi*(n2-1)*sin(phin));
    end
end
for n1 = 1 : 30*K
    phin = -pi/2 + pi/2/Nb + floor((n1-1)/30/K*Nb)*pi/Nb;
    for n2 = 1 : Nb
        FqA530(n2, n1) = 1/sqrt(Nb)*exp(1i*pi*(n2-1)*sin(phin));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nm = 16;
sum_rate = zeros(length(Nbs), 6);
Ncl = 12;
Nray = 20;
sigma_angl = 5 / 180 * pi;
Niten = 1e3;
W = zeros(Nm, K);
am = zeros(Nm, 1);
Bhe = 10;%%%%%%%%%%
Mhe = 1.5;%%%%%%
fc = 30;%%%%%%%%%%%%%
B = 100;%%%%%%%%%%%%%%%%%%
NF = 9;%%%%%%%%%%%%%%%%%%%
noisep = -174+10*log10(B)+NF;%%%%%%%%%%%%%%
noisep = 10^(0.1*noisep);%%%%%%%%%%%%%%%%%%%%%%
rhob = rhob / noisep;
rhom = rhom / noisep;
sigma2 = sigma2 / noisep;
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
    Nb = Nbs(1, 1);
    H1 = zeros(Nb, Nm*K);
    Hxx1 = zeros(Nb, K);
    ab1 = zeros(Nb, 1);
    Nb = Nbs(2, 1);
    H2 = zeros(Nb, Nm*K);
    Hxx2 = zeros(Nb, K);
    ab2 = zeros(Nb, 1);
    Nb = Nbs(3, 1);
    H3 = zeros(Nb, Nm*K);
    Hxx3 = zeros(Nb, K);
    ab3 = zeros(Nb, 1);
    Nb = Nbs(4, 1);
    H4 = zeros(Nb, Nm*K);
    Hxx4 = zeros(Nb, K);
    ab4 = zeros(Nb, 1);
    Nb = Nbs(5, 1);
    H5 = zeros(Nb, Nm*K);
    Hxx5 = zeros(Nb, K);
    ab5 = zeros(Nb, 1);
    for k = 1 : K
        while 1
            x = rand * 100;
            y = rand * 100;
            if x^2+y^2<=10^4 && x^2+y^2>=100 && abs(atan(y/x))<=pi/3
                break;
            else
            end
        end
        d = sqrt(x^2+y^2+(Bhe-Mhe)^2);
        pathloss = 32.4+21*log10(d)+20*log10(fc);
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
                %%%%%%%%%%%%%%%%%%%%%%
                Nb = Nbs(1, 1);
                for n = 1 : Nb
                    ab1(n, 1) = 1/sqrt(Nb)*exp(1i*pi*(n-1)*sin(thetabstt));
                end
                H1(:, (k-1)*Nm+1:k*Nm) = H1(:, (k-1)*Nm+1:k*Nm) + alpha_tt * ab1 * am';
                Hxx1(:, k) = Hxx1(:, k) + alpha_tt * ab1;
                Nb = Nbs(2, 1);
                ab2 = zeros(Nb, 1);
                for n = 1 : Nb
                    ab2(n, 1) = 1/sqrt(Nb)*exp(1i*pi*(n-1)*sin(thetabstt));
                end
                H2(:, (k-1)*Nm+1:k*Nm) = H2(:, (k-1)*Nm+1:k*Nm) + alpha_tt * ab2 * am';
                Hxx2(:, k) = Hxx2(:, k) + alpha_tt * ab2;
                Nb = Nbs(3, 1);
                ab3 = zeros(Nb, 1);
                for n = 1 : Nb
                    ab3(n, 1) = 1/sqrt(Nb)*exp(1i*pi*(n-1)*sin(thetabstt));
                end
                H3(:, (k-1)*Nm+1:k*Nm) = H3(:, (k-1)*Nm+1:k*Nm) + alpha_tt * ab3 * am';
                Hxx3(:, k) = Hxx3(:, k) + alpha_tt * ab3;
                Nb = Nbs(4, 1);
                ab4 = zeros(Nb, 1);
                for n = 1 : Nb
                    ab4(n, 1) = 1/sqrt(Nb)*exp(1i*pi*(n-1)*sin(thetabstt));
                end
                H4(:, (k-1)*Nm+1:k*Nm) = H4(:, (k-1)*Nm+1:k*Nm) + alpha_tt * ab4 * am';
                Hxx4(:, k) = Hxx4(:, k) + alpha_tt * ab4;
                Nb = Nbs(5, 1);
                ab5 = zeros(Nb, 1);
                for n = 1 : Nb
                    ab5(n, 1) = 1/sqrt(Nb)*exp(1i*pi*(n-1)*sin(thetabstt));
                end
                H5(:, (k-1)*Nm+1:k*Nm) = H5(:, (k-1)*Nm+1:k*Nm) + alpha_tt * ab5 * am';
                Hxx5(:, k) = Hxx5(:, k) + alpha_tt * ab5;
            end
        end
        Hxx1(:, k) = Hxx1(:, k) / 10^(0.1*pathloss);
        Hxx2(:, k) = Hxx2(:, k) / 10^(0.1*pathloss);
        Hxx3(:, k) = Hxx3(:, k) / 10^(0.1*pathloss);
        Hxx4(:, k) = Hxx4(:, k) / 10^(0.1*pathloss);
        Hxx5(:, k) = Hxx5(:, k) / 10^(0.1*pathloss);
        Hyy(:, k) = Hyy(:, k) / 10^(0.1*pathloss);
        H1(:, (k-1)*Nm+1:k*Nm) = H1(:, (k-1)*Nm+1:k*Nm) / 10^(0.1*pathloss);
        H2(:, (k-1)*Nm+1:k*Nm) = H2(:, (k-1)*Nm+1:k*Nm) / 10^(0.1*pathloss);
        H3(:, (k-1)*Nm+1:k*Nm) = H3(:, (k-1)*Nm+1:k*Nm) / 10^(0.1*pathloss);
        H4(:, (k-1)*Nm+1:k*Nm) = H4(:, (k-1)*Nm+1:k*Nm) / 10^(0.1*pathloss);
        H5(:, (k-1)*Nm+1:k*Nm) = H5(:, (k-1)*Nm+1:k*Nm) / 10^(0.1*pathloss);
    end
    for nn = 1 : length(Nbs)
        Nb = Nbs(nn, 1);
        switch nn
            case 1
                Fb = Fb1;
                Hxx = Hxx1;
                H = H1;
                FA1 = FqA11;
                FA2 = FqA12;
                FA6 = FqA16;
                FA30 = FqA130;
            case 2
                Fb = Fb2;
                Hxx = Hxx2;
                H = H2;
                FA1 = FqA21;
                FA2 = FqA22;
                FA6 = FqA26;
                FA30 = FqA230;
            case 3
                Fb = Fb3;
                Hxx = Hxx3;
                H = H3;
                FA1 = FqA31;
                FA2 = FqA32;
                FA6 = FqA36;
                FA30 = FqA330;
            case 4
                Fb = Fb4;
                Hxx = Hxx4;
                H = H4;
                FA1 = FqA41;
                FA2 = FqA42;
                FA6 = FqA46;
                FA30 = FqA430;
            case 5
                Fb = Fb5;
                Hxx = Hxx5;
                H = H5;
                FA1 = FqA51;
                FA2 = FqA52;
                FA6 = FqA56;
                FA30 = FqA530;
            otherwise
        end
        FA = zeros(Nb, K);
        H = H * sqrt(Nb*Nm/Ncl/Nray);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        eta_detect = sqrt(Nm*Nb*rhom*0.5/Ncl) / gamma(0.5);
        sqkk = zeros(K, K);
        mm_ind = zeros(K, 1);
        bb_ind = zeros(K, 1);
        Beam_se_indix = zeros(K, 1);
        Fa1sum = zeros(Nb, 1);
        for k = 1 : K
            Fa1sum = Fa1sum + FA1(:, k);
        end
        RQ = zeros(K, K);
        for k = 1 : K%user
            y_k = sqrt(rhob/K) * H(:, (k-1)*Nm+1:k*Nm)' * Fa1sum + (randn(Nm, 1) + 1i * randn(Nm, 1)) / sqrt(2);
            Hyyq = abs(Fm' * y_k);
            [~, mm_ind(k, 1)] = max(Hyyq);
            r_k = sqrt(rhom) * H(:, (k-1)*Nm+1:k*Nm) * Fm(:, mm_ind(k, 1)) + (randn(Nb, 1) + 1i * randn(Nb, 1)) / sqrt(2);
            rq_k = FA1' * r_k;
            RQ(:, k) = rq_k / sqrt(rhom);
            for kk = 1 : K%beam
                if abs(rq_k(kk, 1)) > eta_detect
                    sqkk(k, kk) = 1;
                else
                end
            end
        end
        for k = 1 : K%user
            Hq_tt = zeros(k, k);
            for kk = 1 : k-1%beam
                for kkk = 1 : k%user
                    Hq_tt(kk, kkk) = RQ(bb_ind(kk, 1), kkk);
                end
            end
            sumrate_cc = zeros(K, 1);
            for kxx = 1 : K%beam
                bb_ind(k, 1) = kxx;
                over_indik = zeros(K, 1);
                if Beam_se_indix(kxx, 1) > 0
                    continue;
                else
                end
                kk = k;%beam
                for kkk = 1 : k%user
                    Hq_tt(kk, kkk) = RQ(bb_ind(kk, 1), kkk);
                end
                oversi = 0;
                for kk = k+1 : K%user
                    oversi = oversi + sqkk(kk, kxx);
                    if sqkk(kk, kxx) > 0
                        over_indik(kk, 1) = 1;
                    else
                    end
                end
                if oversi > 0
                    add_user = sum(over_indik);
                    user_add_ind = [(1:k)';zeros(add_user, 1)];
                    add_c = 0;
                    beam_se_temp = Beam_se_indix;
                    for kk = k+1 : K%user
                        if over_indik(kk, 1) > 0
                            add_c = add_c + 1;
                            user_add_ind(k+add_c, 1) = kk;
                            for kyy = 1 : K
                                if  kyy ~= kxx && beam_se_temp(kyy, 1)<1
                                    bb_ind(kk, 1) = kyy;
                                    beam_se_temp(kyy, 1) = 1;
                                    break;
                                else
                                end
                            end
                        else
                        end
                    end
                    add_c = 0;
                    Hq_ttxx = zeros(k+add_user, k+add_user);
                    for kk = 1 : k+add_user
                        for kkk = 1 : k+add_user
                            Hq_ttxx(kk, kkk) = RQ(bb_ind(user_add_ind(kk, 1), 1), user_add_ind(kkk, 1));
                        end
                    end
                    Hqinv_tt = inv(Hq_ttxx);
                    for kk = 1 : k
                        sumrate_cc(kxx, 1) = sumrate_cc(kxx, 1) + log2(1 + abs(sigma2 / (Hqinv_tt(kk, :) * Hqinv_tt(kk, :)')));
                    end
                else
                    Hqinv_tt = inv(Hq_tt);
                    for kk = 1 : k
                        sumrate_cc(kxx, 1) = sumrate_cc(kxx, 1) + log2(1 + abs(sigma2 / (Hqinv_tt(kk, :) * Hqinv_tt(kk, :)')));
                    end
                end
            end
            [~, sum_cc_inc] = max(sumrate_cc);
            bb_ind(k, 1) = sum_cc_inc;
            Beam_se_indix(bb_ind(k, 1), 1) = 1;
        end
        for k = 1 : K
            FA(:, k) = FA1(:, bb_ind(k, 1));
        end
        Hq = zeros(K, K);
        for k = 1 : K
            Hq(:, k) = FA' * H(:, (k-1)*Nm+1:k*Nm) * Fm(:, mm_ind(k, 1));
        end
        Hq_ce = zeros(K, K);
        for kk = 1 : K
            for kkk = 1 : K
                Hq_ce(kk, kkk) = RQ(bb_ind(kk, 1), kkk);
            end
        end
        Hqceinv = inv(Hq_ce);
        Hqinv = inv(Hq);
        for k = 1 : K
            signalpower = sigma2 * abs(Hqceinv(k, :) * Hq(:, k))^2;
            interfpower = Hqinv(k, :) * Hqinv(k, :)';
            for kk = 1 : K
                if kk == k
                else
                    interfpower = interfpower + sigma2 * abs(Hqceinv(k, :) * Hq(:, kk))^2;
                end
            end
            sum_rate(nn, 1) = sum_rate(nn, 1) + log2(1 + abs(signalpower / interfpower));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Beam_sel_n = 2 * K;
        sqkk = zeros(K, Beam_sel_n);
        mm_ind = zeros(K, 1);
        bb_ind = zeros(K, 1);
        Beam_se_indix = zeros(Beam_sel_n, 1);
        Fa1sum = zeros(Nb, 1);
        for k = 1 : K
            Fa1sum = Fa1sum + FA1(:, k);
        end
        RQ = zeros(Beam_sel_n, K);
        for k = 1 : K%user
            y_k = sqrt(rhob/K) * H(:, (k-1)*Nm+1:k*Nm)' * Fa1sum + (randn(Nm, 1) + 1i * randn(Nm, 1)) / sqrt(2);
            Hyyq = abs(Fm' * y_k);
            [~, mm_ind(k, 1)] = max(Hyyq);
            r_k = sqrt(rhom) * H(:, (k-1)*Nm+1:k*Nm) * Fm(:, mm_ind(k, 1)) + (randn(Nb, 1) + 1i * randn(Nb, 1)) / sqrt(2);
            rq_k = FA2' * r_k;
            RQ(:, k) = rq_k / sqrt(rhom);
            for kk = 1 : Beam_sel_n%beam
                if abs(rq_k(kk, 1)) > eta_detect
                    sqkk(k, kk) = 1;
                else
                end
            end
        end
        for k = 1 : K%user
            Hq_tt = zeros(k, k);
            for kk = 1 : k-1%beam
                for kkk = 1 : k%user
                    Hq_tt(kk, kkk) = RQ(bb_ind(kk, 1), kkk);
                end
            end
            sumrate_cc = zeros(Beam_sel_n, 1);
            for kxx = 1 : Beam_sel_n%beam
                bb_ind(k, 1) = kxx;
                over_indik = zeros(K, 1);
                if Beam_se_indix(kxx, 1) > 0
                    continue;
                else
                end
                kk = k;%beam
                for kkk = 1 : k%user
                    Hq_tt(kk, kkk) = RQ(bb_ind(kk, 1), kkk);
                end
                oversi = 0;
                for kk = k+1 : K%user
                    oversi = oversi + sqkk(kk, kxx);
                    if sqkk(kk, kxx) > 0
                        over_indik(kk, 1) = 1;
                    else
                    end
                end
                if oversi > 0
                    add_user = sum(over_indik);
                    user_add_ind = [(1:k)';zeros(add_user, 1)];
                    add_c = 0;
                    beam_se_temp = Beam_se_indix;
                    for kk = k+1 : K%user
                        if over_indik(kk, 1) > 0
                            add_c = add_c + 1;
                            user_add_ind(k+add_c, 1) = kk;
                            for kyy = 1 : Beam_sel_n
                                if  kyy ~= kxx && beam_se_temp(kyy, 1)<1
                                    bb_ind(kk, 1) = kyy;
                                    beam_se_temp(kyy, 1) = 1;
                                    break;
                                else
                                end
                            end
                        else
                        end
                    end
                    add_c = 0;
                    Hq_ttxx = zeros(k+add_user, k+add_user);
                    for kk = 1 : k+add_user
                        for kkk = 1 : k+add_user
                            Hq_ttxx(kk, kkk) = RQ(bb_ind(user_add_ind(kk, 1), 1), user_add_ind(kkk, 1));
                        end
                    end
                    Hqinv_tt = inv(Hq_ttxx);
                    for kk = 1 : k
                        sumrate_cc(kxx, 1) = sumrate_cc(kxx, 1) + log2(1 + abs(sigma2 / (Hqinv_tt(kk, :) * Hqinv_tt(kk, :)')));
                    end
                else
                    Hqinv_tt = inv(Hq_tt);
                    for kk = 1 : k
                        sumrate_cc(kxx, 1) = sumrate_cc(kxx, 1) + log2(1 + abs(sigma2 / (Hqinv_tt(kk, :) * Hqinv_tt(kk, :)')));
                    end
                end
            end
            [~, sum_cc_inc] = max(sumrate_cc);
            bb_ind(k, 1) = sum_cc_inc;
            Beam_se_indix(bb_ind(k, 1), 1) = 1;
        end
        for k = 1 : K
            FA(:, k) = FA2(:, bb_ind(k, 1));
        end
        Hq = zeros(K, K);
        for k = 1 : K
            Hq(:, k) = FA' * H(:, (k-1)*Nm+1:k*Nm) * Fm(:, mm_ind(k, 1));
        end
        Hq_ce = zeros(K, K);
        for kk = 1 : K
            for kkk = 1 : K
                Hq_ce(kk, kkk) = RQ(bb_ind(kk, 1), kkk);
            end
        end
        Hqceinv = inv(Hq_ce);
        Hqinv = inv(Hq);
        for k = 1 : K
            signalpower = sigma2 * abs(Hqceinv(k, :) * Hq(:, k))^2;
            interfpower = Hqinv(k, :) * Hqinv(k, :)';
            for kk = 1 : K
                if kk == k
                else
                    interfpower = interfpower + sigma2 * abs(Hqceinv(k, :) * Hq(:, kk))^2;
                end
            end
            sum_rate(nn, 4) = sum_rate(nn, 4) + log2(1 + abs(signalpower / interfpower));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Beam_sel_n = 6 * K;
        sqkk = zeros(K, Beam_sel_n);
        mm_ind = zeros(K, 1);
        bb_ind = zeros(K, 1);
        Beam_se_indix = zeros(Beam_sel_n, 1);
        Fa1sum = zeros(Nb, 1);
        for k = 1 : K
            Fa1sum = Fa1sum + FA1(:, k);
        end
        RQ = zeros(Beam_sel_n, K);
        for k = 1 : K%user
            y_k = sqrt(rhob/K) * H(:, (k-1)*Nm+1:k*Nm)' * Fa1sum + (randn(Nm, 1) + 1i * randn(Nm, 1)) / sqrt(2);
            Hyyq = abs(Fm' * y_k);
            [~, mm_ind(k, 1)] = max(Hyyq);
            r_k = sqrt(rhom) * H(:, (k-1)*Nm+1:k*Nm) * Fm(:, mm_ind(k, 1)) + (randn(Nb, 1) + 1i * randn(Nb, 1)) / sqrt(2);
            rq_k = FA6' * r_k;
            RQ(:, k) = rq_k / sqrt(rhom);
            for kk = 1 : Beam_sel_n%beam
                if abs(rq_k(kk, 1)) > eta_detect
                    sqkk(k, kk) = 1;
                else
                end
            end
        end
        for k = 1 : K%user
            Hq_tt = zeros(k, k);
            for kk = 1 : k-1%beam
                for kkk = 1 : k%user
                    Hq_tt(kk, kkk) = RQ(bb_ind(kk, 1), kkk);
                end
            end
            sumrate_cc = zeros(Beam_sel_n, 1);
            for kxx = 1 : Beam_sel_n%beam
                bb_ind(k, 1) = kxx;
                over_indik = zeros(K, 1);
                if Beam_se_indix(kxx, 1) > 0
                    continue;
                else
                end
                kk = k;%beam
                for kkk = 1 : k%user
                    Hq_tt(kk, kkk) = RQ(bb_ind(kk, 1), kkk);
                end
                oversi = 0;
                for kk = k+1 : K%user
                    oversi = oversi + sqkk(kk, kxx);
                    if sqkk(kk, kxx) > 0
                        over_indik(kk, 1) = 1;
                    else
                    end
                end
                if oversi > 0
                    add_user = sum(over_indik);
                    user_add_ind = [(1:k)';zeros(add_user, 1)];
                    add_c = 0;
                    beam_se_temp = Beam_se_indix;
                    for kk = k+1 : K%user
                        if over_indik(kk, 1) > 0
                            add_c = add_c + 1;
                            user_add_ind(k+add_c, 1) = kk;
                            for kyy = 1 : Beam_sel_n
                                if  kyy ~= kxx && beam_se_temp(kyy, 1)<1
                                    bb_ind(kk, 1) = kyy;
                                    beam_se_temp(kyy, 1) = 1;
                                    break;
                                else
                                end
                            end
                        else
                        end
                    end
                    add_c = 0;
                    Hq_ttxx = zeros(k+add_user, k+add_user);
                    for kk = 1 : k+add_user
                        for kkk = 1 : k+add_user
                            Hq_ttxx(kk, kkk) = RQ(bb_ind(user_add_ind(kk, 1), 1), user_add_ind(kkk, 1));
                        end
                    end
                    Hqinv_tt = inv(Hq_ttxx);
                    for kk = 1 : k
                        sumrate_cc(kxx, 1) = sumrate_cc(kxx, 1) + log2(1 + abs(sigma2 / (Hqinv_tt(kk, :) * Hqinv_tt(kk, :)')));
                    end
                else
                    Hqinv_tt = inv(Hq_tt);
                    for kk = 1 : k
                        sumrate_cc(kxx, 1) = sumrate_cc(kxx, 1) + log2(1 + abs(sigma2 / (Hqinv_tt(kk, :) * Hqinv_tt(kk, :)')));
                    end
                end
            end
            [~, sum_cc_inc] = max(sumrate_cc);
            bb_ind(k, 1) = sum_cc_inc;
            Beam_se_indix(bb_ind(k, 1), 1) = 1;
        end
        for k = 1 : K
            FA(:, k) = FA6(:, bb_ind(k, 1));
        end
        Hq = zeros(K, K);
        for k = 1 : K
            Hq(:, k) = FA' * H(:, (k-1)*Nm+1:k*Nm) * Fm(:, mm_ind(k, 1));
        end
        Hq_ce = zeros(K, K);
        for kk = 1 : K
            for kkk = 1 : K
                Hq_ce(kk, kkk) = RQ(bb_ind(kk, 1), kkk);
            end
        end
        Hqceinv = inv(Hq_ce);
        Hqinv = inv(Hq);
        for k = 1 : K
            signalpower = sigma2 * abs(Hqceinv(k, :) * Hq(:, k))^2;
            interfpower = Hqinv(k, :) * Hqinv(k, :)';
            for kk = 1 : K
                if kk == k
                else
                    interfpower = interfpower + sigma2 * abs(Hqceinv(k, :) * Hq(:, kk))^2;
                end
            end
            sum_rate(nn, 5) = sum_rate(nn, 5) + log2(1 + abs(signalpower / interfpower));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Beam_sel_n = 30 * K;
        sqkk = zeros(K, Beam_sel_n);
        mm_ind = zeros(K, 1);
        bb_ind = zeros(K, 1);
        Beam_se_indix = zeros(Beam_sel_n, 1);
        Fa1sum = zeros(Nb, 1);
        for k = 1 : K
            Fa1sum = Fa1sum + FA1(:, k);
        end
        RQ = zeros(Beam_sel_n, K);
        for k = 1 : K%user
            y_k = sqrt(rhob/K) * H(:, (k-1)*Nm+1:k*Nm)' * Fa1sum + (randn(Nm, 1) + 1i * randn(Nm, 1)) / sqrt(2);
            Hyyq = abs(Fm' * y_k);
            [~, mm_ind(k, 1)] = max(Hyyq);
            r_k = sqrt(rhom) * H(:, (k-1)*Nm+1:k*Nm) * Fm(:, mm_ind(k, 1)) + (randn(Nb, 1) + 1i * randn(Nb, 1)) / sqrt(2);
            rq_k = FA30' * r_k;
            RQ(:, k) = rq_k / sqrt(rhom);
            for kk = 1 : Beam_sel_n%beam
                if abs(rq_k(kk, 1)) > eta_detect
                    sqkk(k, kk) = 1;
                else
                end
            end
        end
        for k = 1 : K%user
            Hq_tt = zeros(k, k);
            for kk = 1 : k-1%beam
                for kkk = 1 : k%user
                    Hq_tt(kk, kkk) = RQ(bb_ind(kk, 1), kkk);
                end
            end
            sumrate_cc = zeros(Beam_sel_n, 1);
            for kxx = 1 : Beam_sel_n%beam
                bb_ind(k, 1) = kxx;
                over_indik = zeros(K, 1);
                if Beam_se_indix(kxx, 1) > 0
                    continue;
                else
                end
                kk = k;%beam
                for kkk = 1 : k%user
                    Hq_tt(kk, kkk) = RQ(bb_ind(kk, 1), kkk);
                end
                oversi = 0;
                for kk = k+1 : K%user
                    oversi = oversi + sqkk(kk, kxx);
                    if sqkk(kk, kxx) > 0
                        over_indik(kk, 1) = 1;
                    else
                    end
                end
                if oversi > 0
                    add_user = sum(over_indik);
                    user_add_ind = [(1:k)';zeros(add_user, 1)];
                    add_c = 0;
                    beam_se_temp = Beam_se_indix;
                    for kk = k+1 : K%user
                        if over_indik(kk, 1) > 0
                            add_c = add_c + 1;
                            user_add_ind(k+add_c, 1) = kk;
                            for kyy = 1 : Beam_sel_n
                                if  kyy ~= kxx && beam_se_temp(kyy, 1)<1
                                    bb_ind(kk, 1) = kyy;
                                    beam_se_temp(kyy, 1) = 1;
                                    break;
                                else
                                end
                            end
                        else
                        end
                    end
                    add_c = 0;
                    Hq_ttxx = zeros(k+add_user, k+add_user);
                    for kk = 1 : k+add_user
                        for kkk = 1 : k+add_user
                            Hq_ttxx(kk, kkk) = RQ(bb_ind(user_add_ind(kk, 1), 1), user_add_ind(kkk, 1));
                        end
                    end
                    Hqinv_tt = inv(Hq_ttxx);
                    for kk = 1 : k
                        sumrate_cc(kxx, 1) = sumrate_cc(kxx, 1) + log2(1 + abs(sigma2 / (Hqinv_tt(kk, :) * Hqinv_tt(kk, :)')));
                    end
                else
                    Hqinv_tt = inv(Hq_tt);
                    for kk = 1 : k
                        sumrate_cc(kxx, 1) = sumrate_cc(kxx, 1) + log2(1 + abs(sigma2 / (Hqinv_tt(kk, :) * Hqinv_tt(kk, :)')));
                    end
                end
            end
            [~, sum_cc_inc] = max(sumrate_cc);
            bb_ind(k, 1) = sum_cc_inc;
            Beam_se_indix(bb_ind(k, 1), 1) = 1;
        end
        for k = 1 : K
            FA(:, k) = FA30(:, bb_ind(k, 1));
        end
        Hq = zeros(K, K);
        for k = 1 : K
            Hq(:, k) = FA' * H(:, (k-1)*Nm+1:k*Nm) * Fm(:, mm_ind(k, 1));
        end
        Hq_ce = zeros(K, K);
        for kk = 1 : K
            for kkk = 1 : K
                Hq_ce(kk, kkk) = RQ(bb_ind(kk, 1), kkk);
            end
        end
        Hqceinv = inv(Hq_ce);
        Hqinv = inv(Hq);
        for k = 1 : K
            signalpower = sigma2 * abs(Hqceinv(k, :) * Hq(:, k))^2;
            interfpower = Hqinv(k, :) * Hqinv(k, :)';
            for kk = 1 : K
                if kk == k
                else
                    interfpower = interfpower + sigma2 * abs(Hqceinv(k, :) * Hq(:, kk))^2;
                end
            end
            sum_rate(nn, 6) = sum_rate(nn, 6) + log2(1 + abs(signalpower / interfpower));
        end
        %%%%%%%%%%%%%%%%%%%%%%%
        %%%%Fully-digital
        Hyyq = abs(Fm' * Hyy);
        mm_ind = zeros(K, 1);
        for k = 1 : K
            [~, mm_ind(k, 1)] = max(Hyyq(:, k));
        end
        Hqq = zeros(Nb, K);
        for k = 1 : K
            Hqq(:, k) =  H(:, (k-1)*Nm+1:k*Nm) * Fm(:, mm_ind(k, 1));
        end
        Hq = Hqq'*Hqq;
        Hqinv = inv(Hq);
        for k = 1 : K
            sum_rate(nn, 2) = sum_rate(nn, 2) + log2(1 + abs(sigma2 / Hqinv(k,k)));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%SOI
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
        for k = 1 : K
            sum_rate(nn, 3) = sum_rate(nn, 3) + log2(1 + abs(sigma2 / (Hqinv(k, :) * Hqinv(k, :)')));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp([ite, nn])
    end
end
sum_rate = sum_rate / Niten;
%%%%%%%%%%%%%%%%%%%%%%%%%
plot(Nbs, sum_rate(:, 1), 'k-^','LineWidth',1,'MarkerSize',10)
hold on
plot(Nbs, sum_rate(:, 4), 'k-x','LineWidth',1,'MarkerSize',10)
plot(Nbs, sum_rate(:, 5), 'k-s','LineWidth',1,'MarkerSize',10)
plot(Nbs, sum_rate(:, 6), 'k-o','LineWidth',1,'MarkerSize',10)
plot(Nbs, sum_rate(:, 3), 'k-*','LineWidth',1,'MarkerSize',10)
plot(Nbs, sum_rate(:, 2), 'k-','LineWidth',1,'MarkerSize',10)
xlim([min(Nbs), max(Nbs)])
ylim([min(min(sum_rate))-5, max(max(sum_rate))+30])
le = legend('SOI-CE(N_p=K)','SOI-CE(N_p=2K)','SOI-CE(N_p=6K)','SOI-CE(N_p=30K)','SOI-full CSI', 'Fully-digital', 'Location', 'northwest');
set(le,'Fontname','Times')
set(gca,'XTick',Nbs)
xlabel('Number of BS antennas','Fontname','Times')
ylabel('Sum rate (bps/Hz)','Fontname','Times')
grid on%%%%%%%%%%%%%%%%%%%%%%%%
