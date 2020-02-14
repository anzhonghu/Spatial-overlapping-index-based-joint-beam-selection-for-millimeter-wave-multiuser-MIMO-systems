clear;
close all;
Nbx = 16;
Nby = 4;
Nmx = 4;
Nmy = 2;
Nb = Nbx * Nby;
Nm = Nmx * Nmy;
Nr = 8;
K = 30;%30
Nc = 6;
Nsc = 5;
D = 8;%128
N = 32;%512
T = 50;%50
Nite = 1e4;
variable_s = [30; 50; 100; 150; 200];
capacity = zeros(length(variable_s),  5);
jain = zeros(length(variable_s),  5);
Ub = zeros(Nb, Nb);
for nx = 1 : Nbx
    for ny = 1 : Nby
        angle(1, 2) = (-1+2*ny/Nby);%el
        angle(1, 1) = (-1+2*nx/Nbx);%az
        n = (ny - 1) * Nbx + nx;
        for mx = 0 : Nbx-1
            for my = 0 : Nby-1
                m = my * Nbx + 1 + mx;
                Ub(m, n) = exp(1i * pi * (mx * angle(1, 1) + my * angle(1, 2))) / sqrt(Nb);
            end
        end
    end
end
Um = zeros(Nm, Nm);
for nx = 1 : Nmx
    for ny = 1 : Nmy
        angle(1, 2) = (-1+2*ny/Nmy);%el
        angle(1, 1) = (-1+2*nx/Nmx);%az
        n = (ny - 1) * Nmx + nx;
        for mx = 0 : Nmx-1
            for my = 0 : Nmy-1
                m = my * Nmx + 1 + mx;
                Um(m, n) = exp(1i * pi * (mx * angle(1, 1) + my * angle(1, 2))) / sqrt(Nm);
            end
        end
    end
end
for ii = 1 : Nite
    capa_s = zeros(K, length(variable_s)*5);
    H = zeros(Nm, Nb*K);
    H_s = zeros(Nm, Nb*K*N);
    thetab = zeros(K, Nc);
    phib = zeros(K, Nc);
    thetam = zeros(K, Nc);
    phim = zeros(K, Nc);
    thetabts = zeros(K, Nc*Nsc);
    phibts = zeros(K, Nc*Nsc);
    thetamts = zeros(K, Nc*Nsc);
    phimts = zeros(K, Nc*Nsc);
    alphas = zeros(K, Nc*Nsc);
    tau_taus = zeros(K, Nc*Nsc);%/Ts
    for k = 1 : K
        for p = 1 : Nc
            thetab(k, p) = (rand - 0.5) * 2 * pi / 6;
            phib(k, p) = (rand - 0.5) * 2 * pi / 2;
            thetam(k, p) = (rand - 0.5) * 2 * pi * 2;
            phim(k, p) = (rand - 0.5) * 2 * pi / 2;
            for pp = 1 : Nsc
                thetabts(k, (p-1)*Nsc+pp) = thetab(k, p) + (rand - 0.5) * 2 * pi * 0.05;
                phibts(k, (p-1)*Nsc+pp) = phib(k, p) + (rand - 0.5) * 2 * pi * 0.05;
                thetamts(k, (p-1)*Nsc+pp) = thetam(k, p) + (rand - 0.5) * 2 * pi * 0.05;
                phimts(k, (p-1)*Nsc+pp) = phim(k, p) + (rand - 0.5) * 2 * pi * 0.05;
                alphas(k, (p-1)*Nsc+pp) = (randn + 1i * randn) / sqrt(2) * sqrt(Nb*Nm/Nc/Nsc);
                tau_taus(k, (p-1)*Nsc+pp) = rand * D + rand * D * 0.05;
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%
    for n = 1 : N
        for k = 1 : K
            for p = 1 : Nc
                for pp = 1 : Nsc
                    thetabt = thetabts(k, (p-1)*Nsc+pp);
                    phibt = phibts(k, (p-1)*Nsc+pp);
                    thetamt = thetamts(k, (p-1)*Nsc+pp);
                    phimt = phimts(k, (p-1)*Nsc+pp);
                    ab = zeros(Nb, 1);
                    for nx = 0 : Nbx-1
                        for ny = 0 : Nby-1
                            nn = ny * Nbx + 1 + nx;
                            ab(nn, 1) = 1 / sqrt(Nb) * exp(1i * pi * (nx * cos(phibt) * sin(thetabt) + ny * sin(phibt)));
                        end
                    end
                    am = zeros(Nm, 1);
                    for nx = 0 : Nmx-1
                        for ny = 0 : Nmy-1
                            nn = ny * Nmx+ 1 + nx;
                            am(nn, 1) = exp(1i * pi * (nx * cos(phimt) * sin(thetamt) + ny * sin(phimt)));
                        end
                    end
                    beta = 0;
                    for d = 0 : D-1
                        t = d - tau_taus(k, (p-1)*Nsc+pp);
                        if abs(t-0.5)<1e-3 || abs(t+0.5)<1e-3
                            pt = pi / 4 * sin(pi * 0.5) / (pi * 0.5);
                        else
                            pt = sin(pi * t) / (pi * t) * cos(pi * t) / (1-4*t^2);
                        end
                        beta = beta + pt * exp(-1i * 2 * pi * n / N * d);
                    end
                    H(:, Nb*(k-1)+1:Nb*k) = H(:, Nb*(k-1)+1:Nb*k) + alphas(k, (p-1)*Nsc+pp) * beta * am * ab';
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        H_s(:, Nb*K*(n-1)+1:Nb*K*n) = H;
    end
    snr = 1e2;
    P = snr;
    for variable_n = 1 : length(variable_s)
        T = variable_s(variable_n, 1);
        %%Schedule%%%%%%%%%%%%%%%%%%%%
        greedy_appr;%1,WF+greedy
        maxmi_appr;%2,WF+maxmin
        pf_appr;%3,WF+pf
        rr_appr;%4,WF+rr
        sumrate_jain_tra;%%5,Schedule with the proposed approach
        disp([ii, variable_n])
    end
end
capacity = capacity / Nite;
jain = jain / Nite;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
variable_ss = variable_s;
subplot(1,2,1);
plot(variable_ss, capacity(:, 1), 'k--s','LineWidth',1,'MarkerSize',10)
hold on
plot(variable_ss, capacity(:, 2), 'k--*','LineWidth',1,'MarkerSize',10)
plot(variable_ss, capacity(:, 3), 'k--o','LineWidth',1,'MarkerSize',10)
plot(variable_ss, capacity(:, 4), 'k--<','LineWidth',1,'MarkerSize',10)
plot(variable_ss, capacity(:, 5), 'k-^','LineWidth',1,'MarkerSize',10)

xlim([min(variable_ss), max(variable_ss)])
le = legend('Greedy','Max-min','PF','RR','Proposed', 'Location', 'northwest');
set(le,'Fontname','Times')
set(gca,'XTick',variable_ss)
xlabel('Channel coherence interval T','Fontname','Times')
ylabel('Capacity (bps)','Fontname','Times')
grid on
subplot(1,2,2);
plot(variable_ss, jain(:, 1), 'k--s','LineWidth',1,'MarkerSize',10)
hold on
plot(variable_ss, jain(:, 2), 'k--*','LineWidth',1,'MarkerSize',10)
plot(variable_ss, jain(:, 3), 'k--o','LineWidth',1,'MarkerSize',10)
plot(variable_ss, jain(:, 4), 'k--<','LineWidth',1,'MarkerSize',10)
plot(variable_ss, jain(:, 5), 'k-^','LineWidth',1,'MarkerSize',10)
xlim([min(variable_ss), max(variable_ss)])
%le = legend('Greedy','Max-min','PF','RR','Proposed', 'Location', 'northwest');
set(le,'Fontname','Times')
set(gca,'XTick',variable_ss)
xlabel('Channel coherence interval T','Fontname','Times')
ylabel('Jain''s fairness index','Fontname','Times')
grid on