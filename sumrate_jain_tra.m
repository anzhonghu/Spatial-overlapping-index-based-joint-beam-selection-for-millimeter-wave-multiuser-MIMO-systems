%%%%%%%%%%%
%%%Algorithm 2
j = 0;
I = ones(K, 1);
I0 = I;
J = zeros(K, 1);
Gj_max = K;
G_flag = zeros(K, Gj_max);
G_in = zeros(Nr, Gj_max);
indexbklbs = zeros(Nr, Gj_max);
indexbklus = zeros(Nr, Gj_max);
%%%%%%%%%%%%%%%%%
% while sum(I) > 0
%     j = j + 1;
%     J = I;
%     G_forbindden = zeros(Nb, 1);
%     while sum(G_flag(:,j))<Nr && sum(J) > 0
%         exe_for_ori;%%%%%%%%%%%%%%%%%%%%%%
%     end
%     J = I0 - I;
%     while sum(G_flag(:,j))<Nr
%         exe_for_ori;%%%%%%%%%%%%%%%%%%%%%%%%5
%     end
% end
% Nbarr = j;
for j = 1 : K
    J = I;
    k = j;
    G_forbindden = zeros(Nb, 1);
    while sum(G_flag(:,j))<Nr
        exe_for_ori_t;%%%%%%%%%%%%%%%%%%%%%%
    end
end
Nbarr = K;
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Algorithm 1
% capa_cu_s = zeros(K, 1);
% for k = 1 : K
%     p = 1;
%     tfr = 1;
%     indexs(k, tfr) = ceil((cos(phib(k, p)) * sin(thetab(k, p))+1)*Nbx*0.5) + (ceil((sin(phib(k, p))+1)*Nby*0.5)-1)*Nbx;
%     indexus(k, tfr) = ceil((cos(phim(k, p)) * sin(thetam(k, p))+1)*Nmx*0.5) + (ceil((sin(phim(k, p))+1)*Nmy*0.5)-1)*Nmx;
%     wtk = Um(:, indexus(k, tfr));
%     ftr = Ub(:, indexs(k, tfr));
%     for n = 1 : N
%         H = H_s(:, Nb*K*(n-1)+1:Nb*K*n);
%         Hkn = H(:, Nb*(k-1)+1:Nb*k);
%         whf = wtk' * Hkn * ftr;
%         capa_cu_s(k, 1) = capa_cu_s(k, 1) + log2(1 + P*abs(whf)^2);
%     end
% end
% sigma_max = sum(capa_cu_s)*T;
% capa_cu_s_s = sort(capa_cu_s, 'ascend');
% sigma_min = sum(capa_cu_s_s(1:Nr))*T;
cal_P;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aaa = sum(PP);
[aaa_max, aaa_max_ind] = max(aaa);
[aaa_min, aaa_min_ind] = min(aaa);
[aaasortre, aaasortord] = sort(aaa, 'ascend');
sigma_max = aaa_max * T;
sigma_min = aaa_min * T;
ksi = 0.1;
miu = 10;
Nni = 20;
ysigmastar = zeros(Nbarr, Nni+1);
rsigmastar = zeros(K, Nni+1);
J_s = zeros(Nni+1, 1);
fea_flag = 0;
for nni = 1 : Nni+1
    sigma = sigma_max - (nni-1) * (sigma_max-sigma_min) / Nni;
    xx = T/Nbarr*ones(Nbarr, 1);
    if aaa * xx < sigma
        for nnbar_in = 1 : Nbarr
            if aaasortre(1,Nbarr-nnbar_in+1)+1e-5 >= sigma/T
            else
                break;
            end
        end
        nnbar_in = nnbar_in - 1;
        while abs(aaa * xx-sigma)>sigma*1e-3
            t_group = sigma*1e-4 / (mean(aaasortre(1,Nbarr-nnbar_in+1:Nbarr)) - mean(aaasortre(1,1:Nbarr-nnbar_in)));
            xx(aaasortord(1,Nbarr-nnbar_in+1:Nbarr), 1) = xx(aaasortord(1,Nbarr-nnbar_in+1:Nbarr), 1) + t_group / nnbar_in * ones(nnbar_in, 1);
            xx(aaasortord(1,1:Nbarr-nnbar_in), 1) = xx(aaasortord(1,1:Nbarr-nnbar_in), 1) - t_group / (Nbarr-nnbar_in) * ones(Nbarr-nnbar_in, 1);
        end
    else
        for nnbar_in = 1 : Nbarr
            if aaasortre(1,nnbar_in)-1e-5 <= sigma/T
            else
                break;
            end
        end
        nnbar_in = nnbar_in - 1;
        while abs(aaa * xx-sigma)>sigma*1e-3
            t_group = sigma*1e-4 / (mean(aaasortre(1,nnbar_in+1:Nbarr)) - mean(aaasortre(1,1:nnbar_in)));
            xx(aaasortord(1,nnbar_in+1:Nbarr), 1) = xx(aaasortord(1,nnbar_in+1:Nbarr), 1) - t_group / (Nbarr-nnbar_in) * ones(Nbarr-nnbar_in, 1);
            xx(aaasortord(1,1:nnbar_in), 1) = xx(aaasortord(1,1:nnbar_in), 1) + t_group / nnbar_in * ones(nnbar_in, 1);
        end
    end
    %solve23;%%%%%%%%%%%%%%%%%5
    %xx = x_sigma_star;
    %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%
    %     xx = 1e-10*ones(Nbarr, 1);
    %     %     if aaa_max >= sigma/T && aaa_min <= sigma/T
    %     bbb = ones(1, 2);
    %     AA = [aaa_max,aaa_min;bbb];
    %     xxso = (AA.'*AA) \ AA.' * [sigma;T];
    %     xx(aaa_max_ind, 1) = xxso(1, 1);
    %     xx(aaa_min_ind, 1) = xxso(2, 1);
    %     else
    %         continue;
    %     end
    %%%%%%%%%%%%%%%%
    x_plus = xx;
    %     flag = 0;
    %     for kkxx = 1 : Nbarr
    %         if xx(kkxx, 1) < 0
    %             flag = 1;
    %             x_plus(kkxx, 1) = 0;
    %         else
    %         end
    %     end
    %     if flag > 0
    %     s = Nbarr / (x_plus.' * P.' * P * x_plus - xx.' * P.' * P * xx);
%     s = 2;
%     while Nbarr / s >= ksi
        solve24t;%%%%%%%%%%%%%%%%%5
        xx = x_sigma_star;
%         s = s * miu;
%     end
    %     else
    %     end
    ysigmastar(:, nni) = xx;
    rsigmastar(:, nni) = PP * ysigmastar(:, nni);
    J_s(nni) = (sum(rsigmastar(:, nni)))^2 / K / (rsigmastar(:, nni)'*rsigmastar(:, nni));
    if nni > 1 && J_s(nni) <= J_s(nni-1)
        fea_flag = 1;
        break;
    else
    end
end
% if fea_flag < 1
%     d = ones(Nbarr, 1) * floor(T/Nbarr);
%     d(1:rem(T, Nbarr), 1) = d(1:rem(T, Nbarr), 1) + ones(rem(T, Nbarr), 1);
% else
d = ysigmastar(:, nni-1);
[~, d_ind] = sort(d, 'ascend');
%%%%%%%%%%%%%%
%%%%ceil d %%%%%%%
d = ceil(d);
d_c = 1;
while sum(d) > T
    d(d_ind(d_c,1),1) = d(d_ind(d_c,1),1) - 1;
    d_c = d_c + 1;
end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
meth = 5;
capa_jain_cal_a;

