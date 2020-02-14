epsilongsolv = 0.1;
lambdasolv = 10;
lambdasolv_c = 100;
aaa = sum(PP);
bbb = ones(1, Nbarr);
AA = [aaa;bbb];
flag_24 = 0;
index_24 = zeros(Nbarr, 1);
while lambdasolv^2*0.5 > epsilongsolv && abs(lambdasolv-lambdasolv_c)>0.1*lambdasolv
    lambdasolv_c = lambdasolv;
    if flag_24 < 1
        A_temp = 2*PP.'*PP;
        yy = [-2*PP.'*PP*xx; zeros(2, 1)];
        AAA = [A_temp, AA.'; AA, zeros(2, 2)];
        if rank(AAA)<(Nbarr+2)
            deltaxntc = zeros((Nbarr+2), 1);
        else
            deltaxntc = AAA \ yy;
        end
        deltaxnt = deltaxntc(1:Nbarr, 1);
        lambdasolv = sqrt(deltaxnt.' * A_temp * deltaxnt);
        %%%%%%%%%%%%%%
        %backtracking
        alphasolv = 0.1;
        betasolv = 0.5;
        tsolv = 1;
        fxplus = (xx + tsolv * deltaxnt).' * PP.'*PP*(xx + tsolv * deltaxnt) - sum(log(xx + tsolv * deltaxnt));
        fx = xx.' * PP.'*PP*xx - sum(log(xx));
        while fxplus > fx + alphasolv * tsolv * (-yy(1:Nbarr,1)).'*deltaxnt
            tsolv = tsolv * betasolv;
            fxplus = (xx + tsolv * deltaxnt).' * PP.'*PP*(xx + tsolv * deltaxnt) - sum(log(xx + tsolv * deltaxnt));
        end
        %%%%%%%%%%%%%%%%%%%
        xx_temp = xx + tsolv * deltaxnt;
        if min(xx_temp) < 0
            deltaxnt_acc = 0;
            big_c = 0;
            flag_24 = 1;
            inde_str_24 = zeros(Nbarr, 1);
            for nnbar_x = 1 : Nbarr
                if xx_temp(nnbar_x, 1) < 0
                    index_24(nnbar_x, 1) = 1;
                    deltaxnt_acc = deltaxnt_acc - xx_temp(nnbar_x, 1);
                    xx_temp(nnbar_x, 1) = 0;
                else
                    big_c = big_c + 1;
                    inde_str_24(big_c, 1) = nnbar_x;
                end
            end
            xx_temp_ttt = xx_temp;
            for nnbar_x = 1 : Nbarr
                if  xx_temp(nnbar_x, 1) > 0
                    xx_temp(nnbar_x, 1) = xx_temp(nnbar_x, 1) - deltaxnt_acc * xx_temp_ttt(nnbar_x, 1) / sum(xx_temp_ttt);
                else
                end
            end
        else
        end
        xx = xx_temp;
    else
        PP_t = zeros(K, Nbarr-sum(index_24));
        xx_t = zeros(Nbarr-sum(index_24), 1);
        pp_a_i = 1;
        for nnbar_x = 1 : Nbarr
            if index_24(nnbar_x, 1) < 1
                PP_t(:, pp_a_i) = PP(:, nnbar_x);
                xx_t(pp_a_i, 1) = xx(nnbar_x, 1);
                pp_a_i = pp_a_i + 1;
            else
            end
        end
        A_temp = 2*PP_t.'*PP_t;
        yy = [-2*PP_t.'*PP_t*xx_t; zeros(2, 1)];
        aaa_t = sum(PP_t);
        bbb_t = ones(1, Nbarr-sum(index_24));
        AA = [aaa_t;bbb_t];
        AAA = [A_temp, AA.'; AA, zeros(2, 2)];
        if rank(AAA)<(Nbarr-sum(index_24)+2)
            deltaxntc = zeros((Nbarr-sum(index_24)+2), 1);
        else
            deltaxntc = AAA \ yy;
        end
        deltaxnt = deltaxntc(1:Nbarr-sum(index_24), 1);
        lambdasolv = sqrt(deltaxnt.' * A_temp * deltaxnt);
        %%%%%%%%%%%%%%
        %backtracking
        alphasolv = 0.1;
        betasolv = 0.5;
        tsolv = 1;
        fxplus = (xx_t + tsolv * deltaxnt).' * PP_t.'*PP_t*(xx_t + tsolv * deltaxnt);
        fx = xx_t.' * PP_t.'*PP_t*xx_t;
        while fxplus > fx + alphasolv * tsolv * (-yy(1:Nbarr-sum(index_24),1)).'*deltaxnt
            tsolv = tsolv * betasolv;
            fxplus = (xx_t + tsolv * deltaxnt).' * PP_t.'*PP_t*(xx_t + tsolv * deltaxnt);
        end
        %%%%%%%%%%%%%%%%%%%
        xx_temp = xx_t + tsolv * deltaxnt;
        sum_ind = sum(index_24);
        if min(xx_temp) < 0
            deltaxnt_acc = 0;
            big_c = 0;
            inde_str_24_t = zeros(Nbarr, 1);
            for nnbar_x = 1 : Nbarr-sum_ind
                if xx_temp(nnbar_x, 1) < 0
                    index_24(inde_str_24(nnbar_x, 1), 1) = 1;
                    deltaxnt_acc = deltaxnt_acc - xx_temp(nnbar_x, 1);
                    xx_temp(nnbar_x, 1) = 0;
                else
                    big_c = big_c + 1;
                    inde_str_24_t(big_c, 1) = inde_str_24(nnbar_x, 1);
                end
            end
            xx_temp_ttt = xx_temp;
            for nnbar_x = 1 : Nbarr-sum_ind
                if  xx_temp(nnbar_x, 1) > 0
                    xx_temp(nnbar_x, 1) = xx_temp(nnbar_x, 1) - deltaxnt_acc * xx_temp_ttt(nnbar_x, 1) / sum(xx_temp_ttt);
                else
                end
            end
            for nnbar_x = 1 : Nbarr-sum_ind
                xx(inde_str_24(nnbar_x, 1), 1) = xx_temp(nnbar_x,1);
            end
            inde_str_24 = inde_str_24_t;
        else
        end
    end
end
x_sigma_star = xx;