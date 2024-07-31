close all;
clear all;
% Physical Params
perms = [5e-15];
% Reservoir-Sediment Interface Morphology, uniform (u) or heterogeneous (h)
% morph = 'u';
% if morph == 'u'
%     paths = [{'/Users/kmagno/Documents/Final_Sherlock_Data/R50/flat_np21/'},{'/Users/kmagno/Documents/Final_Sherlock_Data/R50/flat_np24/'},...
%         {'/Users/kmagno/Documents/Final_Sherlock_Data/R50/flat_np27/'},{'/Users/kmagno/Documents/Final_Sherlock_Data/R50/flat_np30/'},{'/Users/kmagno/Documents/Final_Sherlock_Data/R50/flat_np33/'}];
% elseif morph == 'h'
%     paths = [{'/Users/kmagno/Documents/Final_Sherlock_Data/R50/np21_R50/'},{'/Users/kmagno/Documents/Final_Sherlock_Data/R50/np24_R50/'},...
%         {'/Users/kmagno/Documents/Final_Sherlock_Data/R50/np27_R50/'},{'/Users/kmagno/Documents/Final_Sherlock_Data/R50/np30_R50/'},{'/Users/kmagno/Documents/Final_Sherlock_Data/R50/np33_R50/'}];
% end
paths =[{'/Users/kmagno/Documents/Final_Sherlock_Data/np30_flat_Phi/'}];

kk=1;
perm_actual = perms(kk);
mu_f = 5e-5;
eta_phi = 5e15;
delta_rhog = 1000*9.81;
tau_c = ((eta_phi*mu_f/perm_actual)^0.5)/delta_rhog;
delta_c = sqrt(perm_actual*eta_phi/mu_f);
tau_years = tau_c/(3600*24*365);
flux_scale = (delta_c/tau_years);
final_area = 0;
final_flux = 0;
for jj = 1:length(paths)
    clearvars -except delta_c tau_years flux_scale paths jj kk perm_actual data_count vol_flux error_volflux
    res = 256;
    %res = resis(jj);
    ny = res*2-1;
    nx = res-1;
    colors_256_3 = [191 211 230; 158 188 218; 140 150 198; 136 86 167; 129 15 124];
    colors_256_3_darker = ([191 211 230; 158 188 218; 140 150 198; 136 86 167; 129 15 124]).*0.90;
    colors_3= colors_256_3./256;
    colors_3_darker = colors_256_3_darker./256;
    FS = 18;
    import_path = string(paths{jj});
    % Load Time and Flux Data and Set Plot Colors
    color_num = jj;
    time_cell = dir(strcat(import_path,'times*.csv'));
    time_fnm = {time_cell(:).name};
    time_fnm_oi = string(time_fnm{1});
    times_t = readmatrix(strcat(import_path,time_fnm_oi));
    A = dir(strcat(import_path,'qDys*.csv'));
    names = {A(:).name};
    B=natsort(names);

    C = dir(strcat(import_path,'Phi*.csv'));
    names = {C(:).name};
    P = natsort(names);

    count = 0;
    ly = 70;
    lx = ly/2;
    dy = ly/(ny) ;
    delta_ny = (200/delta_c)*(ny/ly);
    for ii = 1:length(B)
        T = readmatrix(strcat(import_path,string(B{ii})));
        T = rot90(T);

        Phi = readmatrix(strcat(import_path,string(P{ii})));
        Phi = rot90(Phi);
        if ii == 1
            [rs,~] = find(mean(T,2) > 1.0);
            horz_slice = round(rs(1)-delta_ny);
            disp('horz_slice:')
            disp(horz_slice/ny)
        end
        D = T(horz_slice,:);
        P_D = Phi(horz_slice,:);
        slices(ii,:) = D;
        Phi_slices(ii,:) = P_D;
    end
    [new_dist,mu,sig] = clt(slices(:));
    slices(slices < 1.0)=0;
    % Find where the pockmarks form
    K = find(sum(slices,2) ~=0);
    J = find(sum(Phi_slices,2) ~=0);
    % First row that's non-zero
    start_oi = K(1);
    start_oi_phi = K(1);
    rows_oi = slices(start_oi:end,:);
    phi_oi = Phi_slices(start_oi_phi:end,:);
    times_oi = times_t(start_oi:end);
    time_period = times_oi.*tau_years;
    locs_list = [];
    full_pic = zeros(2,2);
    full_rad = zeros(2,2);
    full_phi = zeros(2,2);
    flux_std = zeros(2,2);
    area_std = zeros(2,2);
    count_loc = 0;

    for rr =1:size(rows_oi,1)
        single_row = rows_oi(rr,:);
        single_phi = phi_oi(rr,:);
        [pks,locs,w,proms] = findpeaks(single_row,'WidthReference','halfheight');
  
        for nn = 1:length(locs)
            %Determine area of each peak for each loc
            promw_test = proms(nn)/w(nn);

            if promw_test > 1.0
                start_pt = floor(locs(nn)-(w(nn)));
                end_pt = floor(locs(nn)+(w(nn)));
                if start_pt <= 0
                    start_pt = 1;
                end
                if end_pt > nx
                    end_pt = nx;
                end
                int_row = single_row(start_pt:end_pt);
                int_phi = single_phi(start_pt:end_pt);
                int_flux = sum(int_row);
                avg_phi  = sum(int_phi);
                pk_radius = w(nn);
                limit_flux = 0;
                if pk_radius == 0
                    continue
                else
                    if ~any(locs(nn) == locs_list) && int_flux > limit_flux
                        if (any(locs(nn) == locs(1:nn-1)+1) || any(locs(nn) == locs(1:nn-1)-1))
                            continue
                        else
                            count_loc = count_loc + 1;
                            locs_list(count_loc) = locs(nn);
                            full_pic(1,count_loc) = int_flux;
                            full_phi(1,count_loc) = avg_phi;
                            flux_std(1,count_loc) = std(int_row);
                            full_rad(1,count_loc) = pk_radius;
                            area_std(1,count_loc) = 1/2;
                            %(lx/(2*nx));
                        end
                    elseif any(locs(nn) == locs_list) && int_flux > limit_flux
                        count_loc_copy = find(locs(nn)==locs_list);
                        add_row = size(full_pic(:,count_loc_copy),1);
                        full_pic(add_row+1,count_loc_copy) = int_flux;
                        full_phi(add_row+1,count_loc_copy) = avg_phi;
                        flux_std(add_row+1,count_loc_copy) = std(int_row);
                        full_rad(add_row+1,count_loc_copy) = pk_radius;
                        area_std(add_row+1,count_loc_copy) = 1/2;
                        %(lx/(2*nx));
                    end
                end
            else
                continue
            end
        end
    end
    % Convert to rad to area
    km_area = (pi*(full_rad*(lx/nx)*delta_c).^2)*1e-6; %km^2
    km3_vol_flux = full_pic.*km_area;
    volflux_tevolv = cumsum(sum(km3_vol_flux,2));
    count_std = 0;
    for gg = 11:size(full_pic,2)
        nz_avg_flux(gg) = mean(nonzeros(full_pic(:,gg)));
        nz_avg_radius(gg) = mean(nonzeros(full_rad(:,gg)));
        nz_avg_phi(gg) = mean(nonzeros(full_phi(:,gg)));
        nz_var_rad(gg) = 0.5./sqrt(length(nonzeros(area_std(:,gg))));
        nz_var_flux(gg) = sqrt(sum((nonzeros(flux_std(:,gg))).^2)/length(nonzeros(flux_std(:,gg))));
    end

    SEM_rad= nz_var_rad;
    SEM_flux= nz_var_flux;
    fin_flux= nz_avg_flux;
    fin_rad = nz_avg_radius;
    fin_phi = nz_avg_phi;

    fin_flux_kmpyr = fin_flux*flux_scale*1e-3; %km/yr
    fin_Sigmaflux = SEM_flux*flux_scale*1e-3; %km/yr

    fin_SigmaArea = (pi*(SEM_rad*(lx/nx)*delta_c).^2)*1e-6;
    fin_area_km2 = (pi*(fin_rad*(lx/nx)*delta_c).^2)*1e-6; %km^2

    figure(kk),set(gcf,'Position',[45 99 1066 663]);
    err1_color = [254 224 182]./255;
    err2_color = [253 184 99]./255;

    scatter(fin_area_km2,fin_flux_kmpyr,150*ones([1 length(fin_flux_kmpyr)]),'filled','LineWidth',60,'MarkerFaceColor',colors_3(jj,:),'MarkerEdgeColor',colors_3_darker(jj,:),'LineWidth',2); grid on;

    hold on;
    %errorbar(fin_area_km2,fin_flux_kmpyr,SEM_rad,'horizontal', 'LineWidth',1,'LineStyle', 'none','Color',err2_color); hold on;
    leg= legend('$\textbf{\boldmath{n$_{k}$=2.1}}$',...
        '$\textbf{\boldmath{n$_{k}$=2.4}}$',...
        '$\textbf{\boldmath{n$_{k}$=2.7}}$',...
        '$\textbf{\boldmath{n$_{k}$=3.0}}$','$\textbf{\boldmath{n$_{k}$=3.3}}$'); ...
    set(leg,'interpreter','latex','Fontsize',FS-3,'FontWeight','bold');
    ax2 = gca;
    ax2.FontSize = FS;
    ax2.TickLabelInterpreter = 'latex';
    xlabel('$\textbf{Pockmark Area  \boldmath{[km$^{2}$]}}$','interpreter','latex','Fontsize',FS);
    ylabel('$\textbf{Fluid Release Rate \boldmath{[km/yr]}}$','interpreter','latex','Fontsize',FS,'FontWeight','bold');
    hold on;

    % Sort sizes
    count_p1 = 0;
    count_p2 = 0;
    count_p3 = 0;
    count_p4 = 0;
    count_p5 = 0;
    count_p6 = 0;
    count_p7 = 0;
    count_p8 = 0;
    count_p9 = 0;
    count_p10 = 0;
    count_p11 = 0;
    count_p12 = 0;
    count_p13 = 0;
    p1_flux = 0;
    p2_flux = 0;
    p3_flux = 0;
    p4_flux = 0;
    p5_flux = 0;
    p6_flux = 0;
    p7_flux = 0;
    p8_flux = 0;
    p9_flux = 0;
    p10_flux = 0;
    p11_flux = 0;
    p12_flux = 0;
    p13_flux = 0;

    p1_area = 0;
    p2_area = 0;
    p3_area = 0;
    p4_area = 0;
    p5_area = 0;
    p6_area = 0;
    p7_area = 0;
    p8_area = 0;
    p9_area = 0;
    p10_area = 0;
    p11_area = 0;
    p12_area = 0;
    p13_area = 0;
    
    p1_phi = 0;
    p2_phi = 0;
    p3_phi = 0;
    p4_phi = 0;
    p5_phi = 0;
    p6_phi = 0;
    p7_phi = 0;
    p8_phi = 0;
    p9_phi = 0;
    p10_phi = 0;
    p11_phi = 0;
    p12_phi = 0;
    p13_phi = 0;

    p1_flux_var = 0;
    p2_flux_var = 0;
    p3_flux_var = 0;
    p4_flux_var = 0;
    p5_flux_var = 0;
    p6_flux_var = 0;
    p7_flux_var = 0;
    p8_flux_var = 0;
    p9_flux_var = 0;
    p10_flux_var = 0;
    p11_flux_var = 0;
    p12_flux_var = 0;
    p13_flux_var = 0;


    p1_area_var = 0;
    p2_area_var = 0;
    p3_area_var = 0;
    p4_area_var = 0;
    p5_area_var = 0;
    p6_area_var = 0;
    p7_area_var = 0;
    p8_area_var = 0;
    p9_area_var = 0;
    p10_area_var = 0;
    p11_area_var = 0;
    p12_area_var = 0;
    p13_area_var = 0;

    for ff = 1:size(fin_area_km2,2)
        area_oi = fin_area_km2(ff);
        if area_oi < 0.05
            continue
        elseif area_oi > 1.45
            continue
        elseif 0.05 < round(area_oi,2) && round(area_oi,2) <= 0.15
            count_p1 = count_p1 + 1;
            p1_area_ind(count_p1) = ff;
            p1_flux(count_p1) = fin_flux_kmpyr(p1_area_ind(count_p1));
            p1_area(count_p1) = fin_area_km2(p1_area_ind(count_p1));
            p1_phi(count_p1) = fin_phi(p1_area_ind(count_p1));
            p1_flux_var(count_p1) = fin_Sigmaflux(p1_area_ind(count_p1));
            p1_area_var(count_p1) = fin_SigmaArea(p1_area_ind(count_p1));
            
        elseif 0.15 < round(area_oi,2) && round(area_oi,2) <= 0.25
            count_p2 = count_p2 + 1;
            p2_area_ind(count_p2) = ff;
            p2_flux(count_p2) = fin_flux_kmpyr(p2_area_ind(count_p2));
            p2_area(count_p2) = fin_area_km2(p2_area_ind(count_p2));
            p2_phi(count_p2) = fin_phi(p2_area_ind(count_p2));
            p2_flux_var(count_p2) = fin_Sigmaflux(p2_area_ind(count_p2));
            p2_area_var(count_p2) = fin_SigmaArea(p2_area_ind(count_p2));

        elseif 0.25 < round(area_oi,2) && round(area_oi,2)<= 0.35
            count_p3 = count_p3 + 1;
            p3_area_ind(count_p3) = ff;
            p3_flux(count_p3) = fin_flux_kmpyr(p3_area_ind(count_p3));
            p3_area(count_p3) = fin_area_km2(p3_area_ind(count_p3));
            p3_phi(count_p3) = fin_phi(p3_area_ind(count_p3));
            p3_flux_var(count_p3) = fin_Sigmaflux(p3_area_ind(count_p3));
            p3_area_var(count_p3) = fin_SigmaArea(p3_area_ind(count_p3));
        elseif 0.35 < round(area_oi,2) && round(area_oi,2) <= 0.45
            count_p4 = count_p4 + 1;
            p4_area_ind(count_p4) = ff;
            p4_flux(count_p4) = fin_flux_kmpyr(p4_area_ind(count_p4));
            p4_area(count_p4) = fin_area_km2(p4_area_ind(count_p4));
            p4_phi(count_p4) = fin_phi(p4_area_ind(count_p4));
            p4_flux_var(count_p4) = fin_Sigmaflux(p4_area_ind(count_p4));
            p4_area_var(count_p4) = fin_SigmaArea(p4_area_ind(count_p4));
        elseif 0.45 < round(area_oi,2) && round(area_oi,2) <= 0.55
            count_p5 = count_p5 + 1;
            p5_area_ind(count_p5) = ff;
            p5_flux(count_p5) = fin_flux_kmpyr(p5_area_ind(count_p5));
            p5_area(count_p5) = fin_area_km2(p5_area_ind(count_p5));
            p5_phi(count_p5) = fin_phi(p5_area_ind(count_p5));
            p5_flux_var(count_p5) = fin_Sigmaflux(p5_area_ind(count_p5));
            p5_area_var(count_p5) = fin_SigmaArea(p5_area_ind(count_p5));
        elseif 0.55 < round(area_oi,2) && round(area_oi,2) <= 0.65
            count_p6 = count_p6 + 1;
            p6_area_ind(count_p6) = ff;
            p6_flux(count_p6) = fin_flux_kmpyr(p6_area_ind(count_p6));
            p6_area(count_p6) = fin_area_km2(p6_area_ind(count_p6));
            p6_phi(count_p6) = fin_phi(p6_area_ind(count_p6));
            p6_flux_var(count_p6) = fin_Sigmaflux(p6_area_ind(count_p6));
            p6_area_var(count_p6) = fin_SigmaArea(p6_area_ind(count_p6));
        elseif 0.65 < round(area_oi,2) && round(area_oi,2) <= 0.75
            count_p7 = count_p7 + 1;
            p7_area_ind(count_p7) = ff;
            p7_flux(count_p7) = fin_flux_kmpyr(p7_area_ind(count_p7));
            p7_area(count_p7) = fin_area_km2(p7_area_ind(count_p7));
            p7_phi(count_p7) = fin_phi(p7_area_ind(count_p7));
            p7_flux_var(count_p7) = fin_Sigmaflux(p7_area_ind(count_p7));
            p7_area_var(count_p7) = fin_SigmaArea(p7_area_ind(count_p7));
        elseif 0.75 < round(area_oi,2) && round(area_oi,2) <= 0.85
            count_p8 = count_p8 + 1;
            p8_area_ind(count_p8) = ff;
            p8_flux(count_p8) = fin_flux_kmpyr(p8_area_ind(count_p8));
            p8_area(count_p8) = fin_area_km2(p8_area_ind(count_p8));
            p8_phi(count_p8) = fin_phi(p8_area_ind(count_p8));
            p8_flux_var(count_p8) = fin_Sigmaflux(p8_area_ind(count_p8));
            p8_area_var(count_p8) = fin_SigmaArea(p8_area_ind(count_p8));
        elseif 0.85 < round(area_oi,2) && round(area_oi,2) <= 0.95
            count_p9 = count_p9 + 1;
            p9_area_ind(count_p9) = ff;
            p9_flux(count_p9) = fin_flux_kmpyr(p9_area_ind(count_p9));
            p9_area(count_p9) = fin_area_km2(p9_area_ind(count_p9));
            p9_phi(count_p9) = fin_phi(p9_area_ind(count_p9));
            p9_flux_var(count_p9) = fin_Sigmaflux(p9_area_ind(count_p9));
            p9_area_var(count_p9) = fin_SigmaArea(p9_area_ind(count_p9));
        elseif 1.05 < round(area_oi,2) && round(area_oi,2) <= 1.15
            count_p10 = count_p10 + 1;
            p10_area_ind(count_p10) = ff;
            p10_flux(count_p10) = fin_flux_kmpyr(p10_area_ind(count_p10));
            p10_area(count_p10) = fin_area_km2(p10_area_ind(count_p10));
            p10_phi(count_p10) = fin_phi(p10_area_ind(count_p10));
            p10_flux_var(count_p10) = fin_Sigmaflux(p10_area_ind(count_p10));
            p10_area_var(count_p10) = fin_SigmaArea(p10_area_ind(count_p10));
        elseif 1.15 < round(area_oi,2) && round(area_oi,2) <=  1.25
            count_p11 = count_p11 + 1;
            p11_area_ind(count_p11) = ff;
            p11_flux(count_p11) = fin_flux_kmpyr(p11_area_ind(count_p11));
            p11_area(count_p11) = fin_area_km2(p11_area_ind(count_p11));
            p11_phi(count_p11) = fin_phi(p11_area_ind(count_p11));
            p11_flux_var(count_p11) = fin_Sigmaflux(p11_area_ind(count_p11));
            p11_area_var(count_p11) = fin_SigmaArea(p11_area_ind(count_p11));
        elseif 1.25 < round(area_oi,2) && round(area_oi,2) <= 1.35
            count_p12 = count_p12 + 1;
            p12_area_ind(count_p12) = ff;
            p12_flux(count_p12) = fin_flux_kmpyr(p12_area_ind(count_p12));
            p12_area(count_p12) = fin_area_km2(p12_area_ind(count_p12));
            p12_phi(count_p12) = fin_phi(p12_area_ind(count_p12));
            p12_flux_var(count_p12) = fin_Sigmaflux(p12_area_ind(count_p12));
            p12_area_var(count_p12) = fin_SigmaArea(p12_area_ind(count_p12));
        elseif 1.35 < round(area_oi,2) && round(area_oi,2) <= 1.45
            count_p13 = count_p13 + 1;
            p13_area_ind(count_p13) = ff;
            p13_flux(count_p13) = fin_flux_kmpyr(p13_area_ind(count_p13));
            p13_area(count_p13) = fin_area_km2(p13_area_ind(count_p13));
            p13_phi(count_p13) = fin_phi(p13_area_ind(count_p13));
            p13_flux_var(count_p13) =fin_Sigmaflux(p13_area_ind(count_p13));
            p13_area_var(count_p13) =fin_SigmaArea(p13_area_ind(count_p13));
        end
    end
    vol_flux(1,jj,kk) = mean(p1_flux.*p1_area.*p1_phi);
    vol_flux(2,jj,kk) = mean(p2_flux.*p2_area.*p2_phi);
    vol_flux(3,jj,kk) = mean(p3_flux.*p3_area.*p3_phi);
    vol_flux(4,jj,kk) = mean(p4_flux.*p4_area.*p4_phi);
    vol_flux(5,jj,kk) = mean(p5_flux.*p5_area.*p5_phi);
    vol_flux(6,jj,kk) = mean(p6_flux.*p6_area.*p6_phi);
    vol_flux(7,jj,kk) = mean(p7_flux.*p7_area.*p7_phi);
    vol_flux(8,jj,kk) = mean(p8_flux.*p8_area.*p8_phi);
    vol_flux(9,jj,kk) = mean(p9_flux.*p9_area.*p9_phi);
    vol_flux(10,jj,kk) = mean(p10_flux.*p10_area.*p10_phi);
    vol_flux(11,jj,kk) = mean(p11_flux.*p11_area.*p11_phi);
    vol_flux(12,jj,kk) = mean(p12_flux.*p12_area.*p12_phi);
    vol_flux(13,jj,kk) = mean(p13_flux.*p13_area.*p13_phi);
    %error for flux
    error_volflux(1,jj,kk) = meanxy_sig_calc(xy_sig_calc(p1_flux_var,p1_flux,p1_area_var,p1_area,p1_flux.*p1_area.*p1_phi));
    error_volflux(2,jj,kk) = meanxy_sig_calc(xy_sig_calc(p2_flux_var,p2_flux,p2_area_var,p2_area,p2_flux.*p2_area.*p2_phi));
    error_volflux(3,jj,kk) = meanxy_sig_calc(xy_sig_calc(p3_flux_var,p3_flux,p3_area_var,p3_area,p3_flux.*p3_area.*p3_phi));
    error_volflux(4,jj,kk) = meanxy_sig_calc(xy_sig_calc(p4_flux_var,p4_flux,p4_area_var,p4_area,p4_flux.*p4_area.*p4_phi));
    error_volflux(5,jj,kk) = meanxy_sig_calc(xy_sig_calc(p5_flux_var,p5_flux,p5_area_var,p5_area,p5_flux.*p5_area.*p5_phi));
    error_volflux(6,jj,kk) = meanxy_sig_calc(xy_sig_calc(p6_flux_var,p6_flux,p6_area_var,p6_area,p6_flux.*p6_area.*p6_phi));
    error_volflux(7,jj,kk) = meanxy_sig_calc(xy_sig_calc(p7_flux_var,p7_flux,p7_area_var,p7_area,p7_flux.*p7_area.*p7_phi));
    error_volflux(8,jj,kk) = meanxy_sig_calc(xy_sig_calc(p8_flux_var,p8_flux,p8_area_var,p8_area,p8_flux.*p8_area.*p8_phi));
    error_volflux(9,jj,kk) = meanxy_sig_calc(xy_sig_calc(p9_flux_var,p9_flux,p9_area_var,p9_area,p9_flux.*p9_area.*p9_phi));
    error_volflux(10,jj,kk) = meanxy_sig_calc(xy_sig_calc(p10_flux_var,p10_flux,p10_area_var,p10_area,p10_flux.*p10_area.*p10_phi));
    error_volflux(11,jj,kk) = meanxy_sig_calc(xy_sig_calc(p11_flux_var,p11_flux,p11_area_var,p11_area,p11_flux.*p11_area.*p11_phi));
    error_volflux(12,jj,kk) = meanxy_sig_calc(xy_sig_calc(p12_flux_var,p12_flux,p12_area_var,p12_area,p12_flux.*p12_area.*p12_phi));
    error_volflux(13,jj,kk) = meanxy_sig_calc(xy_sig_calc(p13_flux_var,p13_flux,p13_area_var,p13_area,p13_flux.*p13_area.*p13_phi));
end

%% Two Permeabilities
rel_k16_flux = vol_flux(:,:,1);
rel_k15_flux = vol_flux(:,:,2);
rel_k16_errflux = error_volflux(:,:,1);
rel_k15_errflux = error_volflux(:,:,2);

A = zeros(13,5);
B = zeros(13,5);
A_ind = find(rel_k16_flux~=0);
B_ind = find(rel_k15_flux~=0);
A(A_ind) = 1;
B(B_ind) = 1;

C = A + B;
C2_ind = find(C==2);
C1_ind = find(C==1);
fin_mean_flux = rel_k16_flux + rel_k15_flux;
% Average flux
fin_mean_flux(C2_ind) = fin_mean_flux(C2_ind)./2;

% Average std of flux
s1 = rel_k16_errflux.^2+rel_k15_errflux.^2;
s2 = rel_k16_errflux + rel_k15_errflux;
fin_mean_sig = zeros([size(C,1) size(C,2)]);
fin_mean_sig(C2_ind) = sqrt(s1(C2_ind)/2);
fin_mean_sig(C1_ind) = s2(C1_ind);
obs_pockmark_n = [40; 114; 136; 89; 52; 23; 11; 1; 4; 2; 1; 1; 2];

cr_flux = fin_mean_flux.*obs_pockmark_n;
cr_fluxsig = fin_mean_sig.*obs_pockmark_n;

final_cr_sig = meanxy_sig_calc(sqrt(sum((cr_fluxsig.^2),1)))
final_cr_flux = mean(sum(cr_flux,1))

disp('DONE')
%% One permeability
obs_pockmark_n = [40; 114; 136; 89; 52; 23; 11; 1; 4; 2; 1; 1; 2];
rel_k15_flux = vol_flux(:,:,1);
rel_k15_errflux = error_volflux(:,:,1);
cr_flux = rel_k15_flux.*obs_pockmark_n;
cr_fluxsig = rel_k15_errflux.*obs_pockmark_n;
final_cr_flux = mean(sum(cr_flux,1))
final_cr_sig = meanxy_sig_calc(sqrt(sum((cr_fluxsig.^2),1)))
% Convert to mass
% liquid CO2
rho_l = 1.1e12;
% gas CO2
rho_g = 1.87e9;

kg2Pg = 1e12; % Convert kg to Pg
mass_rate_liq(1) = final_cr_flux*0.3*rho_l/kg2Pg;
mass_rate_liq(2) = final_cr_flux*0.7*rho_l/kg2Pg;
sig_rates_liq(1) = final_cr_sig*0.3*rho_l/kg2Pg;
sig_rates_liq(2) = final_cr_sig*0.7*rho_l/kg2Pg;

mass_rate_gas(1) = final_cr_flux*0.3*rho_g/kg2Pg;
mass_rate_gas(2) = final_cr_flux*0.7*rho_g/kg2Pg;
sig_rates_gas(1) = final_cr_sig*0.3*rho_g/kg2Pg;
sig_rates_gas(2) = final_cr_sig*0.7*rho_g/kg2Pg;
disp('Gas Phase Estimates: ')
format short
disp(['30%= ' num2str(round(mass_rate_gas(1),3,"significant")) ' +/- ' num2str(round(sig_rates_gas(1),2,"significant"))])
disp(['70%= ' num2str(round(mass_rate_gas(2),3,"significant")) ' +/- ' num2str(round(sig_rates_gas(2),2,"significant"))])

disp('Liquid Phase Estimates: ')
disp(['30%= ' num2str(round(mass_rate_liq(1),3,"significant")) ' +/- ' num2str(round(sig_rates_liq(1),2,"significant"))])
disp(['70%= ' num2str(round(mass_rate_liq(2),3,"significant")) ' +/- ' num2str(round(sig_rates_liq(2),2,"significant"))])

disp('DONE')