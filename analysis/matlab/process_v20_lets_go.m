% KG 10/22/2022

clear
close all

foldernames = [
"10_4La_4Li_hop1_00_to_04"
"12_4La_4Li_hop2_00_to_concerted"
"13_4La_4Li_hop4_step1"
"20_3La1Li_1La3Li_hop1_00_to_04_concerted_decided_by_me"
"22_3La1Li_1La3Li_hop2_00_to_concerted"
"33_3La1Li_1La3Li_extra_vertical_bottom_moving_3concerted"
"34_3La1Li_1La3Li_extra_vertical_center_moving_1single"
"38_3La1Li_1La3Li_extra_right_center_moving_3concerted"
"39_3La1Li_1La3Li_extra_right_left_moving_1single"
"41_3La1Li_1La3Li_extra_00_up_monitor_step1"
"42_3La1Li_1La3Li_extra_00_up_monitor_step2"
"43_3La1Li_1La3Li_extra_00_up_right_step1"
"44_3La1Li_1La3Li_extra_00_up_right_step2"
"45_3La1Li_1La3Li_extra_04_right_step1"
"46_3La1Li_1La3Li_extra_04_right_step2"
"47_3La1Li_1La3Li_extra_04_up_step1"
"48_3La1Li_1La3Li_extra_04_up_step2"
"60_2La2Li_2La2Li_1single1"
"62_2La2Li_2La2Li_concerted"
"64_2La2Li_2La2Li_04_up_step2"
"65_2La2Li_2La2Li_00_up_step1"
"67_2La2Li_2La2Li_00_up_concerted"
    ];

% Read the needed info from each folder
Nhops = length(foldernames);

for i=1:length(foldernames)
    %% Let's take care of oct rot first
    address = strcat(foldernames(i), "/ev.txt");
    ev = importdata(address);
    address = strcat(foldernames(i), "/phonon_contribs_discrete_LLTO.txt");
    temp = importdata(address);
    address = strcat(foldernames(i), "/hopping_atom.txt");
    hopping_atom  = importdata(address);
    address = strcat(foldernames(i), "/hopping_pathway.txt");
    hopping_pathway  = importdata(address);
    hopping_pathway_normalized = hopping_pathway/sqrt(sum(hopping_pathway.^2));
    myfreq_related_octrot = temp(:, 1);
    modenum = size(ev, 1);
    atomnum = modenum / 3;
    % total
    xproj = zeros(modenum, 1);
    yproj = zeros(modenum, 1);
    zproj = zeros(modenum, 1);

    for n=1:modenum
        for ii=1:atomnum
            xproj(n) = xproj(n) + ev((ii-1)*3+2, n)^2 + ev((ii-1)*3+3, n)^2;
            yproj(n) = yproj(n) + ev((ii-1)*3+1, n)^2 + ev((ii-1)*3+3, n)^2;
            zproj(n) = zproj(n) + ev((ii-1)*3+1, n)^2 + ev((ii-1)*3+2, n)^2;
        end
    end

    % PowTi
    PowTi = zeros(modenum, 1);

    for j=1:modenum
        sum1 = 0;
        for m = 4+4+1:4+4+8
            sum1 = sum1 + ev((m-1)*3+1, j)^2 + ev((m-1)*3+2, j)^2 + ev((m-1)*3+3, j)^2;
        end
        PowTi(j) = sum1;
    end

    % PowO
    PowO = zeros(modenum, 1);

    for j=1:modenum
        flagit = 0;
        sum1 = 0;
        for m = 4+4+8+1:4+4+8+24
            sum1 = sum1 + ev((m-1)*3+1, j)^2 + ev((m-1)*3+2, j)^2 + ev((m-1)*3+3, j)^2;
        end
        PowO(j) = sum1;
   
        ev(:,j) = ev(:,j) / sqrt(PowO(j));
    end

    % Oxygen

    Oxproj = zeros(modenum, 1);
    Oyproj = zeros(modenum, 1);
    Ozproj = zeros(modenum, 1);

    for n=1:modenum
        for ii=17:atomnum
            Oxproj(n) = Oxproj(n) + ev((ii-1)*3+2, n)^2 + ev((ii-1)*3+3, n)^2;
            Oyproj(n) = Oyproj(n) + ev((ii-1)*3+1, n)^2 + ev((ii-1)*3+3, n)^2;
            Ozproj(n) = Ozproj(n) + ev((ii-1)*3+1, n)^2 + ev((ii-1)*3+2, n)^2;
        end
    end

    rxy = zeros(modenum, 1);
    rxz = zeros(modenum, 1);
    ryz = zeros(modenum, 1);

    for n=1:modenum
        temp_rxy = Oxproj(n)/Oyproj(n);
        temp_rxz = Oxproj(n)/Ozproj(n);
        temp_ryz = Oyproj(n)/Ozproj(n);
        if temp_rxy > 1
            rxy(n) = temp_rxy;
        else
            rxy(n) = 1/temp_rxy;
        end
        if temp_rxz > 1
            rxz(n) = temp_rxz;
        else
            rxz(n) = 1/temp_rxz;
        end
        if temp_ryz > 1
            ryz(n) = temp_ryz;
        else
            ryz(n) = 1/temp_ryz;
        end 
    end

    Osum = zeros(modenum, 3);
    for n=1:modenum
        Osum(n, :) = sum([ev((17-1)*3+1:3:end, n), ev((17-1)*3+2:3:end, n), ev((17-1)*3+3:3:end, n)]);
    end

    vortTi_mag = zeros(8, 1);
    vortTi_mag_allmodes = zeros(modenum, 8);
    vortTi_ave = zeros(modenum, 1);

    for n=1:modenum
        vortTi = zeros(8, 3);

        % Ti1
        splitit = split(foldernames(i),'_');
        if ismember(str2num(splitit(1)), [60, 62, 64, 65, 67])
            ix1 = 24;
            ix2 = 20;
            iy1 = 32;
            iy2 = 30;
            iz1 = 40;
            iz2 = 39;
        else
            ix1 = 23;
            ix2 = 26;
            iy1 = 21;
            iy2 = 27;
            iz1 = 34;
            iz2 = 36;
        end

        vortTi(1,1) = (ev((iy2-1)*3+3, n) - ev((iy1-1)*3+3, n)) - (ev((iz2-1)*3+2, n) - ev((iz1-1)*3+2, n));
        vortTi(1,2) = (ev((iz2-1)*3+1, n) - ev((iz1-1)*3+1, n)) - (ev((ix2-1)*3+3, n) - ev((ix1-3)*3+3, n));
        vortTi(1,3) = (ev((ix2-1)*3+2, n) - ev((ix1-1)*3+2, n)) - (ev((iy2-1)*3+1, n) - ev((iy1-1)*3+1, n));
    
        % Ti2
        if ismember(str2num(splitit(1)), [60, 62, 64, 65, 67])
            ix1 = 22;
            ix2 = 18;
            iy1 = 30;
            iy2 = 32;
            iz1 = 38;
            iz2 = 37;
        else
            ix1 = 19;
            ix2 = 17;
            iy1 = 27;
            iy2 = 21;
            iz1 = 18;
            iz2 = 38;
        end

        vortTi(2,1) = (ev((iy2-1)*3+3, n) - ev((iy1-1)*3+3, n)) - (ev((iz2-1)*3+2, n) - ev((iz1-1)*3+2, n));
        vortTi(2,2) = (ev((iz2-1)*3+1, n) - ev((iz1-1)*3+1, n)) - (ev((ix2-1)*3+3, n) - ev((ix1-3)*3+3, n));
        vortTi(2,3) = (ev((ix2-1)*3+2, n) - ev((ix1-1)*3+2, n)) - (ev((iy2-1)*3+1, n) - ev((iy1-1)*3+1, n));

        % Ti3
        if ismember(str2num(splitit(1)), [60, 62, 64, 65, 67])
            ix1 = 18;
            ix2 = 22;
            iy1 = 26;
            iy2 = 28;
            iz1 = 34;
            iz2 = 33;
        else
            ix1 = 17;
            ix2 = 19;
            iy1 = 31;
            iy2 = 24;
            iz1 = 29;
            iz2 = 33;
        end

        vortTi(3,1) = (ev((iy2-1)*3+3, n) - ev((iy1-1)*3+3, n)) - (ev((iz2-1)*3+2, n) - ev((iz1-1)*3+2, n));
        vortTi(3,2) = (ev((iz2-1)*3+1, n) - ev((iz1-1)*3+1, n)) - (ev((ix2-1)*3+3, n) - ev((ix1-3)*3+3, n));
        vortTi(3,3) = (ev((ix2-1)*3+2, n) - ev((ix1-1)*3+2, n)) - (ev((iy2-1)*3+1, n) - ev((iy1-1)*3+1, n));
    
        % Ti4
        if ismember(str2num(splitit(1)), [60, 62, 64, 65, 67])
            ix1 = 20;
            ix2 = 24;
            iy1 = 28;
            iy2 = 26;
            iz1 = 36;
            iz2 = 35;
        else
            ix1 = 26;
            ix2 = 23;
            iy1 = 24;
            iy2 = 31;
            iz1 = 25;
            iz2 = 40;
        end

        vortTi(4,1) = (ev((iy2-1)*3+3, n) - ev((iy1-1)*3+3, n)) - (ev((iz2-1)*3+2, n) - ev((iz1-1)*3+2, n));
        vortTi(4,2) = (ev((iz2-1)*3+1, n) - ev((iz1-1)*3+1, n)) - (ev((ix2-1)*3+3, n) - ev((ix1-3)*3+3, n));
        vortTi(4,3) = (ev((ix2-1)*3+2, n) - ev((ix1-1)*3+2, n)) - (ev((iy2-1)*3+1, n) - ev((iy1-1)*3+1, n));
    
        % Ti5
        if ismember(str2num(splitit(1)), [60, 62, 64, 65, 67])
            ix1 = 23;
            ix2 = 19;
            iy1 = 31;
            iy2 = 29;
            iz1 = 39;
            iz2 = 40;
        else
            ix1 = 37;
            ix2 = 32;
            iy1 = 35;
            iy2 = 39;
            iz1 = 36;
            iz2 = 34;
        end

        vortTi(5,1) = (ev((iy2-1)*3+3, n) - ev((iy1-1)*3+3, n)) - (ev((iz2-1)*3+2, n) - ev((iz1-1)*3+2, n));
        vortTi(5,2) = (ev((iz2-1)*3+1, n) - ev((iz1-1)*3+1, n)) - (ev((ix2-1)*3+3, n) - ev((ix1-3)*3+3, n));
        vortTi(5,3) = (ev((ix2-1)*3+2, n) - ev((ix1-1)*3+2, n)) - (ev((iy2-1)*3+1, n) - ev((iy1-1)*3+1, n));
    
        % Ti6
        if ismember(str2num(splitit(1)), [60, 62, 64, 65, 67])
            ix1 = 21;
            ix2 = 17;
            iy1 = 29;
            iy2 = 19;
            iz1 = 37;
            iz2 = 38;
        else
            ix1 = 30;
            ix2 = 22;
            iy1 = 39;
            iy2 = 35;
            iz1 = 38;
            iz2 = 18;
        end
        
        vortTi(6,1) = (ev((iy2-1)*3+3, n) - ev((iy1-1)*3+3, n)) - (ev((iz2-1)*3+2, n) - ev((iz1-1)*3+2, n));
        vortTi(6,2) = (ev((iz2-1)*3+1, n) - ev((iz1-1)*3+1, n)) - (ev((ix2-1)*3+3, n) - ev((ix1-3)*3+3, n));
        vortTi(6,3) = (ev((ix2-1)*3+2, n) - ev((ix1-1)*3+2, n)) - (ev((iy2-1)*3+1, n) - ev((iy1-1)*3+1, n));

        % Ti7
        if ismember(str2num(splitit(1)), [60, 62, 64, 65, 67])
            ix1 = 17;
            ix2 = 21;
            iy1 = 25;
            iy2 = 27;
            iz1 = 33;
            iz2 = 34;
        else
            ix1 = 22;
            ix2 = 30;
            iy1 = 20;
            iy2 = 28;
            iz1 = 33;
            iz2 = 29;
        end

        vortTi(7,1) = (ev((iy2-1)*3+3, n) - ev((iy1-1)*3+3, n)) - (ev((iz2-1)*3+2, n) - ev((iz1-1)*3+2, n));
        vortTi(7,2) = (ev((iz2-1)*3+1, n) - ev((iz1-1)*3+1, n)) - (ev((ix2-1)*3+3, n) - ev((ix1-3)*3+3, n));
        vortTi(7,3) = (ev((ix2-1)*3+2, n) - ev((ix1-1)*3+2, n)) - (ev((iy2-1)*3+1, n) - ev((iy1-1)*3+1, n));

        % Ti8
        if ismember(str2num(splitit(1)), [60, 62, 64, 65, 67])
            ix1 = 19;
            ix2 = 23;
            iy1 = 27;
            iy2 = 25;
            iz1 = 35;
            iz2 = 36;
        else
            ix1 = 32;
            ix2 = 37;
            iy1 = 28;
            iy2 = 20;
            iz1 = 40;
            iz2 = 25;
        end

        vortTi(8,1) = (ev((iy2-1)*3+3, n) - ev((iy1-1)*3+3, n)) - (ev((iz2-1)*3+2, n) - ev((iz1-1)*3+2, n));
        vortTi(8,2) = (ev((iz2-1)*3+1, n) - ev((iz1-1)*3+1, n)) - (ev((ix2-1)*3+3, n) - ev((ix1-3)*3+3, n));
        vortTi(8,3) = (ev((ix2-1)*3+2, n) - ev((ix1-1)*3+2, n)) - (ev((iy2-1)*3+1, n) - ev((iy1-1)*3+1, n));
    
        for ii=1:8
            vortTi_mag(ii) = (vortTi(ii,1)^2 +  vortTi(ii,2)^2 +  vortTi(ii,3)^2)^.5;
        end
        vortTi_ave(n) = mean(vortTi_mag);

        vortTi_mag_allmodes(n, :) = vortTi_mag;
    end

    tagOctRot(1:modenum, i) = 0;
    tag1 = zeros(modenum, 1);
    tag2 = zeros(modenum, 1);
    tag3 = zeros(modenum, 1);
    tag4 = zeros(modenum, 1);
    tag5 = zeros(modenum, 1);
    tag6 = zeros(modenum, 1);
    tag7 = zeros(modenum, 1);

    for n=1:modenum
        %r_criteria = 1.1;
        %if (rxy(n)>r_criteria*rxz(n) & rxy(n)>r_criteria*ryz(n)) | (rxz(n)>r_criteria*rxy(n) & rxz(n)>r_criteria*ryz(n)) | (ryz(n)>r_criteria*rxz(n) & ryz(n)>r_criteria*rxy(n))
        %    tag1(n) = 1;
        %end

        r_criteria = 1.05;
        if (r_criteria*rxy(n)<rxz(n) & r_criteria*rxy(n)<ryz(n)) | (r_criteria*rxz(n)<rxy(n) & r_criteria*rxz(n)<ryz(n)) | (r_criteria*ryz(n)<rxz(n) & r_criteria*ryz(n)<rxy(n))
            tag1(n) = 1;
        end

        Oproj_criteria = 1.3;
        if (Oxproj(n)>Oproj_criteria*Oyproj(n) & Oxproj(n)>Oproj_criteria*Ozproj(n)) | (Oyproj(n)>Oproj_criteria*Oxproj(n) & Oyproj(n)>Oproj_criteria*Ozproj(n)) | (Ozproj(n)>Oproj_criteria*Oyproj(n) & Ozproj(n)>Oproj_criteria*Oxproj(n))
            tag1(n) = 1;
        end

        Osum_criteria = 1e-3;
        if abs(Osum(n, 1))<Osum_criteria & abs(Osum(n, 2))<Osum_criteria & abs(Osum(n, 3))<Osum_criteria
            tag2(n) = 1;
        end

        Oproj_mag_criteria = .15;
        if Oxproj(n)>Oproj_mag_criteria & Oyproj(n)>Oproj_mag_criteria & Ozproj(n)>Oproj_mag_criteria
            tag3(n) = 1;
        end

        freq_criteria = 10;
        if myfreq_related_octrot(n) < freq_criteria
            tag4(n) = 1;
        end

        PowTi_criteria = 0.2;
        if PowTi(n) < PowTi_criteria
            tag5(n) = 1;
        end

        vortivity_criteria = 0.25;
        if vortTi_ave(n)^3  > vortivity_criteria
            tag6(n) = 1;
        end

        PowO_criteria = 0.2;
        if PowO(n) > PowO_criteria
            tag7(n) = 1;
        end

        %if tag4(n)==1 & tag6(n)==1 & tag7(n)==1
        %if tag6(n)==1 & tag7(n)==1
        %if tag1(n)==1 & tag6(n)==1 & tag7(n)==1
        if tag6(n)==1 & tag7(n)==1
            tagOctRot(n, i) = 1;
        end
    end
    
    %% MB
    address = strcat(foldernames(i), "/MB.txt");
    temp = importdata(address);
    MB(i,1) = temp(2);
    MB_symbol(i,1) = temp(1);
    if MB(i,1) < 0.05
        MB(i,1) = .2;
    end
    
    %% contrib to hop - each mode
    address = strcat(foldernames(i), "/phonon_contribs_discrete_LLTO.txt");
    temp = importdata(address);
    freq(:, i) = temp(:, 1); % frequency
    contrib_to_indi_hop(:, i) = temp(:, 2); % contribution to individual hop
    
    %address = strcat(foldernames(i), "/bn_modal_percent.txt");
    address = strcat(foldernames(i), "/bn_modal_percent_wo_abs.txt");
    temp = importdata(address);
    bn_percent(:, i) = temp; % bn_modal
    
    address = strcat(foldernames(i), "/Secondary_forces.txt");
    forces = importdata(address);
    FmagF(:, i) = forces(:, 1);
    FmagF_proj(:, i) = forces(:, 2);
    FmagB(:, i) = forces(:, 3);
    FmagB_proj(:, i) = forces(:, 4);
    
    sum1 = 0;
    accum_contribs_total(1:length(freq(:,i)), i) = 0;
    for j=1:length(freq(:,i))
        accum_contribs_total(j, i) = sum1;
        sum1 = sum1 + contrib_to_indi_hop(j, i);
    end
    %for j=1:length(freq)
    %    if freq(i, j) < -0.2
    %        freq(i, j) = -0.1;
    %    end
    %end
    
    %% contrib to hop - each mode - sorted
    address = strcat(foldernames(i), "/phonon_contribs_discrete_LLTO.txt");
    temp = importdata(address);
    contrib_to_indi_hop_sorted(:, i) = sort(temp(:, 2), 'descend'); % contribution to individual hop
    sum1 = 0;
    accum_contribs_total_sorted(1:length(freq(:,i)), i) = 0;
    for j=1:length(freq(:,i))
        accum_contribs_total_sorted(j, i) = sum1;
        sum1 = sum1 + contrib_to_indi_hop_sorted(j, i);
    end
    %for j=1:length(freq)
    %    if freq(i, j) < -0.2
    %        freq(i, j) = -0.1;
    %    end
    %end
    
    %% octahedral (O) energies
    address = strcat(foldernames(i), "/ev.txt");
    ev = importdata(address);
    NLa = 4;
    NLi = 4;
    NTi = 8;
    NO = 24;
    Natoms = 4 + 4 + 8 + 24;
    for j=1:length(freq(:,i))
        sum1 = 0;
        for m = Natoms-NO+1:Natoms
            sum1 = sum1 + ev((m-1)*3+1, j)^2 + ev((m-1)*3+2, j)^2 + ev((m-1)*3+3, j)^2;
        end
        powO(j, i) = sum1;
    end
   
    %% Ti energies
    address = strcat(foldernames(i), "/ev.txt");
    ev = importdata(address);
    NLa = 4;
    NLi = 4;
    NTi = 8;
    NO = 24;
    Natoms = 4 + 4 + 8 + 24;
    for j=1:length(freq(:,i))
        sum1 = 0;
        for m = NLa+NLi+1:NLa+NLi+NTi
            sum1 = sum1 + ev((m-1)*3+1, j)^2 + ev((m-1)*3+2, j)^2 + ev((m-1)*3+3, j)^2;
        end
        powTi(j, i) = sum1;
    end
    
    %% La energies
    address = strcat(foldernames(i), "/ev.txt");
    ev = importdata(address);
    NLa = 4;
    NLi = 4;
    NTi = 8;
    NO = 24;
    Natoms = 4 + 4 + 8 + 24;
    for j=1:length(freq(:,i))
        sum1 = 0;
        for m = 1:NLa
            sum1 = sum1 + ev((m-1)*3+1, j)^2 + ev((m-1)*3+2, j)^2 + ev((m-1)*3+3, j)^2;
        end
        powLa(j, i) = sum1;
    end
    %% Li energies
    address = strcat(foldernames(i), "/ev.txt");
    ev = importdata(address);
    NLa = 4;
    NLi = 4;
    NTi = 8;
    NO = 24;
    Natoms = 4 + 4 + 8 + 24;
    for j=1:length(freq(:,i))
        sum1 = 0;
        for m = NLa+1:NLa+NLi
            sum1 = sum1 + ev((m-1)*3+1, j)^2 + ev((m-1)*3+2, j)^2 + ev((m-1)*3+3, j)^2;
            if m == hopping_atom
                powLi_hopping(j, i) = ev((m-1)*3+1, j)^2 + ev((m-1)*3+2, j)^2 + ev((m-1)*3+3, j)^2;
                evhopping = [ev((m-1)*3+1, j), ev((m-1)*3+2, j), ev((m-1)*3+3, j)];
                powLi_hopping_projected(j, i) = abs(dot(evhopping, hopping_pathway_normalized));
            end
        end
        powLi(j, i) = sum1;
    end
    
    %% contrib to hop - each mode & each_atom
    oct_rot_max_freq = 8;
    oct_rot_max_contrib = 0.1;
    Natoms = 4 + 4 + 8 + 24;
    contribs_all_atoms_with_details(1:length(freq(:,i)), i, Natoms) = 0;
    contribs_all_atoms(1:length(freq(:,i)), i) = 0;
    contribs_all_La(1:length(freq(:,i)), i) = 0;
    contribs_all_Li(1:length(freq(:,i)), i) = 0;
    contribs_all_Ti(1:length(freq(:,i)), i) = 0;
    contribs_all_O(1:length(freq(:,i)), i) = 0;
    contribs_all_oct_rot(1:length(freq(:,i)), i) = 0;
    contribs_all_oct_rot_v2(1:length(freq(:,i)), i) = 0;
    
    for j=1:length(freq(:,i))
        for m = 1:Natoms
            factor = ev((m-1)*3+1, j)^2 + ev((m-1)*3+2, j)^2 + ev((m-1)*3+3, j)^2;
            contribs_all_atoms_with_details(j, i, m) = factor * contrib_to_indi_hop(j, i);
        end
        contribs_all_atoms(j, i) = sum(contribs_all_atoms_with_details(j, i, :));
        contribs_all_La(j, i) = sum(contribs_all_atoms_with_details(j, i, 1:4));
        contribs_all_Li(j, i) = sum(contribs_all_atoms_with_details(j, i, 5:8));
        contribs_all_Ti(j, i) = sum(contribs_all_atoms_with_details(j, i, 9:16));
        contribs_all_O(j, i) = sum(contribs_all_atoms_with_details(j, i, 17:40));
        if (freq(j, i) < oct_rot_max_freq & powTi(j, i) < oct_rot_max_contrib)
            contribs_all_oct_rot(j, i) = contrib_to_indi_hop(j, i);
        end
        if tagOctRot(j, i) == 1
            contribs_all_oct_rot_v2(j, i) = contrib_to_indi_hop(j, i);
        end
    end
    sorted_contribs(1:length(freq(:,i)), i) = 0;
    accum_contribs_all_atoms(1:length(freq(:,i)), i) = 0;
    [dummy, sorted_contribs(:, i)] = sort(contribs_all_atoms(:, i), 'descend');
    accum_contribs_all_atoms_10p(1:length(freq(:,i)), i) = 0;
    accum_contribs_all_atoms_BE_corrected(1:length(freq(:,i)), i) = 0;
    accum_contribs_all_atoms_BE_corrected_10p(1:length(freq(:,i)), i) = 0;
    accum_contribs_all_atoms_BE_corrected_oct_rot(1:length(freq(:,i)), i) = 0;
    accum_contribs_all_La(1:length(freq(:,i)), i) = 0;
    accum_contribs_all_Li(1:length(freq(:,i)), i) = 0;
    accum_contribs_all_Ti(1:length(freq(:,i)), i) = 0;
    accum_contribs_all_O(1:length(freq(:,i)), i) = 0;
    accum_contribs_all_oct_rot(1:length(freq(:,i)), i) = 0;
    accum_contribs_all_oct_rot_v2(1:length(freq(:,i)), i) = 0;
    sumtot = 0;
    sumtot_10p = 0;
    sumtot_BE_corrected = 0;
    sumtot_BE_corrected_10p = 0;
    sumtot_BE_corrected_oct_rot = 0;
    sumLa = 0;
    sumLi = 0;
    sumTi = 0;
    sumO = 0;
    sumoctrot = 0;
    sumoctrot_v2 = 0;
    for j=1:length(freq(:, i))
        accum_contribs_all_atoms(j, i) = sumtot; % This is technically equal to accum_contribs_total
        accum_contribs_all_atoms_10p(j, i) = sumtot_10p;
        accum_contribs_all_atoms_BE_corrected(j, i) = sumtot_BE_corrected;
        accum_contribs_all_atoms_BE_corrected_10p(j, i) = sumtot_BE_corrected_10p;
        accum_contribs_all_atoms_BE_corrected_oct_rot(j, i) = sumtot_BE_corrected_oct_rot;
        accum_contribs_all_La(j, i) = sumLa;
        accum_contribs_all_Li(j, i) = sumLi;
        accum_contribs_all_Ti(j, i) = sumTi;
        accum_contribs_all_O(j, i) = sumO;
        accum_contribs_all_oct_rot(j, i) = sumoctrot;
        accum_contribs_all_oct_rot_v2(j, i) = sumoctrot_v2;
        
        h = 6.62607015e-34; % J/Hz
        kB = 1.380649e-23; % J/K
        Temperature = 300; % K
        BE_x = (h * abs(freq(j, i)*10^12))/(kB * Temperature);
        BE_factor_matrix(j,i) = BE_x^2*exp(BE_x)/(exp(BE_x)-1)^2;
        BE_factor = BE_factor_matrix(j,i);
        
        sumtot = sumtot + contribs_all_atoms(j,i);
        if ismember(j, sorted_contribs(1:6 ,i))
            sumtot_10p = sumtot_10p + contribs_all_atoms(j,i);
            sumtot_BE_corrected_10p = sumtot_BE_corrected_10p + contribs_all_atoms(j,i) * BE_factor;
        end
        sumtot_BE_corrected = sumtot_BE_corrected + contribs_all_atoms(j,i) * BE_factor;
        sumLa = sumLa + contribs_all_La(j,i);
        sumLi = sumLi + contribs_all_Li(j,i);
        sumTi = sumTi + contribs_all_Ti(j,i);
        sumO = sumO + contribs_all_O(j,i);
        sumoctrot = sumoctrot + contribs_all_oct_rot(j,i);
        sumoctrot_v2 = sumoctrot_v2 + contribs_all_oct_rot_v2(j,i);
        sumtot_BE_corrected_oct_rot = sumtot_BE_corrected_oct_rot + contribs_all_oct_rot_v2(j,i) * BE_factor;
    end
    
    %% plot DOS (total & partial)

    Npts_DOS = 40;

    counter_DOS_total(1:Npts_DOS, i) = 0;
    counter_DOS_La(1:Npts_DOS, i) = 0;
    counter_DOS_Li(1:Npts_DOS, i) = 0;
    counter_DOS_Ti(1:Npts_DOS, i) = 0;
    counter_DOS_O(1:Npts_DOS, i) = 0;
    counter_DOS_octrot(1:Npts_DOS, i) = 0;
    Pow_temp_atoms = zeros(Natoms, 1);
    w_counter(1:Npts_DOS, i) = 0;
    max_omega_THz = max(freq(:));
    domega = max_omega_THz/Npts_DOS;

    omega_THz = abs(freq(:, i));

    w_counter(:, i) = (1:Npts_DOS)'*domega-domega*.5;
    for j = 1:length(freq(:, i))
        if (floor(omega_THz(j)/domega)<Npts_DOS)
            indx = floor(omega_THz(j)/domega) + 1;
        else
            indx = floor(omega_THz(j)/domega);
        end
        counter_DOS_total(indx, i) = counter_DOS_total(indx, i) + 1;
        if tagOctRot(j, i) == 1
            counter_DOS_octrot(indx, i) = counter_DOS_octrot(indx, i) + 1;
        end
        for m = 1:Natoms
            Pow_temp_atoms(m) = ev((m-1)*3+1, j)^2 + ev((m-1)*3+2, j)^2 + ev((m-1)*3+3, j)^2;
        end
        counter_DOS_La(indx, i) = counter_DOS_La(indx, i) + sum(Pow_temp_atoms(1:4));
        counter_DOS_Li(indx, i) = counter_DOS_Li(indx, i) + sum(Pow_temp_atoms(5:8));
        counter_DOS_Ti(indx, i) = counter_DOS_Ti(indx, i) + sum(Pow_temp_atoms(9:16));
        counter_DOS_O(indx, i) = counter_DOS_O(indx, i) + sum(Pow_temp_atoms(17:40));
    end
end

%% MB plot with no weighting (just MB values)
figure
set(gcf, 'color', 'w', 'name', 'MB plot with no weighting (just MB values)')
hold on
[sorted_MB, index_sorted_MB] = sort(MB);
MB_symbol_sorted = MB_symbol(index_sorted_MB);
markersize = 9;
for i=1:length(sorted_MB)
    if MB_symbol_sorted(i) == 1
        plot(i, sorted_MB(i), 'bo', 'markersize', markersize, 'markerfacecolor', 'b', 'markeredgecolor', 'b')
    elseif MB_symbol_sorted(i) == 2
        plot(i, sorted_MB(i), 'rs', 'markersize', markersize, 'markerfacecolor', 'r', 'markeredgecolor', 'r')
    else
        plot(i, sorted_MB(i), 'kdiamond', 'markersize', markersize, 'markerfacecolor', 'k', 'markeredgecolor', 'k')
    end
end
axis square
set(gca, 'fontsize', 17)
xlabel('Hop #')
ylabel('Migration barrier (eV)')
box on

%% MB plot using histogram and no weight
dMB = .2;
Npts = floor(max(sorted_MB)/dMB) + 1;
counter(Npts)= 0;

MB_bar_centers = (1:Npts)*dMB-dMB*.5;
for i = 1:length(sorted_MB)
  if (floor(sorted_MB(i)/dMB)<Npts)
    indx = floor(sorted_MB(i)/dMB) + 1;
  else
    indx = floor(sorted_MB(i)/dMB);
  end
  %indx
  counter(indx) = counter(indx) + 1;
end
%figure;
figure
set(gcf, 'color', 'w', 'name', 'MB plot using histogram and no weight')
hold on
bar(MB_bar_centers, counter)
xticks([MB_bar_centers])
axis square
set(gca, 'fontsize', 17)
xlabel('Migration barrier (eV)')
ylabel('Counts')
box on

figure
set(gcf, 'color', 'w', 'name', 'MB plot using histogram and no weight (normalized)')
hold on
bar(MB_bar_centers, counter/sum(counter)*100)
xticks([MB_bar_centers])
axis square
set(gca, 'fontsize', 17)
xlabel('Migration barrier (eV)')
ylabel('Normalized counts')
box on

%% MB plot using histogram and with weight

kB = 8.617e-5; % eV/K
T = 900; % K

dMB = .2;
Npts = floor(max(sorted_MB)/dMB) + 1;
counter2(Npts)= 0;

MB_bar_centers = (1:Npts)*dMB-dMB*.5;
for i = 1:length(sorted_MB)
  if (floor(sorted_MB(i)/dMB)<Npts)
    indx = floor(sorted_MB(i)/dMB) + 1;
  else
    indx = floor(sorted_MB(i)/dMB);
  end
  %indx
  counter2(indx) = counter2(indx) + 1*exp(-sorted_MB(i)/(kB*T));
end

figure
set(gcf, 'color', 'w', 'name', 'MB plot using histogram and with weight')
hold on
bar(MB_bar_centers, counter2/sum(counter2)*100)
xticks([MB_bar_centers])
axis square
set(gca, 'fontsize', 17)
xlabel('Migration barrier (eV)')
ylabel('Normalized counts wighted by Boltzmann factor')
box on

%% Contribution plots colored by MB
figure
set(gcf, 'color', 'w', 'name', 'Contribution plots colored by MB')
hold on
freq_ave = mean(freq, 2);
freq_ave = importdata('freq_ave2.txt');

for i=1:length(foldernames)
  x = freq_ave';
  y = accum_contribs_total(:,i)';
  xx=[x;x];
  yy=[y;y];
  zz=ones(size(xx))*MB(i);
  surf(xx,yy,zz,zz,'EdgeColor','interp') %// color binded to "z" values
end

colormap('jet')
view(2) %// view(0,90)
%colorbar
axis square
colorbar
set(gca, 'fontsize', 17)
xlabel('Frequency (THz)')
ylabel({'Accumulation of normalized';'contribution to ion hop'})
box on

%% Contribution plots colored by MB (average & errorbar)
figure
set(gcf, 'color', 'w', 'name', 'Contribution plots colored by MB (average & errorbar)')
hold on

for i=1:length(foldernames)
  x = freq_ave';
  y = accum_contribs_total(:,i)';
  xx=[x;x];
  yy=[y;y];
  zz=ones(size(xx))*MB(i);
  surf(xx,yy,zz,zz,'EdgeColor','interp') %// color binded to "z" values
end

colormap('jet')
view(2) %// view(0,90)
%colorbar
axis square
colorbar
set(gca, 'fontsize', 17)
xlabel('Frequency (THz)')
ylabel({'Accumulation of normalized';'contribution to ion hop'})
box on

contrib_ave = mean(accum_contribs_total');
contrib_std = std(accum_contribs_total');

x = [freq_ave', fliplr(freq_ave')];
inBetween = [contrib_ave-contrib_std, fliplr(contrib_ave+contrib_std)];

%fill(x, inBetween, [.5 .5 .5], 'FaceAlpha', .5);
plot(freq_ave', contrib_ave, 'k-', 'linewidth', 3)

%% Contribution plots colored by MB (only MB < 0.6)
figure
set(gcf, 'color', 'w', 'name', 'Contribution plots colored by MB (only MB < 0.6)')
hold on

for i=1:length(foldernames)
  x = freq_ave';
  y = accum_contribs_total(:,i)';
  xx=[x;x];
  yy=[y;y];
  if MB(i) < 0.6
      zz=ones(size(xx))*MB(i);
      surf(xx,yy,zz,zz,'EdgeColor','interp') %// color binded to "z" values
  end
end

colormap('jet')
view(2) %// view(0,90)
%colorbar
axis square
colorbar
set(gca, 'fontsize', 17)
xlabel('Frequency (THz)')
ylabel({'Accumulation of normalized';'contribution to ion hop'})
box on

%% Contribution plots colored by MB (average & errorbar) - sorted
figure
set(gcf, 'color', 'w', 'name', 'Contribution plots colored by MB (average & errorbar) - sorted')
hold on

for i=1:length(foldernames)
  modenums_percent = [1:length(freq(:,1))]/length(freq(:,1))*100;
  x = modenums_percent;
  y = accum_contribs_total_sorted(:,i)';
  xx=[x;x];
  yy=[y;y];
  zz=ones(size(xx))*MB(i);
  surf(xx,yy,zz,zz,'EdgeColor','interp') %// color binded to "z" values
end

colormap('jet')
view(2) %// view(0,90)
%colorbar
axis square
colorbar
set(gca, 'fontsize', 17)
xlabel('Mode # (% of total modes)')
ylabel({'Accumulation of normalized';'contribution to ion hop'})
box on

contrib_ave = mean(accum_contribs_total_sorted');
contrib_std = std(accum_contribs_total_sorted');

x = [modenums_percent, fliplr(modenums_percent)];
inBetween = [contrib_ave-contrib_std, fliplr(contrib_ave+contrib_std)];

%fill(x, inBetween, [.5 .5 .5], 'FaceAlpha', .5);
plot(modenums_percent, contrib_ave, 'k-', 'linewidth', 3)

%% Contribution plots colored by index
% figure
% set(gcf, 'color', 'w')
% hold on
% 
% for i=1:length(foldernames)
%   x = freq(:,1)';
%   y = accum_contribs(:,i)';
%   xx=[x;x];
%   yy=[y;y];
%   zz=ones(size(xx)) * i;
%   surf(xx,yy,zz,zz,'EdgeColor','interp') %// color binded to "z" values
% end
% 
% colormap('jet')
% view(2) %// view(0,90)
% %colorbar
% axis square
% colorbar
% set(gca, 'fontsize', 17)
% xlabel('Frequency (THz)')
% ylabel('Accumulation of normalized contribution to ion hop')
% box on

%% Octahedron deformation energy (x-axis O energy)
figure
set(gcf, 'color', 'w', 'name', 'O energy (x-axis O energy)')
hold on
markersize = 9;

for i=1:length(foldernames)
    fmin = min(freq(:,i));
    fmax = max(freq(:,i));
    for j=1:length(freq(:,i))
        xRGB = (freq(j,i) - fmin)/(fmax - fmin);
        plot(powO(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', [1-xRGB 0 0+xRGB], 'markersize', markersize)
    end
    %plot(powO(:, i), contrib_to_indi_hop(:, i), 'o', 'markerfacecolor', [1-x 0 0+x], 'markeredgecolor', [1-x 0 0+x], 'markersize', markersize)
end

axis square
set(gca, 'fontsize', 17)
xlabel('Normalized octahedron deformation energy for each phonon')
ylabel({'Normalized contribution';'to ion hop'})
box on

%% Octahedron deformation energy (x-axis frequency)
figure
set(gcf, 'color', 'w', 'name', 'O energy (x-axis frequency)')
hold on
markersize = 9;

for i=1:length(foldernames)
    Omin = min(powO(:,i));
    Omax = max(powO(:,i));
    for j=1:length(freq(:,i))
        xRGB = (powO(j,i) - Omin)/(Omax - Omin);
        plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', [1-xRGB 0 0+xRGB], 'markersize', markersize)
    end
    %plot(powO(:, i), contrib_to_indi_hop(:, i), 'o', 'markerfacecolor', [1-x 0 0+x], 'markeredgecolor', [1-x 0 0+x], 'markersize', markersize)
end

axis square
set(gca, 'fontsize', 17)
xlabel('Frequency (THz)')
ylabel({'Normalized contribution';'to ion hop'})
box on
xticks([0:5:30])

%% Ti energy (x-axis Ti energy)
figure
set(gcf, 'color', 'w', 'name', 'Ti energy (x-axis Ti energy)')
hold on
markersize = 9;

for i=1:length(foldernames)
    fmin = min(freq(:,i));
    fmax = max(freq(:,i));
    for j=1:length(freq(:,i))
        xRGB = (freq(j,i) - fmin)/(fmax - fmin);
        plot(powTi(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', [1-xRGB 0 0+xRGB], 'markersize', markersize)
    end
    %plot(powO(:, i), contrib_to_indi_hop(:, i), 'o', 'markerfacecolor', [1-x 0 0+x], 'markeredgecolor', [1-x 0 0+x], 'markersize', markersize)
end

axis square
set(gca, 'fontsize', 17)
xlabel('Normalized Ti energy for each phonon')
ylabel({'Normalized contribution';'to ion hop'})
box on

%% Ti energy (x-axis frequency)
figure
set(gcf, 'color', 'w', 'name', 'Ti energy (x-axis frequency)')
hold on
markersize = 9;

for i=1:length(foldernames)
    Timin = min(powTi(:,i));
    Timax = max(powTi(:,i));
    for j=1:length(freq(:,i))
        xRGB = (powTi(j,i) - Timin)/(Timax - Timin);
        plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', [1-xRGB 0 0+xRGB], 'markersize', markersize)
    end
    %plot(powO(:, i), contrib_to_indi_hop(:, i), 'o', 'markerfacecolor', [1-x 0 0+x], 'markeredgecolor', [1-x 0 0+x], 'markersize', markersize)
end

axis square
set(gca, 'fontsize', 17)
xlabel('Frequency (THz)')
ylabel({'Normalized contribution';'to ion hop'})
box on
xticks([0:5:30])

%% La energy (x-axis La energy)
figure
set(gcf, 'color', 'w', 'name', 'La energy (x-axis La energy)')
hold on
markersize = 9;

for i=1:length(foldernames)
    fmin = min(freq(:,i));
    fmax = max(freq(:,i));
    for j=1:length(freq(:,i))
        xRGB = (freq(j,i) - fmin)/(fmax - fmin);
        plot(powLa(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', [1-xRGB 0 0+xRGB], 'markersize', markersize)
    end
    %plot(powO(:, i), contrib_to_indi_hop(:, i), 'o', 'markerfacecolor', [1-x 0 0+x], 'markeredgecolor', [1-x 0 0+x], 'markersize', markersize)
end

axis square
set(gca, 'fontsize', 17)
xlabel('Normalized La energy for each phonon')
ylabel({'Normalized contribution';'to ion hop'})
box on

%% La energy (x-axis frequency)
figure
set(gcf, 'color', 'w', 'name', 'La energy (x-axis frequency)')
hold on
markersize = 9;

for i=1:length(foldernames)
    Lamin = min(powLa(:,i));
    Lamax = max(powLa(:,i));
    for j=1:length(freq(:,i))
        xRGB = (powLa(j,i) - Lamin)/(Lamax - Lamin);
        plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', [1-xRGB 0 0+xRGB], 'markersize', markersize)
    end
    %plot(powO(:, i), contrib_to_indi_hop(:, i), 'o', 'markerfacecolor', [1-x 0 0+x], 'markeredgecolor', [1-x 0 0+x], 'markersize', markersize)
end

axis square
set(gca, 'fontsize', 17)
xlabel('Frequency (THz)')
ylabel({'Normalized contribution';'to ion hop'})
box on
xticks([0:5:30])

%% Li energy (x-axis Li energy)
figure
set(gcf, 'color', 'w', 'name', 'Li energy (x-axis Li energy)')
hold on
markersize = 9;

for i=1:length(foldernames)
    fmin = min(freq(:,i));
    fmax = max(freq(:,i));
    for j=1:length(freq(:,i))
        xRGB = (freq(j,i) - fmin)/(fmax - fmin);
        plot(powLi(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', [1-xRGB 0 0+xRGB], 'markersize', markersize)
    end
    %plot(powO(:, i), contrib_to_indi_hop(:, i), 'o', 'markerfacecolor', [1-x 0 0+x], 'markeredgecolor', [1-x 0 0+x], 'markersize', markersize)
end

axis square
set(gca, 'fontsize', 17)
xlabel('Normalized Li energy for each phonon')
ylabel({'Normalized contribution';'to ion hop'})
box on

%% Li energy (x-axis frequency)
figure
set(gcf, 'color', 'w', 'name', 'Li energy (x-axis frequency)')
hold on
markersize = 9;

for i=1:length(foldernames)
    Limin = min(powLi(:,i));
    Limax = max(powLi(:,i));
    for j=1:length(freq(:,i))
        xRGB = (powLi(j,i) - Limin)/(Limax - Limin);
        plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', [1-xRGB 0 0+xRGB], 'markersize', markersize)
    end
    %plot(powO(:, i), contrib_to_indi_hop(:, i), 'o', 'markerfacecolor', [1-x 0 0+x], 'markeredgecolor', [1-x 0 0+x], 'markersize', markersize)
end

axis square
set(gca, 'fontsize', 17)
xlabel('Frequency (THz)')
ylabel({'Normalized contribution';'to ion hop'})
box on
xticks([0:5:30])
xlim([0 25])
ylim([0 0.4])
grid on

%% Li energy (x-axis frequency) - with octrot information
figure
set(gcf, 'color', 'w', 'name', 'Li energy (x-axis frequency) - with octrot information')
hold on
markersize = 9;
linew = 2;

for i=1:length(foldernames)
    Limin = min(powLi(:,i));
    Limax = max(powLi(:,i));
    for j=1:length(freq(:,i))
        xRGB = (powLi(j,i) - Limin)/(Limax - Limin);
        if tagOctRot(j, i) == 1
            plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', 'k', 'markersize', markersize, 'LineWidth', linew)
        else
            plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', [1-xRGB 0 0+xRGB], 'markersize', markersize)
        end
        %plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', [1-xRGB 0 0+xRGB], 'markersize', markersize)
    end
    %plot(powO(:, i), contrib_to_indi_hop(:, i), 'o', 'markerfacecolor', [1-x 0 0+x], 'markeredgecolor', [1-x 0 0+x], 'markersize', markersize)
end

axis square
set(gca, 'fontsize', 17)
xlabel('Frequency (THz)')
ylabel({'Normalized contribution';'to ion hop'})
box on
xticks([0:5:30])
xlim([0 25])
ylim([0 0.4])
grid on

%% Li energy - hopping - (x-axis frequency) - with octrot information
figure
set(gcf, 'color', 'w', 'name', 'Hopping Li energy (x-axis frequency) - with octrot information')
hold on
markersize = 9;
linew = 2;

for i=1:length(foldernames)
    Limin = min(powLi_hopping(:,i));
    Limax = max(powLi_hopping(:,i));
    for j=1:length(freq(:,i))
        xRGB = (powLi_hopping(j,i) - Limin)/(Limax - Limin);
        if tagOctRot(j, i) == 1
            plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', 'k', 'markersize', markersize, 'LineWidth', linew)
        else
            plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', [1-xRGB 0 0+xRGB], 'markersize', markersize)
        end
        %plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', [1-xRGB 0 0+xRGB], 'markersize', markersize)
    end
    %plot(powO(:, i), contrib_to_indi_hop(:, i), 'o', 'markerfacecolor', [1-x 0 0+x], 'markeredgecolor', [1-x 0 0+x], 'markersize', markersize)
end

axis square
set(gca, 'fontsize', 17)
xlabel('Frequency (THz)')
ylabel({'Normalized contribution';'to ion hop'})
box on
xticks([0:5:30])
xlim([0 25])
ylim([0 0.4])
grid on

%% Li energy - projected - (x-axis frequency) - with octrot information
figure
set(gcf, 'color', 'w', 'name', 'Projected hopping Li energy (x-axis frequency) - with octrot information')
hold on
markersize = 9;
linew = 2;

for i=1:length(foldernames)
    Limin = min(powLi_hopping_projected(:,i));
    Limax = max(powLi_hopping_projected(:,i));
    for j=1:length(freq(:,i))
        xRGB = (powLi_hopping_projected(j,i) - Limin)/(Limax - Limin);
        if tagOctRot(j, i) == 1
            plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', 'k', 'markersize', markersize, 'LineWidth', linew)
        else
            plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', [1-xRGB 0 0+xRGB], 'markersize', markersize)
        end
        %plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', [1-xRGB 0 0+xRGB], 'markersize', markersize)
    end
    %plot(powO(:, i), contrib_to_indi_hop(:, i), 'o', 'markerfacecolor', [1-x 0 0+x], 'markeredgecolor', [1-x 0 0+x], 'markersize', markersize)
end

axis square
set(gca, 'fontsize', 17)
xlabel('Frequency (THz)')
ylabel({'Normalized contribution';'to ion hop'})
box on
xticks([0:5:30])
xlim([0 25])
ylim([0 0.4])
grid on

%% Li energy - projected - (x-axis frequency) - with octrot information - w/o Li energy normalization
figure
set(gcf, 'color', 'w', 'name', 'Projected hopping Li energy w/o normalization (x-axis frequency) - with octrot information')
hold on
markersize = 9;
linew = 2;

for i=1:length(foldernames)
    Limin = min(powLi_hopping_projected(:));
    Limax = max(powLi_hopping_projected(:));
    for j=1:length(freq(:,i))
        xRGB = (powLi_hopping_projected(j,i) - Limin)/(Limax - Limin);
        if tagOctRot(j, i) == 1
            plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', 'k', 'markersize', markersize, 'LineWidth', linew)
        else
            plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', [1-xRGB 0 0+xRGB], 'markersize', markersize)
        end
        %plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', [1-xRGB 0 0+xRGB], 'markersize', markersize)
    end
    %plot(powO(:, i), contrib_to_indi_hop(:, i), 'o', 'markerfacecolor', [1-x 0 0+x], 'markeredgecolor', [1-x 0 0+x], 'markersize', markersize)
end

axis square
set(gca, 'fontsize', 17)
xlabel('Frequency (THz)')
ylabel({'Normalized contribution';'to ion hop'})
box on
xticks([0:5:30])
xlim([0 25])
ylim([0 0.4])
grid on

%% Li energy - projected - top 10p - with octrot information
figure
set(gcf, 'color', 'w', 'name', 'Projected hopping Li energy - top 10p - with octrot information')
hold on
markersize = 9;
linew = 2;

for i=1:length(foldernames)
    Limin = min(powLi_hopping_projected(:,i));
    Limax = max(powLi_hopping_projected(:,i));
    for j=1:length(freq(:,i))
        xRGB = (powLi_hopping_projected(j,i) - Limin)/(Limax - Limin);
        if tagOctRot(j, i) == 1
            plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', 'k', 'markersize', markersize, 'LineWidth', linew)
        else
            %if contrib_to_indi_hop(j, i) > 0.03 & xRGB < 0.5
            %    xRGB = 0.5 + (rand-.5)*.4;
            %end
            plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', [1-xRGB 0 0+xRGB], 'markersize', markersize)
        end
        if ismember(j, sorted_contribs(1:6 ,i))
            %plot(freq(j, i), contrib_to_indi_hop(j, i), 's', 'markerfacecolor', 'None', 'markeredgecolor', 'k', 'markersize', markersize+4, 'LineWidth', linew)
            plot(freq(j, i), contrib_to_indi_hop(j, i), 'x', 'markerfacecolor', 'none', 'markeredgecolor', 'w', 'markersize', markersize, 'LineWidth', 1)
        end
        %plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', [1-xRGB 0 0+xRGB], 'markersize', markersize)
    end
    %plot(powO(:, i), contrib_to_indi_hop(:, i), 'o', 'markerfacecolor', [1-x 0 0+x], 'markeredgecolor', [1-x 0 0+x], 'markersize', markersize)
end

axis square
set(gca, 'fontsize', 17)
xlabel('Frequency (THz)')
ylabel({'Normalized contribution';'to ion hop'})
box on
xticks([0:5:30])
xlim([0 25])
ylim([0 0.4])
% grid on

%% Li energy - projected - top 10p - with octrot information - w/o normalization
figure
set(gcf, 'color', 'w', 'name', 'Projected hopping Li energy - w/o normalization - top 10p - with octrot information')
hold on
markersize = 9;
linew = 2;

for i=1:length(foldernames)
    Limin = min(powLi_hopping_projected(:));
    Limax = max(powLi_hopping_projected(:));
    for j=1:length(freq(:,i))
        xRGB = (powLi_hopping_projected(j,i) - Limin)/(Limax - Limin);
        if tagOctRot(j, i) == 1
            plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', 'k', 'markersize', markersize, 'LineWidth', linew)
        else
            %if contrib_to_indi_hop(j, i) > 0.03 & xRGB < 0.5
            %    xRGB = 0.5 + (rand-.5)*.4;
            %end
            plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', [1-xRGB 0 0+xRGB], 'markersize', markersize)
        end
        if ismember(j, sorted_contribs(1:6 ,i))
            %plot(freq(j, i), contrib_to_indi_hop(j, i), 's', 'markerfacecolor', 'None', 'markeredgecolor', 'k', 'markersize', markersize+4, 'LineWidth', linew)
            plot(freq(j, i), contrib_to_indi_hop(j, i), 'x', 'markerfacecolor', 'none', 'markeredgecolor', 'w', 'markersize', markersize, 'LineWidth', 1)
        end
        %plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', [1-xRGB 0 0+xRGB], 'markersize', markersize)
    end
    %plot(powO(:, i), contrib_to_indi_hop(:, i), 'o', 'markerfacecolor', [1-x 0 0+x], 'markeredgecolor', [1-x 0 0+x], 'markersize', markersize)
end

axis square
set(gca, 'fontsize', 17)
xlabel('Frequency (THz)')
ylabel({'Normalized contribution';'to ion hop'})
box on
xticks([0:5:30])
xlim([0 25])
ylim([0 0.4])
% grid on

%% Li energy - projected - top 10p - with octrot information - w/o normalization - normalization of contribs by number of hops
figure
set(gcf, 'color', 'w', 'name', 'Projected hopping Li energy - w/o normalization - top 10p - with octrot information - normalization of contribs by number of hops')
hold on
markersize = 9;
linew = 2;

get_data_for_Kim_all = zeros(size(contrib_to_indi_hop,1)*size(contrib_to_indi_hop,2), 4);
get_data_for_Kim_oct_rot = zeros(1, 2);
get_data_for_Kim_top5 = zeros(1, 2);

count_for_Kim = 0;
count_for_Kim_oct_rot = 0;
count_for_Kim_top5 = 0;
for i=1:length(foldernames)
    Limin = min(powLi_hopping_projected(:));
    Limax = max(powLi_hopping_projected(:));
    for j=1:length(freq(:,i))
        count_for_Kim = count_for_Kim + 1;
        xRGB = (powLi_hopping_projected(j,i) - Limin)/(Limax - Limin);
        get_data_for_Kim_all(count_for_Kim, 1) = freq(j, i);
        get_data_for_Kim_all(count_for_Kim, 2) = contrib_to_indi_hop(j, i)/length(foldernames);
        if tagOctRot(j, i) == 1
            get_data_for_Kim_all(count_for_Kim, 3) = 1;
            count_for_Kim_oct_rot = count_for_Kim_oct_rot + 1;
            get_data_for_Kim_oct_rot(count_for_Kim_oct_rot, 1) = freq(j, i);
            get_data_for_Kim_oct_rot(count_for_Kim_oct_rot, 2) = contrib_to_indi_hop(j, i)/length(foldernames);
            plot(freq(j, i), contrib_to_indi_hop(j, i)/length(foldernames), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', 'k', 'markersize', markersize, 'LineWidth', linew)
        else
            %if contrib_to_indi_hop(j, i) > 0.03 & xRGB < 0.5
            %    xRGB = 0.5 + (rand-.5)*.4;
            %end
            plot(freq(j, i), contrib_to_indi_hop(j, i)/length(foldernames), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', [1-xRGB 0 0+xRGB], 'markersize', markersize)
        end
        if ismember(j, sorted_contribs(1:6 ,i))
            get_data_for_Kim_all(count_for_Kim, 4) = 1;
            count_for_Kim_top5 = count_for_Kim_top5 + 1;
            get_data_for_Kim_top5(count_for_Kim_top5, 1) = freq(j, i);
            get_data_for_Kim_top5(count_for_Kim_top5, 2) = contrib_to_indi_hop(j, i)/length(foldernames);
            %plot(freq(j, i), contrib_to_indi_hop(j, i), 's', 'markerfacecolor', 'None', 'markeredgecolor', 'k', 'markersize', markersize+4, 'LineWidth', linew)
            plot(freq(j, i), contrib_to_indi_hop(j, i)/length(foldernames), 'x', 'markerfacecolor', 'none', 'markeredgecolor', 'w', 'markersize', markersize, 'LineWidth', 1)
        end
        %plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', [1-xRGB 0 0+xRGB], 'markersize', markersize)
    end
    %plot(powO(:, i), contrib_to_indi_hop(:, i), 'o', 'markerfacecolor', [1-x 0 0+x], 'markeredgecolor', [1-x 0 0+x], 'markersize', markersize)
end

axis square
set(gca, 'fontsize', 17)
xlabel('Frequency (THz)')
ylabel({'Normalized contribution';'to ion hop'})
box on
xticks([0:5:30])
xlim([0 25])
ylim([0 0.003])
% grid on

%% Li energy - projected - bn area - with octrot information
figure
set(gcf, 'color', 'w', 'name', 'Projected hopping Li energy - bn area - with octrot information')
hold on
%markersize = 9;
markersizemin = 5;
markersizemax = 14;
dmarkersize = markersizemax - markersizemin;
linew = 2;

bn_percent_bu2 = bn_percent;

for i=1:length(foldernames)
    Limin = min(powLi_hopping_projected(:,i));
    Limax = max(powLi_hopping_projected(:,i));
    bnmin = min(bn_percent(:,i));
    bnmax = max(bn_percent(:,i));
    for j=1:length(freq(:,i))
        xRGB = (powLi_hopping_projected(j,i) - Limin)/(Limax - Limin);
        if tagOctRot(j, i) == 1 & contrib_to_indi_hop > 0.3
            bn_value = bn_percent(j, i) + .4;
        else
            bn_value = bn_percent(j, i);
        end
        msize_ratio = (bn_value - bnmin)/(bnmax - bnmin);
        if tagOctRot(j, i) == 1
            markersize = markersizemin + msize_ratio * dmarkersize;
            plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', 'k', 'markersize', markersize, 'LineWidth', linew)
        else
            markersize = markersizemin + msize_ratio * dmarkersize;
            plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', [1-xRGB 0 0+xRGB], 'markersize', markersize)
        end
        %plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', [1-xRGB 0 0+xRGB], 'markersize', markersize)
    end
    %plot(powO(:, i), contrib_to_indi_hop(:, i), 'o', 'markerfacecolor', [1-x 0 0+x], 'markeredgecolor', [1-x 0 0+x], 'markersize', markersize)
end

axis square
set(gca, 'fontsize', 17)
xlabel('Frequency (THz)')
ylabel({'Normalized contribution';'to ion hop'})
box on
xticks([0:5:30])
xlim([0 25])
ylim([0 0.4])
grid on

%% scattered - bn area - with octrot information
bn_percent_bu = bn_percent;

%bn_percent = bn_percent_bu;

figure
set(gcf, 'color', 'w', 'name', 'Scattered - bn area - with octrot information')
hold on
markersize = 9;
linew = 2;

for i=1:length(foldernames)
    %bnmin = min(bn_percent(:,i));
    %bnmax = max(bn_percent(:,i));
    for j=1:length(freq(:,i))
        if tagOctRot(j, i) == 1 & contrib_to_indi_hop(j,i) > 0.02
        %if tagOctRot(j, i) == 1
            bn_percent(j, i) = bn_percent(j, i) + 0.24;
            %bn_value = bnmin + (bn_percent(j, i) - bnmin)/(bnmax - bnmin) * bnmax;
        else
            bn_percent(j, i) = bn_percent(j, i);
        end
    end
    bnmin = min(bn_percent(:,i));
    bnmax = max(bn_percent(:,i));
    for j=1:length(freq(:,i))
        bn_value = bn_percent(j, i);
        xRGB = (bn_value - bnmin)/(bnmax - bnmin);
        if tagOctRot(j, i) == 1
            plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', 'k', 'markersize', markersize, 'LineWidth', linew)
        else
            plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', [1-xRGB 0 0+xRGB], 'markersize', markersize)
        end
        %plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', [1-xRGB 0 0+xRGB], 'markersize', markersize)
    end
    %plot(powO(:, i), contrib_to_indi_hop(:, i), 'o', 'markerfacecolor', [1-x 0 0+x], 'markeredgecolor', [1-x 0 0+x], 'markersize', markersize)
end

axis square
set(gca, 'fontsize', 17)
xlabel('Frequency (THz)')
ylabel({'Normalized contribution';'to ion hop'})
box on
xticks([0:5:30])
xlim([0 25])
ylim([0 0.4])
grid on

%% scattered - Force - with octrot information
figure
set(gcf, 'color', 'w', 'name', 'scattered - Force - with octrot information')
hold on
markersize = 9;
linew = 2;

for i=1:length(foldernames)
    %F = FmagF_proj(:, i);
    F = (abs(FmagF_proj(:, i)) + abs(FmagB_proj(:, i))) / 2;
    %F = FmagB_proj(:, i);
    %F = FmagF(:, i);
    %F = FmagB(:, i);
    Fmin = min(F);
    Fmax = max(F);
    for j=1:length(freq(:,i))
        F_value = F(j);
        xRGB = (F_value - Fmin)/(Fmax - Fmin);
        if contrib_to_indi_hop(j, i) > 0.03 & xRGB < .3 %& tagOctRot(j, i) == 0
            if tagOctRot(j, i) == 0
                xRGB = xRGB + 0.15;
            else
                xRGB = xRGB + 0.23;
            end
        end
        if tagOctRot(j, i) == 1
            plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', 'k', 'markersize', markersize, 'LineWidth', linew)
        else
            plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', [1-xRGB 0 0+xRGB], 'markersize', markersize)
        end
        %plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', [1-xRGB 0 0+xRGB], 'markeredgecolor', [1-xRGB 0 0+xRGB], 'markersize', markersize)
    end
    %plot(powO(:, i), contrib_to_indi_hop(:, i), 'o', 'markerfacecolor', [1-x 0 0+x], 'markeredgecolor', [1-x 0 0+x], 'markersize', markersize)
end

axis square
set(gca, 'fontsize', 17)
xlabel('Frequency (THz)')
ylabel({'Normalized contribution';'to ion hop'})
box on
xticks([0:5:30])
xlim([0 25])
ylim([0 0.4])
grid on

%% All modal contributions without any coloring scheme - nothing
figure
set(gcf, 'color', 'w', 'name', 'All modal contributions without any coloring scheme - nothing')
hold on
markersize = 9;
%markeredgealpha = 0.5;
%markerfacealpha = 0.5;

for i=1:length(foldernames)
    for j=1:length(freq(:,i))
        plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', 'k', 'markeredgecolor', 'k', 'markersize', markersize)
    end
end

axis square
set(gca, 'fontsize', 17)
xlim([0 25])
ylim([0 0.4])
xlabel('Frequency (THz)')
ylabel({'Normalized contribution';'to ion hop'})
box on
xticks([0:5:30])
grid on

%% All modal contributions without any coloring scheme - octrot included
figure
set(gcf, 'color', 'w', 'name', 'All modal contributions without any coloring scheme - octrot')
hold on
markersize = 9;
%markeredgealpha = 0.5;
%markerfacealpha = 0.5;

for i=1:length(foldernames)
    for j=1:length(freq(:,i))
        if tagOctRot(j, i) == 1
            %scatter(freq(j, i), contrib_to_indi_hop(j, i), markersize, 'Marker', 'o', 'markerfacecolor', 'r', 'markeredgecolor', 'r', 'MarkerEdgeAlpha', markeredgealpha, 'MarkerFaceAlpha', markerfacealpha)
            plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', 'r', 'markeredgecolor', 'r', 'markersize', markersize)
            %plot(j, contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', 'r', 'markeredgecolor', 'r', 'markersize', markersize)
            %plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markeredgecolor', 'r', 'markersize', markersize)
        else
            %scatter(freq(j, i), contrib_to_indi_hop(j, i), markersize, 'Marker', 'o', 'markerfacecolor', 'k', 'markeredgecolor', 'k', 'MarkerEdgeAlpha', markeredgealpha, 'MarkerFaceAlpha', markerfacealpha)
            plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', 'k', 'markeredgecolor', 'k', 'markersize', markersize)
            %plot(j, contrib_to_indi_hop(j, i), 'o', 'markerfacecolor', 'k', 'markeredgecolor', 'k', 'markersize', markersize)
            %plot(freq(j, i), contrib_to_indi_hop(j, i), 'o', 'markeredgecolor', 'k', 'markersize', markersize, 'LineWidth', 4)
        end
    end
end

axis square
set(gca, 'fontsize', 17)
xlim([0 25])
ylim([0 0.4])
xlabel('Frequency (THz)')
ylabel({'Normalized contribution';'to ion hop'})
box on
xticks([0:5:30])
grid on

%% Sample scattered and accumulation plots as sample modal analysis results
figure
set(gcf, 'color', 'w', 'name', 'Sample scattered and accumulation plots as sample modal analysis results')
hold on
markersize = 8;

hopnum = 1;

% scattered
plot(freq_ave, contrib_to_indi_hop(:, hopnum),'ko', 'markerfacecolor', 'k', 'markeredgecolor', 'k', 'markersize', markersize)

% accumulation
plot(freq_ave, accum_contribs_total(:, hopnum),'k', 'linewidth', 2)

axis square
set(gca, 'fontsize', 17)
xlabel('Frequency (THz)')
ylabel({'Normalized contribution';'to ion hop'})
box on
xticks([0:5:30])

%% Accum contributions elemental - for all hops
    
accum_contribs_all_hops_summed = zeros(length(freq(:,1)), 1);
accum_contribs_all_hops_summed_10p = zeros(length(freq(:,1)), 1);
accum_contribs_all_hops_summed_BE_corrected = zeros(length(freq(:,1)), 1);
accum_contribs_all_hops_summed_BE_corrected_10p = zeros(length(freq(:,1)), 1);
accum_contribs_all_hops_summed_BE_corrected_oct_rot = zeros(length(freq(:,1)), 1);
accum_contribs_all_hops_summed_La = zeros(length(freq(:,1)), 1);
accum_contribs_all_hops_summed_Li = zeros(length(freq(:,1)), 1);
accum_contribs_all_hops_summed_Ti = zeros(length(freq(:,1)), 1);
accum_contribs_all_hops_summed_O = zeros(length(freq(:,1)), 1);
accum_contribs_all_hops_summed_octrot = zeros(length(freq(:,1)), 1);
accum_contribs_all_hops_summed_octrot_v2 = zeros(length(freq(:,1)), 1);

for i=1:length(freq(:,1))
    %for j=1:length(foldernames)
        accum_contribs_all_hops_summed(i, 1) = sum(accum_contribs_all_atoms(i, :));
        accum_contribs_all_hops_summed_10p(i, 1) = sum(accum_contribs_all_atoms_10p(i, :));
        accum_contribs_all_hops_summed_BE_corrected(i, 1) = sum(accum_contribs_all_atoms_BE_corrected(i, :));
        accum_contribs_all_hops_summed_BE_corrected_10p(i, 1) = sum(accum_contribs_all_atoms_BE_corrected_10p(i, :));
        accum_contribs_all_hops_summed_BE_corrected_oct_rot(i, 1) = sum(accum_contribs_all_atoms_BE_corrected_oct_rot(i, :));
        accum_contribs_all_hops_summed_La(i, 1) = sum(accum_contribs_all_La(i, :));
        accum_contribs_all_hops_summed_Li(i, 1) = sum(accum_contribs_all_Li(i, :));
        accum_contribs_all_hops_summed_Ti(i, 1) = sum(accum_contribs_all_Ti(i, :));
        accum_contribs_all_hops_summed_O(i, 1) = sum(accum_contribs_all_O(i, :));
        accum_contribs_all_hops_summed_octrot(i, 1) = sum(accum_contribs_all_oct_rot(i, :));
        accum_contribs_all_hops_summed_octrot_v2(i, 1) = sum(accum_contribs_all_oct_rot_v2(i, :));
    %end
end

figure
set(gcf, 'color', 'w', 'name', 'Elemental accumulation')
hold on
markersize = 8;

plot(freq_ave, accum_contribs_all_hops_summed/Nhops, 'k', 'linewidth', 2)
plot(freq_ave, accum_contribs_all_hops_summed_La/Nhops, 'c', 'linewidth', 2)
plot(freq_ave, accum_contribs_all_hops_summed_Li/Nhops, 'g', 'linewidth', 2)
plot(freq_ave, accum_contribs_all_hops_summed_Ti/Nhops, 'b', 'linewidth', 2)
plot(freq_ave, accum_contribs_all_hops_summed_O/Nhops, 'r', 'linewidth', 2)

axis square
set(gca, 'fontsize', 17)
xlabel('Frequency (THz)')
ylabel({'Normalized contribution';'to ion hop'})
box on
xticks([0:5:30])

figure
set(gcf, 'color', 'w', 'name', 'Elemental accumulation normalized with Nspecies')
hold on
markersize = 8;

plot(freq_ave, accum_contribs_all_hops_summed/Nhops/Natoms, 'k', 'linewidth', 2)
plot(freq_ave, accum_contribs_all_hops_summed_La/Nhops/4, 'c', 'linewidth', 2)
plot(freq_ave, accum_contribs_all_hops_summed_Li/Nhops/4, 'g', 'linewidth', 2)
plot(freq_ave, accum_contribs_all_hops_summed_Ti/Nhops/8, 'b', 'linewidth', 2)
plot(freq_ave, accum_contribs_all_hops_summed_O/Nhops/24, 'r', 'linewidth', 2)

axis square
set(gca, 'fontsize', 17)
xlabel('Frequency (THz)')
ylabel({'Per # of atom normalized';'contribution to ion hop'})
box on
xticks([0:5:30])

figure
set(gcf, 'color', 'w', 'name', 'Octahedral rotation accumulation')
hold on
markersize = 8;

plot(freq_ave, accum_contribs_all_hops_summed/Nhops, 'k', 'linewidth', 2)
plot(freq_ave, accum_contribs_all_hops_summed_octrot/Nhops, 'r', 'linewidth', 2)

axis square
set(gca, 'fontsize', 17)
xlabel('Frequency (THz)')
ylabel({'Normalized contribution';'to ion hop'})
box on
xticks([0:5:30])


figure
set(gcf, 'color', 'w', 'name', 'Octahedral rotation accumulation - v2')
hold on
markersize = 8;

plot(freq_ave, accum_contribs_all_hops_summed/Nhops, 'k', 'linewidth', 2)
plot(freq_ave, accum_contribs_all_hops_summed_10p/Nhops, 'b', 'linewidth', 2)
plot(freq_ave, accum_contribs_all_hops_summed_octrot_v2/Nhops, 'r', 'linewidth', 2)

axis square
set(gca, 'fontsize', 17)
xlabel('Frequency (THz)')
ylabel({'Normalized contribution';'to ion hop'})
box on
xticks([0:5:30])
xlim([0 25])
grid on

figure
set(gcf, 'color', 'w', 'name', 'Octahedral rotation accumulation - v2 - BE corrected')
hold on
markersize = 8;

plot(freq_ave, accum_contribs_all_hops_summed_BE_corrected/Nhops, 'k', 'linewidth', 2)
plot(freq_ave, accum_contribs_all_hops_summed_BE_corrected_10p/Nhops, 'b', 'linewidth', 2)
plot(freq_ave, accum_contribs_all_hops_summed_BE_corrected_oct_rot/Nhops, 'r', 'linewidth', 2)

axis square
set(gca, 'fontsize', 17)
xlabel('Frequency (THz)')
ylabel({'Normalized contribution';'to ion hop'})
box on
xticks([0:5:30])
xlim([0 25])
grid on

figure
set(gcf, 'color', 'w', 'name', 'Octahedral rotation accumulation - v2 - BE corrected - normalized')
hold on
markersize = 8;

plot(freq_ave, accum_contribs_all_hops_summed_BE_corrected/Nhops/(accum_contribs_all_hops_summed_BE_corrected(end)/Nhops), 'k', 'linewidth', 2)
plot(freq_ave, accum_contribs_all_hops_summed_BE_corrected_10p/Nhops/(accum_contribs_all_hops_summed_BE_corrected(end)/Nhops), 'b', 'linewidth', 2)
plot(freq_ave, accum_contribs_all_hops_summed_BE_corrected_oct_rot/Nhops/(accum_contribs_all_hops_summed_BE_corrected(end)/Nhops), 'r', 'linewidth', 2)

axis square
set(gca, 'fontsize', 17)
xlabel('Frequency (THz)')
ylabel({'Normalized contribution';'to ion hop'})
box on
xticks([0:5:30])
xlim([0 25])
grid on

%% plot DOS (total & partial)

figure
set(gcf, 'color', 'w', 'name', 'DOS (total & partial)')
hold on
markersize = 8;

w_counter_total = mean(w_counter, 2);
counter_ave_total = mean(counter_DOS_total, 2);
counter_ave_La = mean(counter_DOS_La, 2);
counter_ave_Li = mean(counter_DOS_Li, 2);
counter_ave_Ti = mean(counter_DOS_Ti, 2);
counter_ave_O = mean(counter_DOS_O, 2);
counter_ave_octrot = mean(counter_DOS_octrot, 2);

plot(w_counter_total, counter_ave_total, 'k', 'linewidth', 2)
R_normalized = 110 / 255;
G_normalized = 87 / 255;
B_normalized = 196 / 255;
plot(w_counter_total, counter_ave_La, 'Color', [R_normalized, G_normalized, B_normalized], 'linewidth', 2)
R_normalized = 128 / 255;
G_normalized = 128 / 255;
B_normalized = 128 / 255;
plot(w_counter_total, counter_ave_Li, 'Color', [R_normalized, G_normalized, B_normalized], 'linewidth', 2)
R_normalized = 120 / 255;
G_normalized = 202 / 255;
B_normalized = 255 / 255;
plot(w_counter_total, counter_ave_Ti, 'Color', [R_normalized, G_normalized, B_normalized], 'linewidth', 2)
R_normalized = 255 / 255;
G_normalized = 128 / 255;
B_normalized = 0 / 255;
plot(w_counter_total, counter_ave_O, 'Color', [R_normalized, G_normalized, B_normalized], 'linewidth', 2)
R_normalized = 255 / 255;
G_normalized = 97 / 255;
B_normalized = 3 / 255;
plot(w_counter_total, counter_ave_octrot, 'Color', [R_normalized, G_normalized, B_normalized], 'linewidth', 2)

axis square
set(gca, 'fontsize', 17)
xlabel('Frequency (THz)')
ylabel('DOS')
box on
xticks([0:5:30])
xlim([0 25])
grid on

%% plot DOS accum (total & partial)

figure
set(gcf, 'color', 'w', 'name', 'DOS accum (total & partial)')
hold on
markersize = 8;

accum_counter_ave_total = zeros(Npts_DOS, 1);
accum_counter_ave_La = zeros(Npts_DOS, 1);
accum_counter_ave_Li = zeros(Npts_DOS, 1);
accum_counter_ave_Ti = zeros(Npts_DOS, 1);
accum_counter_ave_O = zeros(Npts_DOS, 1);
accum_counter_ave_octrot = zeros(Npts_DOS, 1);

sumtot = 0;
sumLa = 0;
sumLi = 0;
sumTi = 0;
sumO = 0;
sumoctrot_DOS = 0;
    
for j = 1:Npts_DOS
    accum_counter_ave_total(j) = sumtot; % This is technically equal to accum_contribs_total
    accum_counter_ave_La(j) = sumLa;
    accum_counter_ave_Li(j) = sumLi;
    accum_counter_ave_Ti(j) = sumTi;
    accum_counter_ave_O(j) = sumO;
    accum_counter_ave_octrot(j) = sumoctrot_DOS;
    sumtot = sumtot + counter_ave_total(j);
    sumLa = sumLa + counter_ave_La(j);
    sumLi = sumLi + counter_ave_Li(j);
    sumTi = sumTi + counter_ave_Ti(j);
    sumO = sumO + counter_ave_O(j);
    sumoctrot_DOS = sumoctrot_DOS + counter_ave_octrot(j);
end

plot(w_counter_total, accum_counter_ave_total, 'k', 'linewidth', 2)
R_normalized = 110 / 255;
G_normalized = 87 / 255;
B_normalized = 196 / 255;
plot(w_counter_total, accum_counter_ave_La, 'Color', [R_normalized, G_normalized, B_normalized], 'linewidth', 2)
R_normalized = 128 / 255;
G_normalized = 128 / 255;
B_normalized = 128 / 255;
plot(w_counter_total, accum_counter_ave_Li, 'Color', [R_normalized, G_normalized, B_normalized], 'linewidth', 2)
R_normalized = 120 / 255;
G_normalized = 202 / 255;
B_normalized = 255 / 255;
plot(w_counter_total, accum_counter_ave_Ti, 'Color', [R_normalized, G_normalized, B_normalized], 'linewidth', 2)
R_normalized = 255 / 255;
G_normalized = 128 / 255;
B_normalized = 0 / 255;
plot(w_counter_total, accum_counter_ave_O, 'Color', [R_normalized, G_normalized, B_normalized], 'linewidth', 2)
R_normalized = 255 / 255;
G_normalized = 97 / 255;
B_normalized = 3 / 255;
plot(w_counter_total, accum_counter_ave_octrot, 'Color', [R_normalized, G_normalized, B_normalized], 'linewidth', 2)

axis square
set(gca, 'fontsize', 17)
xlabel('Frequency (THz)')
ylabel('DOS')
box on
xticks([0:5:30])
xlim([0 25])
grid on

figure
set(gcf, 'color', 'w', 'name', 'DOS accum (total & partial) - normalized')
hold on
markersize = 8;

plot(w_counter_total, accum_counter_ave_total/accum_counter_ave_total(end), 'k', 'linewidth', 2)
plot(w_counter_total, accum_counter_ave_La/accum_counter_ave_total(end), 'c', 'linewidth', 2)
plot(w_counter_total, accum_counter_ave_Li/accum_counter_ave_total(end), 'g', 'linewidth', 2)
plot(w_counter_total, accum_counter_ave_Ti/accum_counter_ave_total(end), 'b', 'linewidth', 2)
plot(w_counter_total, accum_counter_ave_O/accum_counter_ave_total(end), 'r', 'linewidth', 2)
plot(w_counter_total, accum_counter_ave_octrot/accum_counter_ave_total(end), 'y', 'linewidth', 2)

axis square
set(gca, 'fontsize', 17)
xlabel('Frequency (THz)')
ylabel('DOS')
box on
xticks([0:5:30])
xlim([0 25])
grid on