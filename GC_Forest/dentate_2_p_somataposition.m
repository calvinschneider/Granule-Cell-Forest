%% Determines Whether Somata are Deep/Superficial and Infra/Suprapyramidal and Discards Those Outside GCL
% input ranges from 1 to 257 for current GCL

function dentate_2_p_somataposition(input,directory)

% Load somata and define GCL points
load('./Outputs/Somata_all.mat')
[x_o,y_o,z_o]   = layer_eq_GCL_2(0);
[x_i,y_i,z_i]   = layer_eq_GCL_2(-1.95);
M_o             = [x_o,y_o,z_o];
M_i             = [x_i,y_i,z_i];

% Define somata range based on input
subset_size = 5000;
input2      = str2double(input);
startsoma   = (input2-1)*subset_size + 1;
endsoma     = input2*subset_size;
if endsoma > length(Somata)
    endsoma = length(Somata);
end

% Split GCL points into 100 micron bins for speedup
xmax    = ceil(max(x_o)/100)*100;
xmin    = floor(min(x_o)/100)*100;
n_bins  = (xmax-xmin)/100;
GCL_o   = cell(1,n_bins);
GCL_i   = cell(1,n_bins);
for bin = 1:n_bins
    GCL_o{bin} = M_o(M_o(:,1)>=(xmin+(bin-1)*100)& M_o(:,1)<=(xmin+bin*100),:);
    GCL_i{bin} = M_i(M_i(:,1)>=(xmin+(bin-1)*100)& M_i(:,1)<=(xmin+bin*100),:);
end

% Define granule cell layer parameters from layer_eq_GCL
u_params   = [pi*1/100,pi*98/100,2000];
v_params   = [pi*-23/100,pi*142.5/100,1000];

% Define infra/suprapyramidal split as halfway around C-shape
supra_infra_border = v_params(1,1) + (v_params(1,2) - v_params(1,1))/2;

% Define somata radius (used to eliminate cells with portions outside GCL)
radius = 6.27;

% Allocate variables
Somata_pos      = zeros(((endsoma-startsoma)+1),6);
Somata_trimmed  = zeros(((endsoma-startsoma)+1),3);
counter         = 1;

for i = startsoma:endsoma
    % Limit GCL points tested to those near soma
    bin_number_soma     = ceil((Somata(i,1)-xmin)/100);
    bin_start_surface   = bin_number_soma - 3;
    bin_end_surface     = bin_number_soma + 3;
    if bin_start_surface < 1
        bin_start_surface = 1;
    end
    if bin_end_surface > n_bins
        bin_end_surface = n_bins;
    end
    GCL_o_current = vertcat(GCL_o{bin_start_surface:bin_end_surface});
    GCL_i_current = vertcat(GCL_i{bin_start_surface:bin_end_surface});
    
    % Find closest point on the inner and outer GCL
    [k_o,d_o]   = dsearchn(GCL_o_current,Somata(i,:));
    [k_i,d_i]   = dsearchn(GCL_i_current,Somata(i,:));
    index_o     = find(M_o(:,1) == GCL_o_current(k_o,1) & M_o(:,2) == GCL_o_current(k_o,2) & M_o(:,3) == GCL_o_current(k_o,3));
    index_i     = find(M_i(:,1) == GCL_i_current(k_i,1) & M_i(:,2) == GCL_i_current(k_i,2) & M_i(:,3) == GCL_i_current(k_i,3));
    
    % Find u and v coordinates from closest points
    u_bin_o     = ceil(index_o/v_params(1,3));
    u_bin_i     = ceil(index_i/v_params(1,3));
    u_o         = u_params(1,1) + (u_bin_o - 1) * ((u_params(1,2)-u_params(1,1))/(u_params(1,3)-1));
    %u_i         = u_params(1,1) + (u_bin_i - 1) * ((u_params(1,2)-u_params(1,1))/(u_params(1,3)-1));
    v_bin_o     = index_o - ((u_bin_o - 1) * v_params(1,3));
    v_bin_i     = index_i - ((u_bin_i - 1) * v_params(1,3));
    v_o         = v_params(1,1) + (v_bin_o - 1) * ((v_params(1,2)-v_params(1,1))/(v_params(1,3)-1));
    v_i         = v_params(1,1) + (v_bin_i - 1) * ((v_params(1,2)-v_params(1,1))/(v_params(1,3)-1));
    
    % Write out position and keep only if completely inside GCL
    if d_o > radius && d_i > radius
        if d_o > d_i
            % Write out deep
            Somata_pos(counter,1) = 0;
        else
            % Write out superficial
            Somata_pos(counter,1) = 1;
        end

        if v_i > supra_infra_border
            % Write out infrapyramidal
            Somata_pos(counter,2) = 0;
        else
            % Write out suprapyramidal
            Somata_pos(counter,2) = 1;
        end
        Somata_trimmed(counter,:)   = Somata(i,:);
        Somata_pos(counter,3)       = u_o;
        Somata_pos(counter,4)       = v_o;
        Somata_pos(counter,5)       = d_i;
        Somata_pos(counter,6)       = d_o;
        counter                     = counter + 1;
    end
end

Somata      = Somata_trimmed(1:(counter-1),:);
Somata_pos  = Somata_pos(1:(counter-1),:);

save(sprintf('%s/Soma_Position/%s.mat',directory,input),'Somata_pos','-v7.3')
save(sprintf('%s/Soma_Trimmed/%s.mat',directory,input),'Somata','-v7.3')
