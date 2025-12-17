% AUTHOR: Julien Schmidt
% DATE: Fall 2022

clear

method = 'default'; %RUN FOR PART A
%method = 'c'; %RUN FOR PART B : CHAINS
%method = 'r'; %RUN FOR PART C : RUNGS

if method == 'c'
    wall_width = 1;
    v_wall = 5;
    N_chain = 42;
    N_ring = 0;
    r_S = 0.1;
elseif method == 'r'
    wall_width = 1;
    v_wall = 0;
    N_chain = 0;
    N_ring = 10;
    r_S = 0.1;
else
    wall_width = 0;
    v_wall = 0;
    N_chain = 0;
    N_ring = 0;
    r_S = 0.3;
end

L = 15;
sigma = 1;
gamma = 4.5;
r_c = 1;
rho = 4;
K_S = 100;
kBT = sigma^2/(2*gamma);

%we save our particle positions, velocities and accelerations in these three structs
pos = struct();
vels = struct();
accs = struct();

N2 = (rho * L^2);
N = (rho * L^2) - (N_chain * 7) - (N_ring * 9);
particles = complex(rand(N,1)*L,rand(N,1)*L);
particlesWall = ((real(particles) >= 0 & real(particles) <= wall_width) | (real(particles) >= L-wall_width & real(particles) <= L));

%init fluid and wall positions
pos.fluid = particles(~particlesWall);
%pos.fluid = [];
pos.wall = particles(particlesWall);
pos.A = [];
pos.B = [];

%init chains and chain particle positions
chains = zeros(N_chain,7); %here this holds the index of particles A and B from pos.A and pos.B. chains(:,1:2) holds indicies of A and chains(:,3:7) holds indices of B
%build chains
idx_A = 1;
idx_B = 1;
for c=1:N_chain
    start_pos = complex(rand(1),rand(1))*L;
    angle_of_build = complex(rand(1),rand(1));
    angle_of_build = angle_of_build/abs(angle_of_build);
    dist_modifier = r_S;

    molecule_idxs = [];
    for a=0:1
        pos.A(idx_A,1) = start_pos + angle_of_build*dist_modifier*a;
        pos.A(idx_A,1) = complex(mod(real(pos.A(idx_A)),L),mod(imag(pos.A(idx_A)),L));
        molecule_idxs = [molecule_idxs; idx_A];
        idx_A = idx_A + 1;
    end
    for b=2:6
        pos.B(idx_B,1) = start_pos + angle_of_build*dist_modifier*b;
        pos.B(idx_B,1) = complex(mod(real(pos.B(idx_B)),L),mod(imag(pos.B(idx_B)),L));
        molecule_idxs = [molecule_idxs; idx_B];
        idx_B = idx_B + 1;
    end 
    
    chains(c,:) = molecule_idxs;
end

%init rings and ring particle positions
rings = zeros(N_ring,9);
for r=1:N_ring
    start_pos = complex(rand(1),rand(1))*L;
    thetas = linspace(0,2*pi,10);
    molecule_idxs = [];
    for a=1:9
        pos.A(idx_A,1) = start_pos + complex(sin(thetas(a)),cos(thetas(a)))/2.2;
        molecule_idxs = [molecule_idxs; idx_A];
        idx_A = idx_A + 1;
    end
    rings(r,:) = molecule_idxs;
end

%init velocities
vels.fluid = zeros(length(pos.fluid),1);
vels.wall = complex(zeros(length(pos.wall),1),ones(length(pos.wall),1)*v_wall);
vels.A = zeros(length(pos.A),1);
vels.B = zeros(length(pos.B),1);

%init accelerations
accs.fluid = zeros(length(pos.fluid),1);
accs.wall = zeros(length(pos.wall),1);
accs.A = zeros(length(pos.A),1);
accs.B = zeros(length(pos.B),1);

a_ij = [50 25 25 200; 25 1 300 200; 25 300 25 200; 200 200 200 0];

dt = 0.01;

saved_data = struct();
saved_data.Mnt = [];
saved_data.T = [];

f = figure(100);
f.Position = [0 0 1000 1000];

%uncomment needed sections to save as video if desired
%writerObj = VideoWriter('PartC');
%writerObj.FrameRate = 30;
%open(writerObj)

for step=1:1000 %number of iterations
    
    %update particle positions
    pos.wall(pos.wall < L/2) = pos.wall(pos.wall < L/2) + (vels.wall(pos.wall < L/2) .* dt);
    pos.wall(pos.wall >= L/2) = pos.wall(pos.wall >= L/2) - (vels.wall(pos.wall >= L/2) .* dt);
    pos.fluid = pos.fluid + (vels.fluid .* dt);
    pos.A = pos.A + (vels.A .* dt);
    pos.B = pos.B + (vels.B .* dt);
    

    %force potitions to be periodic
    pos.wall = complex(mod(real(pos.wall),L),mod(imag(pos.wall),L));
    pos.fluid = complex(mod(real(pos.fluid),L),mod(imag(pos.fluid),L));
    pos.A = complex(mod(real(pos.A),L),mod(imag(pos.A),L));
    pos.B = complex(mod(real(pos.B),L),mod(imag(pos.B),L));
    
    %update accelerations
    accs_prev = accs;
    accs = calculate_forces(pos,vels,L,r_c,gamma,sigma,a_ij);
    molecule_forces = calculate_molecule_forces(pos, chains, rings, r_S, K_S, L, method);
    if ~isempty(molecule_forces.chainsA) && ~isempty(molecule_forces.chainsB)
        accs.A = accs.A + molecule_forces.chainsA;
        accs.B = accs.B + molecule_forces.chainsB;
    end
    if ~isempty(molecule_forces.ringsA)
        accs.A = accs.A + molecule_forces.ringsA;
    end
    if method == 'r'
        accs.A = accs.A + complex(0,0.3);
        accs.B = accs.B + complex(0,0.3);
        accs.fluid = accs.fluid + complex(0,0.3);
    end

    %update velocities
    vels.fluid = vels.fluid + ((accs.fluid) .* dt);
    vels.A = vels.A + ((accs.A) .* dt);
    vels.B = vels.B + ((accs.B) .* dt);
    
    subplot(4,4,[2 3 6 7])
    cla
    hold on
    %plot(pos.fluid,'o','color','#d49817','MarkerFaceColor','#d49817','MarkerSize',4)
    %plot(pos.wall,'o','color','#323ea8','MarkerFaceColor','#323ea8','MarkerSize',4)
    %plot(pos.A,'o','color','#8c17d4','MarkerFaceColor','#8c17d4','MarkerSize',4)
    %plot(pos.B,'o','color','#1fab98','MarkerFaceColor','#1fab98','MarkerSize',4)
    plot(pos.fluid,'o','color','#d49817','MarkerFaceColor','#d49817')
    plot(pos.wall,'o','color','#323ea8','MarkerFaceColor','#323ea8')
    plot(pos.A,'o','color','#8c17d4','MarkerFaceColor','#8c17d4')
    plot(pos.B,'o','color','#1fab98','MarkerFaceColor','#1fab98')
    xlim([0,L])
    ylim([0,L])
    title('Simulation','fontsize',16)

    saved_data.Mnt = [saved_data.Mnt, abs(sum((vels.fluid)) + (sum((vels.wall))) + (sum((vels.A))) + (sum((vels.B))))];
    all_vels = [vels.fluid; vels.wall; vels.A; vels.B];
    saved_data.T = [saved_data.T, (sum(abs(all_vels).^2) / (2)) / N2];
    if mod(step,1) == 0

        %plotting momentum
        subplot(4,4,[9 12])
        cla
        hold on
        plot(saved_data.Mnt,'b-','LineWidth',2)
        %ylim([-1e-11,1e-11])
        title('Total Momentum','fontsize',16)

        %plotting thermostat
        subplot(4,4,[13 16])
        cla
        hold on
        plot(saved_data.T,'r-','LineWidth',2)
        xlabel('Iteration Number','FontSize',14)
        title('Temperature','fontsize',16)




    end

    drawnow
    %pause(0.01)
    %frame = getframe(gcf); %get frame
    %writeVideo(writerObj, frame);
    
end
%close(writerObj)

function accs = calculate_forces(pos,vels,L,r_c,gamma,sigma,a_ij)

    accs = struct();
    accs.fluid = zeros(length(pos.fluid),1);
    accs.wall = zeros(length(pos.wall),1);
    accs.A = zeros(length(pos.A),1);
    accs.B = zeros(length(pos.B),1);
    
    fn = fieldnames(pos);
    fnv = fieldnames(vels);
    for i=1:numel(fn)
        for j=i:numel(fn)
            %here we compare one set of particles to another set. i.e.
            %fluid to A.. or fluid to walls.. or walls to B..
            posi = pos.(fn{i}); %positions of one set of particles
            posj = pos.(fn{j}); %positions of other set of particles
            if ~isempty(posi) && ~isempty(posj)
                
                velsi = vels.(fnv{i});
                velsj = vels.(fnv{j});
                
                ifs = [1,2,3,4]; %below, ii and jj are retured as the correct indicies for a_ij entry
                ii = ifs([strcmp(fn{i},'A'), strcmp(fn{i},'B'), strcmp(fn{i},'fluid'), strcmp(fn{i},'wall')]);
                jj = ifs([strcmp(fn{j},'A'), strcmp(fn{j},'B'), strcmp(fn{j},'fluid'), strcmp(fn{j},'wall')]);
                aa = a_ij(ii,jj); %get the proper a_ij entry for these two partcile types interaction
    
                dists = posi - posj.'; %get distance and apply periodic accounting
                dists(real(dists) > L/2) = dists(real(dists) > L/2) - complex(L,0);
                dists(real(dists) < -L/2) = dists(real(dists) < -L/2) + complex(L,0);
                dists(imag(dists) > L/2) = dists(imag(dists) > L/2) - complex(0,L);
                dists(imag(dists) < -L/2) = dists(imag(dists) < -L/2) + complex(0,L);

                r = abs(dists);
                r_hat = dists ./ r;
    
                v = velsi - velsj.';
    
                w_R = zeros(length(posi),length(posj));
                w_R(r<r_c) = 1 - r(r<r_c)./r_c;
    
                w_D = zeros(length(posi),length(posj));
                w_D(r<r_c) = w_R(r<r_c).^2;
                
                epsilon = randn(length(posi),length(posj));
                %epsilon = ones(length(posi),length(posj));
    
                F_C = zeros(length(posi),length(posj));
                F_C(r<r_c) = aa.*(1-r(r<r_c)./r_c).*r_hat(r<r_c); %this may not work.. im not sure
    
                F_D = zeros(length(posi),length(posj));
                F_D(r<r_c) = -gamma*w_D(r<r_c) .* (real(r_hat(r<r_c)) .* real(v(r<r_c)) + imag(r_hat(r<r_c)) .* imag(v(r<r_c))) .* r_hat(r<r_c);      %this may not work either..... remember v(r<r_c)=0 forall at first iteration 
    
                F_R = zeros(length(posi),length(posj));
                F_R(r<r_c) = sigma .* w_R(r<r_c) .* epsilon(r<r_c) .* r_hat(r<r_c);
                F_R(r<r_c) = sigma .* w_R(r<r_c) .* r_hat(r<r_c);
    
                F_all = F_C + F_D + F_R;
                
                F_i = sum(F_all, 2, 'omitnan');
                F_j = sum(F_all, 1, 'omitnan').';
    
                %save new acceleration forces
                accs.(fn{i}) = accs.(fn{i}) + F_i;
                if i~=j
                    accs.(fn{j}) = accs.(fn{j}) - F_j;
                end
            end
            clear posi posj
        end
    end
end


function forces = calculate_molecule_forces(pos, chains, rings, r_S, K_S, L, method)
    n_c = length(chains(:,1));
    n_r = length(rings(:,1));

    forces = struct();
    forces.chainsA = zeros(2*n_c,1);
    forces.chainsB = zeros(5*n_c,1);
    forces.ringsA = zeros(9*n_r,1);

    % IF CHAINS METHOD
    if method == 'c'
        ud = circshift(eye(7),1);
        ud(1,:) = 0;
        ld = circshift(eye(7),-1);
        ld(end,:) = 0;
        dd = ud+ld;
        for c=1:n_c
            posc = [pos.A(chains(c,1));pos.A(chains(c,2));pos.B(chains(c,3));pos.B(chains(c,4));pos.B(chains(c,5));pos.B(chains(c,6));pos.B(chains(c,7))];
            
            dists = posc - posc.';
            dists(real(dists) > L/2) = dists(real(dists) > L/2) - complex(L,0);
            dists(real(dists) < -L/2) = dists(real(dists) < -L/2) + complex(L,0);
            dists(imag(dists) > L/2) = dists(imag(dists) > L/2) - complex(0,L);
            dists(imag(dists) < -L/2) = dists(imag(dists) < -L/2) + complex(0,L);
    
            r = abs(dists);
            r_hat = dists ./ r;
            
            F_S = K_S .* (1 - (r ./ r_S)) .* r_hat;
            
            F_S = F_S .* dd;
            F_S = sum(F_S,1,'omitnan').';
    
            forces.chainsA(((c-1)*2)+1:((c-1)*2)+2,1) = -F_S(1:2);
            forces.chainsB(((c-1)*5)+1:((c-1)*5)+5,1) = -F_S(3:7);
    
        end
    % IF RINGS METHOD
    elseif method == 'r'
        for c=1:n_r
            posr = [pos.A(rings(c,1));pos.A(rings(c,2));pos.A(rings(c,3));pos.A(rings(c,4));pos.A(rings(c,5));pos.A(rings(c,6));pos.A(rings(c,7));pos.A(rings(c,8));pos.A(rings(c,9))];
            
            dists = posr - posr.';
            dists(real(dists) > L/2) = dists(real(dists) > L/2) - complex(L,0);
            dists(real(dists) < -L/2) = dists(real(dists) < -L/2) + complex(L,0);
            dists(imag(dists) > L/2) = dists(imag(dists) > L/2) - complex(0,L);
            dists(imag(dists) < -L/2) = dists(imag(dists) < -L/2) + complex(0,L);
    
            r = abs(dists);
            r_hat = dists ./ r;
            
            F_S = K_S .* (1 - (r ./ r_S)) .* r_hat;
            
            ud = circshift(eye(length(F_S)),1);
            ld = circshift(eye(length(F_S)),-1);
            dd = ud+ld;
            
            F_S = F_S .* dd;
            F_S = sum(F_S,1,'omitnan');
    
            forces.ringsA(((c-1)*9)+1:((c-1)*9)+9,1) = -F_S(1:end).'; 
        end
    else

    end
end


