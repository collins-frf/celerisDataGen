%% required input 
% x, y, b, dx_target, water_level_change, bc, jj
% script to plot and write bathy

x=cast(x, 'double');
y=cast(y, 'double');
B=cast(B, 'double');
instruments=cast(instruments, 'double');
xoffset=cast(xoffset, 'double');
yoffset=cast(yoffset, 'double');

xo=x;
yo=y;
Bo=B;

range_x=max(xo)-min(xo);

range_y=max(yo)-min(yo);

dy_target=dx_target;
x=[min(xo):dx_target:max(xo)];
y=[min(yo):dy_target:max(yo)];

nx=length(x);
ny=length(y);
[myxo, myyo] = meshgrid(xo, yo);

[myx, myy] = meshgrid(x,y);

B=interp2(myxo, myyo, Bo', myx, myy);  % interpolate to target
B=B'-water_level_change;  % shift bath/topo for water level datum

% set depth along wave boundary constant value
irrbc_ind_all=find(bc==4);
hbc=[];
for ii=1:length(irrbc_ind_all)   % create spectrum input file if any bc = 4
    irrbc_ind=irrbc_ind_all(ii);
    if irrbc_ind==1
        hbc=[hbc -B(1,:)];
    elseif irrbc_ind==2
        hbc=[hbc -B(nx,:)];
    elseif irrbc_ind==3
        hbc=[hbc -B(:,1)'];
    elseif irrbc_ind==4
        hbc=[hbc -B(:,ny)'];
    end
end

%        hbc=-mean(hbc);
hbc=-8; % force generation area depth to "wavedepth", since this is what is provided by CDIP data points

strip_len=10;                      % number of points to smooth
irrbc_ind_all=find(bc==4);
for ii=1:length(irrbc_ind_all)  %
    irrbc_ind=irrbc_ind_all(ii);
    if irrbc_ind==1
        B(1:strip_len,:)=B(1:strip_len,:)*0+hbc;
        for j=1:ny
            slope=(hbc-B(strip_len*2,j))/strip_len;
            for i=strip_len+1:strip_len*2
                B(i,j)=B(i-1,j)-slope;
            end
        end
        
    elseif irrbc_ind==2
        B(nx-strip_len:nx,:)=B(nx-strip_len:nx,:)*0+hbc;
        for j=1:ny
            slope=(hbc-B(nx-strip_len*2,j))/strip_len;
            for i=nx-strip_len-1:-1:nx-2.*strip_len
                B(i,j)=B(i+1,j)-slope;
            end
        end
        
    elseif irrbc_ind==3
        B(:,1:strip_len)=B(:,1:strip_len)*0+hbc;
        for i=1:nx
            slope=(hbc-B(i,strip_len*2))/strip_len;
            for j=strip_len+1:strip_len*2
                B(i,j)=B(i,j-1)-slope;
            end
        end
        
    elseif irrbc_ind==4
        B(:,ny-strip_len:ny)=B(:,ny-strip_len:ny)*0+hbc;
        for i=1:nx
            slope=(hbc-B(i,ny-strip_len*2))/strip_len;
            for j=ny-strip_len-1:-1:ny-2.*strip_len
                B(i,j)=B(i,j+1)-slope;
            end
        end
        
    end
end


% plot bathy/topo with boundary condition info
hf2=figure(2);
clf
subplot('Position',[0.05 0.05 .9 .9])
pcolor(x,y,B')
hold on
shading interp
xlabel('East-West, x (m)','FontSize',5)
ylabel('North-South, y (m)','FontSize',5)
title(['Bathy/Topo Grid Using Grid Size (m): ' num2str(dx_target)],'FontSize',5)
set(gca,'fontsize',5)
axis([-Inf Inf -Inf Inf])
caxis([min(min(B)) 15])
colorbar
file_name_cbf = 'matlab_launch.cbf';
m_width=abs(x(nx)-x(1));
m_length=abs(y(ny)-y(1));
if nx*ny>8.0e6
    num_nodes= 0.1*(round(nx*ny/1e5));
    disp(['WARNING: Total number of nodes is large (' num2str(num_nodes) 'M), your graphics card may run out of memory'])
end
write_bathy(file_name_cbf,nx,ny,B)

% Display input info on bathy image
for i=1:4
    if i==1
        cboundry='westBoundary';
        xc=x(1)+m_width/50;
        yc=y(round(ny/4));
    elseif i==2
        cboundry='eastBoundary';
        xc=x(nx)-m_width/50;
        yc=y(round(ny/4));
    elseif i==3
        cboundry='southBoundary';
        xc=x(round(nx/4));
        yc=y(1)+m_length/50.;
    elseif i==4
        cboundry='northBoundary';
        xc=x(round(nx/4));
        yc=y(ny)-m_length/50;
    end
    
    if bc(i)==1 %solid wall
        linec=['Solid'];
    elseif bc(i)==2 % sponge
        linec=['Sponge'];
    elseif bc(i)==3 % input wave
        linec=['SineWave'];
    elseif bc(i)==4 % input wave
        linec=['RandomWave'];
    end
    
    % plot boundary info text
    %     txtx=[cboundry ' is ' linec];
    %     htext=text(xc,yc,txtx,'color','k','FontSize',5);
    %     if i<=2
    %         set(htext,'Rotation',90);
    %     end
    
    % plots sponge layers
    if bc(i)==2 % sponge
        if i==1
            xspg=[x(sponge(i)) x(sponge(i))];
            yspg=[y(1) y(ny)];
        elseif i==2
            xspg=[x(nx-sponge(i)) x(nx-sponge(i))];
            yspg=[y(1) y(ny)];
        elseif i==3
            yspg=[y(sponge(i)) y(sponge(i))];
            xspg=[x(1) x(nx)];
        elseif i==4
            yspg=[y(ny-sponge(i)) y(ny-sponge(i))];
            xspg=[x(1) x(nx)];
        end
        
        plot(xspg,yspg,'w--')
    end
end
