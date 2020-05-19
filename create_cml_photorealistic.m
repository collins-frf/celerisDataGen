function create_cml(file_name_cml,timestep,Mannings_n,m_width,m_length,nx,ny,wave_type,H,T,theta,xc,yc,bc,sponge,min_depth,brk1,brk2,inst_ind,range_x,range_y,Courant_c)


fileID = fopen(file_name_cml,'w');

line1='<?xml version="1.0" ?>\n';

line2='<Experiment>\n';
line3='    <name>USC Coastal Waves Forecast</name>\n';
line4='    <!-- Settings for Model -->\n';
line5='    <model type = "BSNQ">\n';
line6=['        <parameters epsilon = ' num2str(min_depth) ' Theta = 2.0 correctionStepsNum = 0 timestep = ' num2str(timestep) '  adaptive = true></parameters>\n'];
line7=['        <friction type = "Quadratic" coef = ' num2str(Mannings_n) '></friction>\n'];
line8=['    </model>\n'];

line9=['\n'];

line10=['	<!-- Settings for Solution field -->\n'];
line11=['	<fieldDimensions width = ' num2str(m_width) ' length = ' num2str(m_length) ' stillWaterElevation = 0></fieldDimensions>\n'];
line12=['	<gridSize nx = ' num2str(nx) ' ny = ' num2str(ny) '></gridSize>\n'];
cdir=['matlab_launch.cbf'];
line13=['	<bathymetryFilePath> ' cdir ' </bathymetryFilePath>\n'];
   
line14=['   <tideSurgeSLR auto = false min = -10 max = 10 set = 0.0></tideSurgeSLR> \n'];

line15=['	<!-- Settings for Initial Condition -->\n'];
if wave_type==1
    line16=['	<solitaryWave H = ' num2str(H) ' theta = ' num2str(theta) ' xc = ' num2str(xc) ' yc = ' num2str(yc) '></solitaryWave>\n'];
else
    line16=['\n'];
end

line17=['\n'];

line18=['	<!-- Settings for Boundaries-->\n'];
for i=1:4
    if i==1
        cboundry='westBoundary';
        bc_c=bc(1);
    elseif i==2
        cboundry='eastBoundary';
        bc_c=bc(2);
    elseif i==3
        cboundry='southBoundary';
        bc_c=bc(3);
    elseif i==4
        cboundry='northBoundary';
        bc_c=bc(4);
    end
    
    if bc_c==1 %solid wall
        linec=['	<' cboundry ' type = "Solid"  seaLevel = 0 widthNum = 2></' cboundry '>\n'];
    elseif bc_c==2 % sponge
        linec=['	<' cboundry ' type = "Sponge" seaLevel = 0 widthNum = ' num2str(sponge(i)) '></' cboundry '>\n'];
    elseif bc_c==3 % single harmonic wave
        linec=['	<' cboundry ' type = "SineWave"  seaLevel = 0 widthNum = 2><sineWave amplitude = ' num2str(H/2) ' period = ' num2str(T) ' theta = ' num2str(theta) '></sineWave></' cboundry '>\n'];
    elseif bc_c==4 % irregular wave
        linec=['	<' cboundry ' type = "IrregularWaves" seaLevel = 0 widthNum = 5><filePath> irrWaves.txt </filePath></' cboundry '>\n'];
    elseif bc_c==5 % time series input
        linec=['	<' cboundry ' type = "UniformTimeSeries" seaLevel = 0 widthNum =><filePath> InputTS.txt </filePath></' cboundry '>\n'];
    end
    
    eval(['line' num2str(18+i) '=linec;'])
end


line23=['\n'];

line24=['	<!-- Settings for Logging Data-->\n'];

logStep=20; % provides enough points
line25=['	<logData doLog = true logStep = ' num2str(logStep) '>\n'];
line26=['		<logPath> </logPath>\n'];

line27=['		<range filename = "array"><bottomLeft x = ' num2str(range_x(1)) ' y = ' num2str(range_y) '></bottomLeft><topRight x = ' num2str(range_x(2)) ' y = ' num2str(range_y) '></topRight></range>\n'];

%gage_str=[];
%for i=1:length(inst_ind(:,1))
%    gage_str=[gage_str num2str(inst_ind(i,1)) ',' num2str(inst_ind(i,2)) ','];
%end
%gage_str=gage_str(1:length(gage_str)-1);
%line28=['		<gauges filename = "gauges">' gage_str '</gauges>\n'];

line28=[' \n'];

line29=['	</logData>\n'];

line30=['\n'];

line31=['	<!-- Settings for graphics are optional-->\n'];

line32=['	<graphics>\n'];
line33=['			<vertical scale = 1></vertical>\n'];
line34=['			<!-- Photorealistic = 0, PARULA = 1, JET = 2 -->\n'];
line35=['			<surfaceShading type = 0>\n'];
			
line36=['				<!-- Eta = 0, U = 1, V = 2, |U+V| = 3, Vorticity = 4 -->\n'];
line37=['				<shadingVariable value =  0></shadingVariable>\n'];

line38=['				<!-- -minMaxValue < value < +minMaxValue -->\n'];
line39=['				<colormap auto = false min = ' num2str(-1.0*H) ' max = ' num2str(2.0*H) '></colormap>\n'];
line40=['			    <drylandDepthOfInundation auto = true show = true value = 0.1 max = 1></drylandDepthOfInundation>\n'];
line41=['			    <dissipationIntensity show = true threshold = ' num2str(brk1) ' decay = ' num2str(brk2) '></dissipationIntensity>\n'];
line42=['			</surfaceShading>\n'];
		
line43=['			<!-- Textures = 0 to 5, COLORMAP = 6, CONTOURS = 7 -->\n'];
line44=['			<terrainTexture type = 6><colormap auto = true min = 0 max = 3></colormap>\n'];
line45=['			</terrainTexture>\n'];
		
line46=['			<lighting ambient = 1.1  sunAltitude = 70. sunAzimuth = 70></lighting>\n'];
line47=['			<fresnel coef = 0.23 refractive_index = 3 attenuation_1 = 0 attenuation_2 = 0.29></fresnel>\n'];
line48=['			<camera auto = false FOV = 70 x = 0 y = 300 z = ' num2str(sqrt(m_width^2+m_length^2)/1.75) ' pitch = -20 yaw = 0></camera>\n'];
line49=['			<window maximized = true></window>\n'];
line50=['			<skybox type = 2></skybox>\n'];
line51=['			<grid show = false scale = 1></grid>\n'];

line52=['		</graphics>\n'];

line53=['</Experiment>\n'];


numlines=53;
for i=1:numlines
    eval(['cline = line' num2str(i) ';'])
    fprintf(fileID,cline);	
end
    


fclose(fileID);
