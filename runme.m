steps=[4];

%We are going to need to split the area into several exp files, here is the first one:
domain = 'Domain1';

org=organizer('repository',['./' domain],'prefix',['Model_' domain '_'],'steps',steps); clear steps;

if perform(org,'ProcessData')% {{{

	%Domain Boundary
	disp('Getting domain boundaries');
	offset=20*10^3;
	A=expread(['./' domain '/' domain '.exp']);
	domainxlim = [min(A(1).x)-offset max(A(1).x)+offset];
	domainylim = [min(A(1).y)-offset max(A(1).y)+offset];

	disp(' => Loading ice thickness data');

	%Read everything
	fid  = fopen('data/radar_5tracks.csv');
	data = textscan(fid,'%f%f%d%d%f%f%f','CommentStyle','#','Delimiter',',','Headerlines',1);
	fid  = fclose(fid);
	x = data{1}; y = data{2}; thickness = data{7};

	%Select only the data that are in the domain
	pos = find(x>domainxlim(1) & x<domainxlim(2) & y>domainylim(1) & y<domainylim(2));
	x=x(pos); y=y(pos); thickness=thickness(pos);

	disp('Saving');
	save(['./' domain '/ProcessedTracks'],'x','y','thickness');
end %}}}
if perform(org,'Mesh')% {{{

	hmax=400;

	%disp('Reading tracks');
	%load([domain '/ProcessedTracks.mat']);
	%numpoints=length(x);
	%tracks=double([x y [1:1:numpoints]']);

	%disp('Removing points outside of the domain');
	%A=expread(['./' domain '/' domain '.exp']);
	%flags=ContourTest(A(1).x,A(1).y,tracks(:,1),tracks(:,2));
	%for i=2:length(A),
	%	flags=(flags & ~ContourTest(A(i).x,A(i).y,tracks(:,1),tracks(:,2)));
	%end
	%tracks=tracks(find(flags),:);

	%disp('Coarsening observations for mesh');
	%tic
	%pos=GroupObs(tracks(:,1),tracks(:,2),tracks(:,3),hmin);
	%tracks=tracks(pos,:);
	%toc
	%clear x y thickness;

	disp('Generating first constrained mesh');
	%md=bamg(model,'domain',['./' domain '/' domain '.exp'],'hmax',hmax,'RequiredVertices',tracks,'MaxCornerAngle',0.001,'splitcorners',0);
	md=bamg(model,'domain',['./' domain '/' domain '.exp'],'hmax',hmax);
	md.mesh.epsg=8833;
	md.miscellaneous.name=domain;

	savemodel(org,md);
end %}}}
if perform(org,'Param')% {{{

	md=loadmodel(org,'Mesh');
	md=setmask(md,'','');
	md=setflowequation(md,'SSA','all');

	%read and interpolate velocities
	ncdata = 'data/velxy_NCAA02_120m.nc';
	x  = ncread(ncdata, 'easting');
	y  = ncread(ncdata, 'northing');
	vx = ncread(ncdata, 'vx')';
   vy = ncread(ncdata, 'vy')';
   md.inversion.vx_obs = InterpFromGridToMesh(x, flipud(y), flipud(vx), md.mesh.x, md.mesh.y, 0.);
   md.inversion.vy_obs = InterpFromGridToMesh(x, flipud(y), flipud(vy), md.mesh.x, md.mesh.y, 0.);

   md.inversion.vel_obs  = sqrt(md.inversion.vx_obs.^2+md.inversion.vy_obs.^2);
   md.initialization.vx  = md.inversion.vx_obs;
   md.initialization.vy  = md.inversion.vy_obs;
   md.initialization.vel = md.inversion.vel_obs;

	%surface elevation
	ncdata = 'data/arcticdem_NCAA02_2f100m.nc';
	x   = ncread(ncdata, 'easting');
	y   = ncread(ncdata, 'northing');
	dat = ncread(ncdata, 'elev')';
	md.geometry.surface = InterpFromGridToMesh(x, flipud(y), flipud(dat), md.mesh.x, md.mesh.y, 0.);

	%thickness
	ncdata = 'data/thickness_NCAA02_100m.nc';
	x   = ncread(ncdata, 'easting');
	y   = ncread(ncdata, 'northing');
	dat = ncread(ncdata, 'thick')';
	md.geometry.thickness = InterpFromGridToMesh(x, flipud(y), flipud(dat), md.mesh.x, md.mesh.y, 0.);
	pos = find(md.geometry.thickness<1);
	md.geometry.thickness(pos) = 1;

	md.geometry.base = md.geometry.surface-md.geometry.thickness;
	md.geometry.bed  = md.geometry.base;

	%SMB
	disp('Converting SMB from mm/yr water eq to m/yr ice eq.');
	ncdata = 'data/SMB_NCAA_2f1km.nc';
	x   = ncread(ncdata, 'easting');
	y   = ncread(ncdata, 'northing');
	dat = ncread(ncdata, 'SMB_rec')' * 1e-3*md.materials.rho_freshwater/md.materials.rho_ice;
	md.smb.mass_balance = InterpFromGridToMesh(x, flipud(y), flipud(dat), md.mesh.x, md.mesh.y, 0.);

	%DHDT
	disp('Converting dHdt from mm/yr to m/yr');
	ncdata = 'data/dhdt_NCAA02_2010-2020_100m.nc';
	x   = ncread(ncdata, 'easting');
	y   = ncread(ncdata, 'northing');
	dat = ncread(ncdata, 'dhdt')' / 1000;
	md.balancethickness.thickening_rate = InterpFromGridToMesh(x, flipud(y), flipud(dat), md.mesh.x, md.mesh.y, 0.);

	%Deal with other params
	md.basalforcings.floatingice_melting_rate = zeros(md.mesh.numberofvertices, 1);
	md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices, 1);

	md=SetIceSheetBC(md);

	savemodel(org,md);
end %}}}

if perform(org,'ReconstructThickness')% {{{

	md=loadmodel(org,'Param');

	%Only run mass transport!
	md.transient.isstressbalance = 0;
	md.transient.isthermal = 0;

	%run one step only
	md.timestepping.final_time=md.timestepping.time_step;

	Hinit = md.geometry.thickness;
	for i=1:20
		md=solve(md, 'tr');

		dHdt = (md.results.TransientSolution(end).Thickness - md.geometry.thickness)/(md.timestepping.final_time - md.timestepping.start_time);
		md.geometry.thickness = md.geometry.thickness + dHdt*1;
		disp(['Mean correction: ' num2str(mean(abs(dHdt)))]);

		pos = find(md.geometry.thickness<1);
		md.geometry.thickness(pos) = 1;
		md.geometry.base = md.geometry.surface-md.geometry.thickness;
		md.geometry.bed  = md.geometry.base;
	end

	%savemodel(org,md);
end %}}}
