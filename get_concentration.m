function concentration = get_concentration(pollution_source, sensor_location, sensor_time, Para)
% pollution_source:slt

U = Para.U;
D = Para.D;
A = Para.A;
k = Para.k; 

C1 = (sensor_location-pollution_source(2)-U*(sensor_time - pollution_source(3))).^2;
C2 = 4*D*(sensor_time - pollution_source(3));
C3 = -k*(sensor_time - pollution_source(3));
concentration =  pollution_source(1)./(A*sqrt(4*pi*D*(sensor_time - pollution_source(3)))+ 1e-18).*exp(-C1./C2+C3);

concentration(sensor_time < pollution_source(3)) = 0; 
concentration(sensor_location < pollution_source(2)) = 0; 
