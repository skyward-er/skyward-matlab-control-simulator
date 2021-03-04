function [alpha_degree] = controlAlgorithm(sensorData, z, vz, vMod)

sd = sensorData;
sd.Yf.z = z;
sd.Yf.vz = vz;
sd.Yf.vMod = vMod;

serialbridge("Write", structToSingles(sd));
alpha_degree = serialbridge("Read", 1);

end