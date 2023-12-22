function [diams, density] = simulatePSDs(parameters, settings)
%simulatePSDs  Simulates particle size distributions.
%
%   simulatePSDs(parameters, settings) simulates particle size
%   distributions based on model parameters and other settings.
%       *   parameters is FIXME
%       *   settings is FIXME



    % Name of the executable
    exec_name = settings.executable;

    % Input file with parameters that c++ program will read
    inp_file = settings.inp_file;
    % Write parameters to this file
    fileID = fopen(inp_file,'w');
    fprintf(fileID,'%.40f\n',parameters);
    fclose(fileID);

    % Output file with diameters
    diam_out = settings.diameter_file;

    % Output file with densities for each time
    dens_out = settings.density_file;

    % Put the command that needs to be run together
    command = exec_name + " " + inp_file + " " + diam_out + " " + dens_out;

    % Call c++ program
    [status, cmdout] = system(command);
    if status ~= 0
        msg = "Error: Running the command `" + command + "' resulting in error code " + num2str(status);
        error(msg)
    end
    
    if settings.verbose
        fprintf(cmdout)
    end

    % Read results to matlab variables
    diams = readmatrix(diam_out);
    density = readmatrix(dens_out);
    density = density';
    area = trapz(diams, density);
    density = density./area;
end