% Exports the Kobuk Lake geometry in Kobuk.mat (from SubzeroDebugging/bonds_coastline/Initialize_Model)
% to CSV files read by examples/kobuk_lake.jl. Each row is `polygon id, x, y` (m), and each
% polygon's ring is listed once without repeating its first vertex.
%   kobuk_coast.csv: boundaryFloes, the land around the lake, used as topography
%   kobuk_ice.csv:   KobukIce, the region the initial ice field fills
% Usage: export_kobuk(path_to_Kobuk_mat)
function export_kobuk(mat_path)
    load(mat_path, "boundaryFloes", "KobukIce");
    out_dir = fileparts(mfilename("fullpath"));
    write_polys(fullfile(out_dir, "kobuk_coast.csv"), boundaryFloes);
    write_polys(fullfile(out_dir, "kobuk_ice.csv"), KobukIce);
end

function write_polys(fn, polys)
    rows = [];
    id = 0;
    for p = polys
        if area(p) <= 0, continue; end
        for r = regions(p)'
            id = id + 1;
            v = r.Vertices;
            v = v(all(~isnan(v), 2), :);  % holes would be NaN-separated; these polygons have none
            rows = [rows; id * ones(size(v, 1), 1), v];
        end
    end
    writematrix(rows, fn);
    fprintf("%s: %d polygons, %d vertices\n", fn, id, size(rows, 1));
end
