% GEN_REFS  Generate the MATLAB reference the Python tests compare against.
%
% Run from anywhere:
%
%   matlab -batch "run('<repo>/python/tests/matlab/gen_refs.m')"
%
% Single-threaded, because MATLAB is only bit-reproducible at a fixed thread
% count. Writes python/tests/data/ref_2D.mat (v7.3, readable by h5py).
%
% The output is the same stage-by-stage snapshot the MATLAB regression
% harness uses, so both suites are checking against identical arrays.

maxNumCompThreads(1);

here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));   % .../HARP-I

addpath(fullfile(repo, 'tools'));
harpi_paths();
addpath(fullfile(repo, 'tools', 'regression'));

global HARPI_STRICT_ARGS; HARPI_STRICT_ARGS = true; %#ok<GVMIS>

outdir = fullfile(repo, 'python', 'tests', 'data');
if ~isfolder(outdir); mkdir(outdir); end

capture_reference('Output', fullfile(outdir, 'ref_2D.mat'));

fprintf('\nReference written. Run the Python tests with:\n');
fprintf('  cd %s/python && PYTHONPATH=src python3 -m pytest tests -q\n', repo);
