function assertNoUnknownArgs(fname, other_args)
%ASSERTNOUNKNOWNARGS Warn about parameters a function did not recognise.
%
%   assertNoUnknownArgs(mfilename, other_args)
%
%   PARSEINPUTS returns unrecognised parameter/value pairs in its second
%   output, which every caller in this toolbox used to discard. A misspelled
%   or renamed parameter therefore did nothing at all and nobody was told --
%   the cause of several long-lived bugs, e.g. getSegmentation passing
%   'Contours'/'Resolution' to a getMask that expected 'Contour'/'ContourRes'
%   (the mask came back empty every time), and 'Frames' being silently
%   ignored by getStrainBySegments.
%
%   Emits a warning by default. Set the global HARPI_STRICT_ARGS to true to
%   turn it into an error, which is what the test suite does:
%
%       global HARPI_STRICT_ARGS; HARPI_STRICT_ARGS = true;
%
%   Note this is deliberately NOT built into parseinputs: that function
%   shadows an identical copy in DENSEanalysis, whose callers (unwrap2,
%   imcline, ...) pass extra arguments on purpose.
%
%   See also PARSEINPUTS.

    if isempty(other_args)
        return
    end

    if isstruct(other_args)
        names = fieldnames(other_args);
    elseif iscell(other_args)
        names = other_args(1:2:end);
    else
        names = {};
    end
    if isempty(names)
        return
    end

    names = names(cellfun(@ischar, names) | cellfun(@isstring, names));
    if isempty(names)
        return
    end

    msg = sprintf('%s: unrecognised parameter(s) ignored: %s', ...
                  fname, strjoin(cellfun(@char, names(:)', 'uni', 0), ', '));

    global HARPI_STRICT_ARGS %#ok<GVMIS>
    if ~isempty(HARPI_STRICT_ARGS) && HARPI_STRICT_ARGS
        error([fname ':unknownParameter'], '%s', msg);
    else
        warning([fname ':unknownParameter'], '%s', msg);
    end
end
