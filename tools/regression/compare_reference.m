function ok = compare_reference(fileA, fileB, varargin)
%COMPARE_REFERENCE Compare two HARP-I reference snapshots field by field.
%
%   ok = compare_reference(refA, refB)
%   ok = compare_reference(refA, refB, 'Tol', 1e-12, 'Verbose', true)
%
%   Prints max|d|, rms(d) and the relative error for every shared field, and
%   for the unwrapped-phase fields additionally checks that the integer wrap
%   counts round((a-b)/2pi) are identically zero -- a phase that differs by a
%   whole number of cycles is a different answer, not a rounding difference.
%
%   Returns false (and, under -batch, exits non-zero) if any field exceeds
%   the tolerance, so it can gate a commit.
%
%   See also CAPTURE_REFERENCE.

    p = struct('Tol', 1e-12, 'Verbose', true, 'Exit', false);
    for k = 1:2:numel(varargin)
        p.(validatestring(varargin{k}, fieldnames(p))) = varargin{k+1};
    end

    A = load(fileA);
    B = load(fileB);

    names = intersect(fieldnames(A), fieldnames(B));
    names = setdiff(names, {'meta'});

    onlyA = setdiff(fieldnames(A), fieldnames(B));
    onlyB = setdiff(fieldnames(B), fieldnames(A));
    if ~isempty(onlyA); fprintf(2, 'only in A: %s\n', strjoin(onlyA', ' ')); end
    if ~isempty(onlyB); fprintf(2, 'only in B: %s\n', strjoin(onlyB', ' ')); end

    if p.Verbose
        fprintf('%-16s %12s %12s %12s   %s\n', ...
                'field', 'max|d|', 'rms(d)', 'rel', 'status');
        fprintf('%s\n', repmat('-', 1, 70));
    end

    ok = true;
    for k = 1:numel(names)
        n = names{k};
        if isstruct(A.(n))
            [ok_f, rows] = compare_struct(A.(n), B.(n), n, p.Tol);
        else
            [ok_f, rows] = compare_array(A.(n), B.(n), n, p.Tol);
        end
        ok = ok && ok_f;
        if p.Verbose
            for r = 1:numel(rows); fprintf('%s\n', rows{r}); end
        end
    end

    % integer wrap-count check on the unwrapped phases
    for n = {'uXpha', 'uYpha', 'unwrapped'}
        if isfield(A, n{1}) && isfield(B, n{1})
            d = A.(n{1}) - B.(n{1});
            d = d(isfinite(d));
            if ~isempty(d)
                nwrap = max(abs(round(d/(2*pi))));
                if nwrap > 0
                    fprintf(2, '  WRAP MISMATCH %s: up to %d cycles\n', n{1}, nwrap);
                    ok = false;
                else
                    fprintf('  wrap counts %s: identical\n', n{1});
                end
            end
        end
    end

    if isfield(A,'meta') && isfield(B,'meta')
        fprintf('\nA: sha %s  threads %d  %s\n', A.meta.git_sha, A.meta.threads, A.meta.timestamp);
        fprintf('B: sha %s  threads %d  %s\n', B.meta.git_sha, B.meta.threads, B.meta.timestamp);
        if A.meta.threads ~= B.meta.threads
            fprintf(2, 'WARNING: thread counts differ; results are not bit-reproducible.\n');
        end
    end

    if ok
        fprintf('\nRESULT: identical within tol %g\n', p.Tol);
    else
        fprintf(2, '\nRESULT: DIFFERENCES EXCEED TOLERANCE\n');
    end

    if p.Exit && ~ok; exit(1); end
end

% =========================================================================
function [ok, rows] = compare_struct(a, b, prefix, tol)
    ok = true; rows = {};
    fn = intersect(fieldnames(a), fieldnames(b));
    for k = 1:numel(fn)
        va = a.(fn{k}); vb = b.(fn{k});
        if isnumeric(va) || islogical(va)
            [o, r] = compare_array(va, vb, [prefix '.' fn{k}], tol);
            ok = ok && o; rows = [rows, r]; %#ok<AGROW>
        end
    end
end

% =========================================================================
function [ok, rows] = compare_array(a, b, name, tol)
    rows = {};
    if ~isequal(size(a), size(b))
        rows{end+1} = sprintf('%-16s %s vs %s   SIZE MISMATCH', name, ...
                              mat2str(size(a)), mat2str(size(b)));
        ok = false; return
    end

    a = double(a(:)); b = double(b(:));
    finite = isfinite(a) & isfinite(b);

    % NaN patterns must agree -- a NaN moving is a real change
    if ~isequal(isfinite(a), isfinite(b))
        rows{end+1} = sprintf('%-16s %12s %12s %12s   NaN PATTERN DIFFERS', ...
                              name, '-', '-', '-');
        ok = false; return
    end

    d = a(finite) - b(finite);
    if isempty(d)
        rows{end+1} = sprintf('%-16s %12s %12s %12s   empty', name, '-', '-', '-');
        ok = true; return
    end

    mx  = max(abs(d));
    rms = sqrt(mean(d.^2));
    scale = max(abs(a(finite)));
    rel = mx / max(scale, eps);

    if mx == 0
        status = 'bitwise';
    elseif mx <= tol
        status = 'ok';
    else
        status = 'DIFFERS';
    end
    ok = mx <= tol;

    rows{end+1} = sprintf('%-16s %12.4g %12.4g %12.4g   %s', name, mx, rms, rel, status);
end
