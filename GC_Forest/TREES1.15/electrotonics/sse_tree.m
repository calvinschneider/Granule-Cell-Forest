% SSE_TREE   steady-state electrotonic signature of a tree.
% (trees package)
%
% sse = sse_tree (intree, I, options)
% -----------------------------------
%
% calculates the steady state matrix describing the electrotonic properties
% of the neuron in the trees format. Each column i is the potential
% distribution during injection of current into compartment i and the
% diagonal is therefore the local input resistances in each compartment.
% sse is then symmetric.
% If input current I is not identity matrix then H columns in sse
% correspond to potential distributions in separate experiments
% corresponding to the input current distribution in that column. Note that
% sse is obtained by inverse matrix calculation and therefore goes very
% quickly but takes memory. In special cases it is advisable to split calls
% in several I input matrices
% 
% Input
% -----
% - intree::integer:index of tree in trees or structured tree
% - I::NxH matrix or value:(optional) current injection vector
%     if I is a number, then 1 nA is injected in position I)
%     if I is omitted I is the identity matrix {DEFAULT}
% - options::string: {DEFAULT: ''}
%     '-s' : show - full matrix if I is left empty (full sse)
%                 - tree distribution if I is Nx1 vector
%                 - other Is first column
%
% Output
% ------
% - sse::NxH matrix: electrotonic signature matrix
%
% Examples
% --------
% sse_tree (sample_tree, [], '-s')
% sse_tree (sample_tree, 11, '-s')
%
% See also M_tree
% Uses M_tree ver_tree
%
% the TREES toolbox: edit, visualize and analyze neuronal trees
% Copyright (C) 2009  Hermann Cuntz

function sse = sse_tree (intree, I, options)

% trees : contains the tree structures in the trees package
global trees

if (nargin < 1)||isempty(intree),
    intree = length (trees);
end;

ver_tree (intree);

if (nargin < 3)||isempty(options),
    options = '';
end

M = M_tree (intree);

if (nargin<2)||isempty(I),
    sse = full (inv (M));
else
    if numel(I)==1,
        dI = I;
        I  = sparse (size (M, 1), 1); I (dI) = 1;
    end
    sse = full (M \ I);
end

if strfind (options, '-s'),
    if numel (M) == numel (sse)
        clf; imagesc (sse); colorbar; axis image;
        xlabel ('node #'); ylabel ('node #');
        title  ('potential distribution [mV]');
    else
        clf; shine; hold on; plot_tree (intree, sse (:, 1)); colorbar;
        title  ('potential distribution [mV]');
        xlabel ('x [\mum]'); ylabel ('y [\mum]'); zlabel ('z [\mum]');
        view (2); grid on; axis image;
    end
end
