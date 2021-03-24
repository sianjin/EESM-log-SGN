function W = tgaxMMSEFilter(varargin)
%tgaxMMSEFilter MMSE linear filter weights
%   W = tgaxMMSEFilter(H,P,N0) returns linear MMSE receive filter given the
%   channel matrix and precoding matrix used at the transmitter.
%
%   P is a Nst-by-Nsts-by-Nlinks matrix containing the MMSE coefficients.
%   Nst is the number of subcarriers, Nsts is the number of space-time
%   nstreams, and Nlinks is the number of links.
%
%   H is a Nst-by-Nt-by-Nr-by-Nlinks array containing the channel. Nt is
%   the number of transmit antennas and Nr is the number of receive
%   antennas.
%
%   N0 is the receiver noise power in Watts. It is assumed to be the same
%   for all links.
%
%   W = tgaxMMSEFilter(HE,N0) returns linear MMSE receive filter where HE
%   is a Nst-by-Nsts-by-Nrx-by-Nlinks array containing the channel
%   response with precoding (spatial mapping + cyclic shift) included. 
%
%   % Example: Generate MMSE equalizer weights
%   H = rand(242,2); % Channel of interest
%   P = ones(242,1,2); % Precoding matrix
%   N0 = -120; % dBW
%   W = tgaxMMSEFilter(H,P,db2pow(N0));
%
%   See also calculateSINR.

%   Copyright 2019 The MathWorks, Inc.

%#codegen

if nargin==3
    % tgaxMMSEFilter(H,P,N0)
    % H is a Nst-by-Ntx-by-Nrx-by-Nlinks array containing the channel
    % response
    % P is a Nst-by-Nsts-by-Ntx-by-Nlinks array containing the precoding
    % matrix (spatial mapping + cyclic shift)
    % N0 is the noise variance
    [H,P,N0] = varargin{:};
    [Nst,Ntx,Nrx,Nlinks] = size(H);
    [Nst_w,Nsts,Ntx_w,Nlinks_w] = size(P);
    assert(all([Nst Ntx Nlinks]==[Nst_w Ntx_w Nlinks_w]))
    
    % The combined channel matrix estimated by the receiver includes
    % precoding and the channel response. Create this as it will be used by
    % the linear receiver
    Hp = permute(H,[1 4 2 3]); % Nst-by-Nlinks-by-Ntx-by-Nrx
    Pp = permute(P,[1 4 3 2]); % Nst-by-Nlinks-by-Ntx-by-Nsts
    Hc = coder.nullcopy(complex(zeros(Nst,Nlinks,Nrx,Nsts)));
    for i = 1:Nrx
        for j = 1:Nsts
            Hc(:,:,i,j) = sum(Hp(:,:,:,i).*Pp(:,:,:,j),3);
        end
    end
    Hcp = permute(Hc,[4 3 1 2]); % Nsts-by-Nrx-by-Nst-by-Nlinks
else
    % tgaxMMSEFilter(HE,N0)
    % HE is a Nst-by-Nsts-by-Nrx-by-Nlinks array containing the channel
    % response with precoding (spatial mapping + cyclic shift) included. N0
    % is the noise variance
    [Hc,N0] = varargin{:};
    [Nst,Nsts,Nrx,Nlinks] = size(Hc);
    Hcp = permute(Hc,[2 3 1 4]); % Nsts-by-Nrx-by-Nst-by-Nlinks
end

% Calculate the linear receiver to be used at the receiver
W = coder.nullcopy(complex(zeros(Nst,Nrx,Nsts,Nlinks)));
eyensts = eye(Nsts);
noiseesteye = N0*eyensts;
for i = 1:Nlinks
    for idx = 1:Nst
        % Calculate the linear receiver matrix
        Hsc = Hcp(:,:,idx,i);
        invH = (Hsc*Hsc'+noiseesteye)\eyensts;
        Wsc = Hsc' * invH;  
        W(idx,1:Nrx,1:Nsts,i) = Wsc;
    end
end
end