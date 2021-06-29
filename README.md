# EESM-log-SGN
## Overview
For achieving higher network spectral efficiencies, the main wireless network families (3GPP mobile cellular networks and IEEE 802.11 WLANs) are implementing small cells (higher node density) and more advanced physical (PHY) layer technologies. For example, IEEE 802.11ax dense deployment scenarios include enterprise, multi-dwelling establishments, and crowded public hotspots (malls, airports, sports stadiums, etc.); the number of antennas per device scales to 8 for Multi-Input Multi-Output (MIMO) operation and channel bandwidths scale up to 160 MHz. In addition, new modes of Multi-User (MU) transmission are introduced, based on Orthogonal Frequency Division Multiple Access (OFDMA) and Multi-User MIMO (MU-MIMO). Performance evaluation of such (even modest scale) networks in various operational scenarios of interest via pure mathematical analyses is infeasible and field trials are costly. This leaves system-level (network) simulations as the most feasible option, rendering network simulators like ns-3 indispensable for investigating wireless network performance as a function of network dimensionality.

With increasing node density and complexity in the PHY layer, associated computational complexity is a major concern for network simulators. System-level simulations quantify the network performance (e.g., throughput) from the packet-level performance metrics such as instantaneous (per-packet) Packet Error Ratio (PER). In contrast, full PHY (link) simulations, typically run on link simulators, obtain packet performance at a single receiver, using runtime-costly symbol-level simulation of the PHY layer. Running a full PHY simulation that involves generating channel realizations and transceiver signal processing is impractical within a network simulator. For example, the average runtimes of a full PHY simulation for 40000 PHY packets on a single Orthogonal Frequency Division Multiplexing (OFDM) MIMO link require order of hour runtimes. Thus, system-level simulation must incorporate suitable PHY layer abstractions that represent the PHY layer performance with sufficient accuracy to be incorporated within network simulators at runtime. The PHY layer abstraction in a network simulator is the packet error model that produces a decision on whether the packet is successfully received or not based on the PHY layer setup including the received (RX) signal/interference/noise power, the modulation and coding type, the channel models, etc.

![alt text](https://depts.washington.edu/funlab/wp-content/uploads/2021/04/traditionalPhyAbs.png)
Figure 1: Flow chart of traditional PHY layer abstraction

Traditional PHY layer abstraction suggested by IEEE TGax group is shown in Figure 1. This PHY layer abstraction compress channel, OFDM MIMO setups into a single metric called effective SINR, which is then mapped into instantaneous PER for system simulations. Our contribution finds while these earlier PHY layer abstractions do reduce runtimes compared with full PHY simulation, they still suffer from scaling with MIMO dimensions (the number of transmit antennas and the number of receive antennas), MU dimensions (the number of users that are simultaneously served), bandwidth, and the number of interferers. The reason is that it requires generating channel, precoding and the decoding matrices, and calculating post-MIMO processing SINR matrices online. These operations involve expensive matrix calculations that scale with MIMO dimensions, MU dimensions, bandwidth, and the number of interferes. This motivates us to create a new PHY layer abstraction whose runtime is insensitive to the above dimensionality variations.

![alt text](https://depts.washington.edu/funlab/wp-content/uploads/2021/04/validation2.png)

Figure 2: Approximating effective SINR by log-SGN distribution under OFDMA allocation with 52 subcarriers, 8 × 2 MIMO with 2 streams, TGax channel model-D, MCS4. For the interference case, the RX INR is 20dB lower than the RX SNR

Our approach is driven by a key underlying question: since the link performance (instantaneous PER) only needs the effective SINR, can we bypass channel generation, precoder calculation, decoder calculation as well as post-MIMO processing SINR matrix calculation steps, and directly model effective SINR? Under frequency-selective independent and identically distributed (i.i.d.) block fading channels, we discover that the effective SINR distribution under Exponential Effective SINR Mapping (EESM) L2S mapping can be well approximated by a 4-parameter distribution called log-SGN distribution. This result is shown in Figure 2. Therefore, under EESM L2S mapping, we only need to store 4 log-SGN parameters under a specific PHY layer setup (e.g., channel type, MIMO dimension, MCS, etc.) and draw effective SINR samples from this specific log-SGN distribution directly. We call such an approach the EESM-log-SGN PHY layer abstraction, whose flow chart is shown in Figure 2. The payoff is significant: the runtime of the traditional PHY layer abstraction in Figure 1 scales as system dimensionality increases, while the runtime of the EESM-log-SGN PHY layer abstraction in Figure 2 is much reduced and is insensitive to the dimensionality change. We show the runtime results in Table 1.

![alt text](https://depts.washington.edu/funlab/wp-content/uploads/2021/04/EESMlogSGNAbs.png)
Figure 3: Flow chart of the proposed EESM-log-SGN PHY layer abstraction

![alt text](https://depts.washington.edu/funlab/wp-content/uploads/2021/04/runtimeComp2.png)

Table 1: Average runtime comparison between the traditional EESM PHY layer abstraction and the proposed EESM-log-SGN PHY layer abstraction for running a 40000-packet simulation at a specific RX SNR in MATLAB

For more detailed principles of EESM-log-SGN, please see the IEEE TCOM paper: “Efficient PHY Layer Abstraction for Fast Simulations in Complex System Environments”. 

## About this code
This code provides the MATLAB implementation of EESM-log-SGN. The prerequisites are MATLAB 2020b (or later version) and MATLAB WLAN toolbox in MATLAB 2020b (or later version). The document named "EESMlogSGN_quickStart" shows a quick start guide of our code. Multiple step-by-step MATLAB examples are also provided as tutorials for implementing of EESM-log-SGN in MATLAB.

## Code tutorial video and slides
Video: https://vimeo.com/567675867

Slides: https://www.nsnam.org/wp-content/uploads/2021/tutorials/EESM-log-SGN-tutorial.pptx

