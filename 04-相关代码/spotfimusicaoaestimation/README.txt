MATLAB code for Angle of arrival estimation algorithm in SpotFi paper [1].

Setup is as described in section 4 of [1].
Instructions:
1. Calibrate the receiver radio chains as described in [2]
2. reshape CSI into 90x1 vector where first 30 elements correspond to subcarriers for first rx antenna, second 30 correspond to CSI from next 30 subcarriers and last 30 elements correspond to CSI from subcarriers of last receiver antenna and run main.m

[1] Kotaru, Manikanta, et al. "Spotfi: Decimeter level localization using wifi." ACM SIGCOMM Computer Communication Review. Vol. 45. No. 4. ACM, 2015.
[2] Xiong, Jie, and Kyle Jamieson. "ArrayTrack: A Fine-Grained Indoor Location System." NSDI. 2013.
