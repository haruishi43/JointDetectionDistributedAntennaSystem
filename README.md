# Antenna and User Scheduling Simulation using MATLAB

This repo is for experiments used in my paper ["System Throughput Analysis on Usin Joint Detection in Disributed Antenna System"](https://ieeexplore.ieee.org/document/8690889) for VTC2018 and ["Effect of Joint Detection on System Throughput in Distributed Antenna Network"](https://www.jstage.jst.go.jp/article/transcom/E102.B/3/E102.B_2018EBP3077/_article) for IEICE Journal.
If this repo was helpful for your research, please cite this work as:
```
@article{JointDistributedAntennaNetwork,
  title={System Throughput Analysis on Usin Joint Detection in Disributed Antenna System},
  author={Ishikawa, Haruya and Sanada, Yukitoshi},
  journal={IEEE 88th Vehicular Technology Conference (VTC-Fall)},
  year={2018}
}
@article{JointDistributedAntennaNetworkIEICE,
  title={Effect of Joint Detection on System Throughput in Distributed Antenna Network},
  author={Ishikawa, Haruya and Sanada, Yukitoshi},
  journal={IEICE},
  year={2019}
}
```

## About this project:

Throughput analysis of the usage of joint maximum likelihood detection (MLD) in a distributed antenna system (DAS) with multiple mobile terminal (MT) scheduling is evaluated in this paper. In the DAS, mobile terminals are closer to desired antennas, which improves system throughput and increases frequency utilization efficiency. However, since there are more antennas in a cell, interference can occur between transmitted signals from the antennas when multiple MTs are scheduled. Therefore, we propose a novel MT scheduling method and the use of joint-MLD to reduce the effects of interference. A system level simulation shows that the usage of joint-MLD in densely packed DAS provides better system throughput. 

### Distributed Antenna System (DAS)

- One of the solution being proposed in 5G to accommodate for the explosive traffic growth in recent years.
- In a cell, there exists multiple DAs which are controlled by a central BBU.
- It reduces power consumption and increases system capacity (mitigates channel impacts such as shadowing and fading).
- Channel capacity could be increased with frequency reuse, but the co-channel interference increases as well.

### Joint Maximum-Likelihood Detection (Joint-MLD)

![img1]("https://github.com/Toraudonn/Switching_BS_Matlab/edit/master/Graphs/Screen\ Shot\ 2020-08-14\ at\ 11.18.14.png")

Joint MLD treats the received signals as the superposition of the desired and the interference signals and regards them as a signal with a larger number of constellation points.

### Cell Model and Scheduling

![img2]("./Graphs/Screen\ Shot\ 2020-08-14\ at\ 11.20.34.png")


### Results

![img3]("./Graphs/Screen\ Shot\ 2020-08-14\ at\ 11.22.11.png")

