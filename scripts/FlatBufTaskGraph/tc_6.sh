#!/bin/bash

tc qdisc del dev pp1eth parent ffff: >/dev/null 2>&1
tc qdisc add dev pp1eth ingress
tc qdisc del dev pp2eth parent ffff: >/dev/null 2>&1
tc qdisc add dev pp2eth ingress
tc qdisc del dev pp3eth parent ffff: >/dev/null 2>&1
tc qdisc add dev pp3eth ingress
tc qdisc del dev pp4eth parent ffff: >/dev/null 2>&1
tc qdisc add dev pp4eth ingress
tc filter add dev pp2eth prio 0 protocol ip parent ffff: flower skip_hw  dst_ip 10.100.2.0/24 action pedit ex munge eth dst set 14:02:ec:ca:a8:ce pipe action mirred egress redirect dev pp3rdma
tc filter add dev pp1eth prio 0 protocol 802.1Q parent ffff: flower skip_hw vlan_ethtype ip dst_ip 10.100.4.0/24 action pedit ex munge eth dst set 34:80:0d:bc:3d:1d pipe action mirred egress redirect dev pp2rdma
tc filter add dev pp2eth prio 0 protocol ip parent ffff: flower skip_hw  dst_ip 10.100.1.0/24 action pedit ex munge eth dst set 14:02:ec:ca:e8:8a pipe action mirred egress redirect dev pp1rdma
tc filter add dev pp1eth prio 0 protocol 802.1Q parent ffff: flower skip_hw vlan_ethtype ip dst_ip 10.100.5.0/24 action pedit ex munge eth dst set 34:80:0d:bc:3d:19 pipe action mirred egress redirect dev pp2rdma
tc filter add dev pp3eth prio 0 protocol 802.1Q parent ffff: flower skip_hw vlan_ethtype ip dst_ip 10.100.5.0/24 action pedit ex munge eth dst set 34:80:0d:bc:3d:19 pipe action mirred egress redirect dev pp2rdma
tc filter add dev pp2eth prio 0 protocol 802.1Q parent ffff: flower skip_hw vlan_ethtype ip dst_ip 10.100.2.0/24 action pedit ex munge eth dst set 14:02:ec:ca:a8:ce pipe action mirred egress redirect dev pp3rdma
tc filter add dev pp2eth prio 0 protocol 802.1Q parent ffff: flower skip_hw vlan_ethtype ip dst_ip 10.100.7.0/24 action pedit ex munge eth dst set 14:02:ec:ca:a8:ca pipe action mirred egress redirect dev pp3rdma
tc filter add dev pp1eth prio 0 protocol ip parent ffff: flower skip_hw  dst_ip 10.100.4.0/24 action pedit ex munge eth dst set 34:80:0d:bc:3d:1d pipe action mirred egress redirect dev pp2rdma
tc filter add dev pp2eth prio 0 protocol ip parent ffff: flower skip_hw  dst_ip 10.100.7.0/24 action pedit ex munge eth dst set 14:02:ec:ca:a8:ca pipe action mirred egress redirect dev pp3rdma
tc filter add dev pp2eth prio 0 protocol 802.1Q parent ffff: flower skip_hw vlan_ethtype ip dst_ip 10.100.1.0/24 action pedit ex munge eth dst set 14:02:ec:ca:e8:8a pipe action mirred egress redirect dev pp1rdma
tc filter add dev pp3eth prio 0 protocol ip parent ffff: flower skip_hw  dst_ip 10.100.5.0/24 action pedit ex munge eth dst set 34:80:0d:bc:3d:19 pipe action mirred egress redirect dev pp2rdma
tc filter add dev pp1eth prio 0 protocol ip parent ffff: flower skip_hw  dst_ip 10.100.5.0/24 action pedit ex munge eth dst set 34:80:0d:bc:3d:19 pipe action mirred egress redirect dev pp2rdma
arp -d 10.100.0.4 > /dev/null 2>&1
ip r d 10.100.0.4/32 > /dev/null 2>&1
arp -s 10.100.0.4 14:02:ec:ca:a8:ce -i pp3rdma
ip r a 10.100.0.4/32 src 10.100.6.3 dev pp3rdma
arp -d 10.100.1.3 > /dev/null 2>&1
ip r d 10.100.1.3/32 > /dev/null 2>&1
arp -s 10.100.1.3 14:02:ec:ca:e8:8a -i pp1rdma
ip r a 10.100.1.3/32 src 10.100.6.1 dev pp1rdma
arp -d 10.100.2.3 > /dev/null 2>&1
ip r d 10.100.2.3/32 > /dev/null 2>&1
arp -s 10.100.2.3 14:02:ec:ca:a8:ce -i pp3rdma
ip r a 10.100.2.3/32 src 10.100.6.3 dev pp3rdma
arp -d 10.100.3.3 > /dev/null 2>&1
ip r d 10.100.3.3/32 > /dev/null 2>&1
arp -s 10.100.3.3 14:02:ec:ca:a8:ce -i pp3rdma
ip r a 10.100.3.3/32 src 10.100.6.3 dev pp3rdma
arp -d 10.100.4.2 > /dev/null 2>&1
ip r d 10.100.4.2/32 > /dev/null 2>&1
arp -s 10.100.4.2 34:80:0d:bc:3d:1d -i pp2rdma
ip r a 10.100.4.2/32 src 10.100.6.2 dev pp2rdma
arp -d 10.100.5.3 > /dev/null 2>&1
ip r d 10.100.5.3/32 > /dev/null 2>&1
arp -s 10.100.5.3 34:80:0d:bc:3d:19 -i pp2rdma
ip r a 10.100.5.3/32 src 10.100.6.2 dev pp2rdma
arp -d 10.100.7.3 > /dev/null 2>&1
ip r d 10.100.7.3/32 > /dev/null 2>&1
arp -s 10.100.7.3 14:02:ec:ca:a8:ca -i pp3rdma
ip r a 10.100.7.3/32 src 10.100.6.3 dev pp3rdma
arp -d 10.100.8.3 > /dev/null 2>&1
ip r d 10.100.8.3/32 > /dev/null 2>&1
arp -s 10.100.8.3 14:02:ec:ca:a8:ce -i pp3rdma
ip r a 10.100.8.3/32 src 10.100.6.3 dev pp3rdma
arp -d 10.100.9.4 > /dev/null 2>&1
ip r d 10.100.9.4/32 > /dev/null 2>&1
arp -s 10.100.9.4 34:80:0d:bc:3d:c7 -i pp4rdma
ip r a 10.100.9.4/32 src 10.100.6.4 dev pp4rdma
arp -d 10.100.10.4 > /dev/null 2>&1
ip r d 10.100.10.4/32 > /dev/null 2>&1
arp -s 10.100.10.4 34:80:0d:bc:3d:c7 -i pp4rdma
ip r a 10.100.10.4/32 src 10.100.6.4 dev pp4rdma
arp -d 10.100.11.3 > /dev/null 2>&1
ip r d 10.100.11.3/32 > /dev/null 2>&1
arp -s 10.100.11.3 34:80:0d:bc:3d:c3 -i pp4rdma
ip r a 10.100.11.3/32 src 10.100.6.4 dev pp4rdma
