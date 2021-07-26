#!/usr/bin/expect -f
set timeout -1
spawn cu -l /dev/ttyACM0 -s 115200
expect "DONE"
send "~.\r"
