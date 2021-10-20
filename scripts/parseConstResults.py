import sys

with open("kem") as f:
    lines = f.read().splitlines()
out = sys.argv[1]
out = str(out) + ", " + str(lines[4]) + ", " + str(lines[6]) + ", " + str(lines[8])


with open("pke") as f:
    lines = f.read().splitlines()
out = str(out) + ", " + str(lines[4]) + ", " + str(lines[6]) + ", " + str(lines[8]) + "\n"

file = open("../report/const.csv", "a")
file.write(out)
file.close()