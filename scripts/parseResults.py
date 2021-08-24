import sys

with open("temp") as f:
    lines = f.read().splitlines()

out = sys.argv[1]
out = str(out) + ", " + str(lines[4]) + ", " + str(lines[7]) + ", " + str(lines[10]) + "\n"

file = open("../report/timing.csv", "a")
file.write(out)
file.close()
