import sys

with open("kem") as f:
    lines = f.read().splitlines()

out = sys.argv[1]
out = str(out) + ", " + str(lines[4]) + ", " + str(lines[7])

with open("pke") as f:
    lines = f.read().splitlines()
out = out + ", " + str(lines[4]) + ", " + str(lines[7]) + ", " + str(lines[10]) + "\n"

print(out)
#file = open("../report/timing.csv", "a")
#file.write(out)
#file.close()