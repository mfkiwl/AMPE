#!/usr/bin/env python
import sys
import subprocess

print("Test ChenShen98...")

#prepare initial conditions file
subprocess.call(["python3", "../../utils/make_sphere.py", "-x", "256", "-y", "256", "-z", "1", "-r", "100", "--qlen", "1", "initial256.nc"])


mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
exe = sys.argv[4]
inp = sys.argv[5]

#run AMPE
command = "{} {} {}".format(mpicmd,exe,inp)
output = subprocess.check_output(command,shell=True)

#analyse AMPE standard output
lines=output.split(b'\n')

end_reached = False
for line in lines:
  num_matches = line.count(b'cycle')
  if num_matches:
    print(line)
    words=line.split()
    if eval(words[6])>3000.:
      end_reached = True
      dt=words[10]
      if (eval(dt)-5.68)>1.e-2:
        print("Wrong dt")
        sys.exit(1)
  num_matches = line.count(b'fraction')
  if num_matches:
    print(line)
    if end_reached:
      words=line.split()
      if abs(eval(words[6])-0.21)>1.e-2:
        print("Wrong solid fraction")
        sys.exit(1)

if end_reached:
  sys.exit(0)
else:
  sys.exit(1)
