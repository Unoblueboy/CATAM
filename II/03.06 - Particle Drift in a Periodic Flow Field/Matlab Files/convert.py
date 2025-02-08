import os
import subprocess
print(os. getcwd())
for filename in os.listdir():
    if filename[-3:] == "eps":
        arg = "--outfile=" + filename[:-4] + "-eps-converted-to.pdf"
        subprocess.run(["epstopdf",arg,filename])
