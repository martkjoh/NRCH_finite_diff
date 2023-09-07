import subprocess

folder = "article/animate/"
for i in range(1, 6):
    print(i)
    subprocess.run("python3" + " " + folder + "animate_"+str(i)+".py", shell=True)
