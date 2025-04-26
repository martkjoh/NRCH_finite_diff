import subprocess
from multiprocessing import Pool


folder = "article/animate/"
def make_anim(i):
    subprocess.run("python3" + " " + folder + "animate_"+str(i)+".py", shell=True)

with Pool(10) as pool:
    pool.map(make_anim, range(1, 6))
